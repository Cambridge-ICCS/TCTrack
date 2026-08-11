"""Module providing a class to run a machine-learning tracking algorithm.

References
----------
- `Instructions for adding a new tracker.
  <https://tctrack.readthedocs.io/en/latest/developer/adding_algorithms.html>`__
"""

from collections import deque
from collections.abc import Sequence
from dataclasses import dataclass

import warnings

import cf
import h5py
import numpy as np
import torch
from cftime import date2num

from tctrack.core import (
    TCTrackerMetadata,
    Trajectory,
)
from tctrack.core.ml_tracker import TCMLParameters, TCMLTracker


def _point_variables(candidate: dict) -> dict:
    """Strip a candidate down to the numeric variables of a trajectory point.

    Candidates carry the timestep they were found at, but
    :meth:`~tctrack.core.Trajectory.add_point` takes the time separately and
    requires every remaining value to be numeric, so it is dropped here.
    """
    return {key: value for key, value in candidate.items() if key != "time"}


@dataclass(repr=False)
class MLParameters(TCMLParameters):
    """Dataclass containing values for parameters used by MLTracker.

    Extends :class:`~tctrack.core.ml_tracker.TCMLParameters` with the input
    file path and sets the default HuggingFace repository for this model.

    See Also
    --------
    MLTracker : The tracker class that uses these parameters.
    """

    input_file: str = ""
    """CF-NetCDF file containing the input variables for the model."""

    hf_repo_id: str = "surbhigoel456/cyclone-TC-ML"
    """HuggingFace Hub repository ID for the cyclone detection model."""

    pressure_variables: Sequence[str] = (
        "relative_humidity",
        "air_temperature",
        "eastward_wind",
        "northward_wind",
        "atmosphere_relative_vorticity",
    )
    """CF identities of the pressure-level input variables, in model channel order.

    Correspond to ``pressure_var_codes`` in the reference training pipeline
    (`cyclone-track-ml <https://github.com/robert-edwin-rouse/cyclone-track-ml>`__),
    but given as the actual CF ``standard_name`` each variable is converted to
    (via cfgrib) in ERA5 output, which differs from the CDS API request name
    used there for several variables (e.g. ``"temperature"`` in the reference
    pipeline is ``"air_temperature"`` here).
    """

    pressure_levels: Sequence[float] = (1000, 750, 500)
    """Pressure levels (hPa) to extract for each of :attr:`pressure_variables`."""

    sst_variable: str = "ncvar%sst"
    """CF identity of the sea surface temperature variable.

    ERA5 sea surface temperature has no ``standard_name`` (cfgrib leaves it
    ``"unknown"``), so it must be selected by NetCDF variable name instead.
    """

    t2m_variable: str = "ncvar%t2m"
    """CF identity of the 2-metre temperature variable, used to gap-fill SST over land."""

    normalisation_stats_path: str | None = None
    """Path to a ``.nc`` file with per-channel ``mean`` and ``range`` variables
    (each shape ``(channel,)``, ordered as in :meth:`MLTracker.preprocess`),
    used to normalise input as ``(x - mean) / range``.

    These must be the statistics computed over the model's *training* set;
    they cannot be recomputed from a single inference file. Matches the
    ``data/normalisation_parameters.nc`` file saved by ``compiler.py`` in the
    reference training pipeline. If ``None``, :meth:`MLTracker.preprocess`
    raises rather than silently passing unnormalised data to the model.
    """

    merge_distance_deg: float = 2.0
    """Distance (degrees lat/lon) within which candidates in the same timestep
    are merged, keeping only the highest-scoring one.

    The model often splits a single storm into several adjacent clusters -
    either side of a lifecycle-class boundary, or where confidence dips below
    :attr:`~TCMLParameters.threshold` mid-storm - which would otherwise each
    become a separate track. Mirrors ``merge_dist`` in TempestExtremes'
    ``DetectNodes``, which resolves the same problem for its own candidates.
    Set to ``0`` to disable merging.
    """

    stitch_max_distance_deg: float = 3.0
    """Maximum distance (degrees lat/lon) between a track's last point and a
    candidate in the next timestep for them to be linked as the same storm.

    Not part of the reference training pipeline, which only trains
    per-pixel classification and has no timestep-to-timestep linking step.
    The default assumes 6-hourly input (as used for training) and typical
    cyclone translation speeds; reduce for higher-frequency input.
    """

    stitch_max_gap: int = 1
    """Number of consecutive timesteps a track may go unmatched before it is
    closed, tolerating brief drops below :attr:`~TCMLParameters.threshold`.
    """

    stitch_min_length: int = 2
    """Minimum number of observations a track must have to be kept, filtering
    out single-timestep detections likely to be noise.
    """


class MLTracker(TCMLTracker):
    """Class used to run the machine-learning cyclone tracking algorithm.

    The model weights are loaded from a HuggingFace Hub repository on first
    use and cached locally by ``huggingface_hub``.

    Attributes
    ----------
    parameters : MLParameters
        Class containing the parameters for the tracking algorithm.
    model : torch.jit.ScriptModule
        The loaded TorchScript model in evaluation mode.

    Examples
    --------
    >>> params = MLParameters(input_file="era5_2020.nc")
    >>> tracker = MLTracker(params)
    >>> tracker.run_tracker("trajectories.nc")
    """

    @property
    def _hf_filename(self) -> str:
        """Filename of the model weights in the HuggingFace repository."""
        return "cyclone-detect-ml-scripted.pt"

    def __init__(self, parameters: MLParameters):
        """Construct the MLTracker and load the model.

        Parameters
        ----------
        parameters : MLParameters
            Class containing the parameters for the tracking algorithm.

        Raises
        ------
        OSError
            If ``parameters.model_path`` is given but the file does not exist.
        huggingface_hub.errors.RepositoryNotFoundError
            If the HuggingFace repository cannot be found or accessed.
        """
        self.parameters: MLParameters = parameters
        self._trajectories: list[Trajectory] = []
        self._scores: list[dict] = []
        # TC locations found by detect(), consumed by stitch().
        self._candidates: list[dict] = []
        # Grid/time coordinates of the input file, populated by preprocess().
        self._lats: np.ndarray = np.array([])
        self._lons: np.ndarray = np.array([])
        self._times: np.ndarray = np.array([])
        self._load_model(parameters)

    @property
    def _parameters(self) -> list[MLParameters]:
        """A list of the parameter objects accessible from the base class."""
        return [self.parameters]

    def read_trajectories(self) -> list[Trajectory]:
        """Parse tracker outputs into a list of :class:`tctrack.core.Trajectory`.

        Returns
        -------
        list[Trajectory]
            A list of :class:`tctrack.core.Trajectory` objects.
        """
        return self._trajectories

    def _set_metadata(self) -> None:
        """Set the time and variable metadata attributes.

        Raises
        ------
        ValueError
            If a variable with time coordinate is not found in the input file.
        """
        vars_with_time = cf.read(self.parameters.input_file).select_by_construct("time")  # type: ignore[operator]
        if len(vars_with_time) == 0:
            msg = "No variable with 'time' coordinate found in the input file."
            raise ValueError(msg)
        time_coord = vars_with_time[0].coordinate("time")
        self._time_metadata = {
            "calendar": time_coord.calendar,
            "units": time_coord.units,
            "start_time": time_coord.data[0],
            "end_time": time_coord.data[-1],
        }

        self._variable_metadata = {}
        self._variable_metadata["variable_name"] = TCTrackerMetadata(
            properties={
                "standard_name": "...",
                "long_name": "...",
                "units": "...",
            }
        )

    def _load_normalisation_stats(self) -> tuple[np.ndarray, np.ndarray]:
        """Load per-channel normalisation statistics.

        Reads the ``mean``/``range`` variables from the ``.nc`` file produced
        by ``compiler.py`` in the reference training pipeline
        (``data/normalisation_parameters.nc``). Uses ``h5py`` rather than
        ``cf.read`` - this file's variables carry inherited GRIB/cfgrib
        attributes (e.g. an auxiliary ``number`` coordinate) that make
        ``cfdm`` build them through a data-creation path requiring a lazy
        ``scipy`` import; in this environment that import fails once
        ``torch`` has already been loaded (a ``libstdc++`` version conflict
        between the two), which ``h5py`` - already a TCTrack dependency -
        avoids entirely by reading the plain arrays directly.

        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray]
            ``(mean, range)`` arrays, each of shape ``(channel,)``, ordered to
            match the channel stack built by :meth:`preprocess`.

        Raises
        ------
        ValueError
            If :attr:`parameters.normalisation_stats_path` is not set.
        """
        if self.parameters.normalisation_stats_path is None:
            msg = (
                "parameters.normalisation_stats_path must be set to normalise "
                "model input - see MLTracker.preprocess."
            )
            raise ValueError(msg)
        with h5py.File(self.parameters.normalisation_stats_path, "r") as stats:
            mean = stats["mean"][:]
            value_range = stats["range"][:]
        return mean, value_range

    def preprocess(self) -> torch.Tensor:
        """Build the model input tensor from the configured ERA5 variables.

        Reproduces the channel layout the model was trained on: each of
        parameters.pressure_variables at each of
        parameters.pressure_levels,
        followed by a land/sea mask, followed by sea surface
        temperature gap-filled with 2-metre temperature over land. Channels
        are then normalised using _load_normalisation_stats.

        Returns
        -------
        torch.Tensor
            Tensor of shape ``(channel, time, lat, lon)``, normalised.

        Raises
        ------
        ValueError
            If a configured variable cannot be found in the input file, or if
            :attr:`parameters.normalisation_stats_path` is not set.
        """
        fields = cf.read(self.parameters.input_file)

        pressure_channels = []
        for variable in self.parameters.pressure_variables:
            field = fields.select_field(variable)
            for level in self.parameters.pressure_levels:
                pressure_channels.append(field.subspace(Z=level).squeeze("Z").array)

        sst = fields.select_field(self.parameters.sst_variable).array
        t2m = fields.select_field(self.parameters.t2m_variable).array

        # Static land/sea mask: SST is undefined over land. Taken from the first
        # timestep and held fixed, matching the reference preprocessing.
        land_mask = np.broadcast_to(
            np.ma.getmaskarray(sst)[0].astype(np.float32), sst.shape
        )

        # Gap-fill SST over land with 2m air temperature, as in training.
        sst_filled = np.where(np.ma.getmaskarray(sst), t2m, np.ma.filled(sst, 0.0))

        channels = [*pressure_channels, land_mask, sst_filled]

        self._lats = fields[0].coordinate("latitude").array
        self._lons = fields[0].coordinate("longitude").array
        self._times = fields[0].coordinate("time").datetime_array

        data = np.stack(channels, axis=0)
        mean, value_range = self._load_normalisation_stats()
        data = (data - mean[:, None, None, None]) / value_range[:, None, None, None]

        return torch.from_numpy(data).float()

    @property
    def _channel_names(self) -> list[str]:
        """Names for the model's input channels, in the order they are stacked.

        Matches the loop in :meth:`preprocess`: each pressure variable at each
        pressure level (variable-major), then the land/sea mask, then sea
        surface temperature.
        """
        names = [
            f"{variable}_{int(level)}"
            for variable in self.parameters.pressure_variables
            for level in self.parameters.pressure_levels
        ]
        return [*names, "land_mask", "sea_surface_temperature"]

    def _colocated_variables(
        self,
        frame: np.ndarray,
        candidate: dict,
        mean: np.ndarray,
        value_range: np.ndarray,
    ) -> dict:
        """Sample the input variables at a candidate's location.

        Values are taken at the grid point nearest the candidate's reported
        position and converted back to physical units, undoing the
        ``(x - mean) / range`` scaling applied in :meth:`preprocess` so that
        what is reported is e.g. a temperature in Kelvin rather than a
        normalised number.

        Parameters
        ----------
        frame : numpy.ndarray
            This timestep's normalised input, shape ``(channel, lat, lon)``.
        candidate : dict
            A candidate with ``lat``/``lon`` keys.
        mean : numpy.ndarray
            Per-channel means used to normalise, shape ``(channel,)``.
        value_range : numpy.ndarray
            Per-channel ranges used to normalise, shape ``(channel,)``.

        Returns
        -------
        dict
            One entry per input channel, keyed by :attr:`_channel_names`.
        """
        i = int(np.argmin(np.abs(self._lats - candidate["lat"])))
        j = int(np.argmin(np.abs(self._lons - candidate["lon"])))
        physical = frame[:, i, j] * value_range + mean
        return dict(zip(self._channel_names, physical.tolist(), strict=True))

    def detect(self) -> None:
        """Run the U-Net over the grid and locate the tropical cyclones in it.

        For each timestep the model is run, its per-pixel class probabilities
        are thresholded, and the surviving pixels are grouped into discrete
        candidate storms (see :meth:`_cluster_candidates` and
        :meth:`_merge_nearby_candidates`). The result is a list of TC
        locations on :attr:`_candidates`, which :meth:`stitch` then links into
        trajectories - mirroring how the other TCTrack trackers split the work,
        where detection emits candidate points and stitching consumes them.

        The full per-pixel probabilities are also kept on :attr:`_scores`,
        which retains the per-class detail that thresholding discards (useful
        for evaluation, e.g. one-vs-rest ROC per lifecycle class).

        Raises
        ------
        ValueError
            If the model's output spatial shape doesn't match the input
            file's latitude/longitude grid, it would indicate a bug.
        """
        data = self.preprocess()  # (channel, time, lat, lon)
        _, n_time, _, _ = data.shape

        # Reused to convert the normalised input back to physical units when
        # recording the variables co-located with each detection.
        mean, value_range = self._load_normalisation_stats()

        self._scores = []
        self._candidates = []
        with torch.no_grad():
            for t in range(n_time):
                frame = data[:, t, :, :].unsqueeze(0)  # (1, channel, lat, lon)
                output = self.model(frame)  # (1, 5, lat', lon')
                probs = torch.softmax(output, dim=1)[0].numpy()  # (5, lat', lon')

                if probs.shape[1:] != (len(self._lats), len(self._lons)):
                    msg = (
                        f"Model output shape {probs.shape[1:]} does not match "
                        f"the input grid {(len(self._lats), len(self._lons))}."
                    )
                    raise ValueError(msg)

                for i, lat in enumerate(self._lats):
                    for j, lon in enumerate(self._lons):
                        self._scores.append(
                            {
                                "time": self._times[t],
                                "lat": float(lat),
                                "lon": float(lon),
                                "probs": probs[:, i, j].tolist(),
                            }
                        )

                # Reduce this timestep's pixels to discrete storm locations.
                class_idx = probs.argmax(axis=0)
                class_prob = probs.max(axis=0)
                is_storm = (class_idx != 0) & (class_prob >= self.parameters.threshold)

                candidates = self._cluster_candidates(is_storm, class_idx, class_prob)
                candidates = self._merge_nearby_candidates(candidates)

                frame_data = data[:, t, :, :].numpy()  # (channel, lat, lon)
                for candidate in candidates:
                    candidate["time"] = self._times[t]
                    candidate.update(
                        self._colocated_variables(
                            frame_data, candidate, mean, value_range
                        )
                    )
                self._candidates.extend(candidates)

    def _cluster_candidates(
        self,
        is_storm: np.ndarray,
        class_idx: np.ndarray,
        class_prob: np.ndarray,
    ) -> list[dict]:
        """Group adjacent same-class detected pixels into single point candidates.

        Parameters
        ----------
        is_storm : numpy.ndarray
            Boolean array, shape ``(lat, lon)``, True for detected pixels.
        class_idx : numpy.ndarray
            Winning class index per pixel, shape ``(lat, lon)``.
        class_prob : numpy.ndarray
            Winning class probability per pixel, shape ``(lat, lon)``.

        Returns
        -------
        list[dict]
            One dict per cluster with keys ``lat``, ``lon``, ``class_index``,
            and ``score`` - the probability-weighted centroid, class, and
            peak probability of each connected group of same-class pixels.
        """
        unvisited = set(zip(*np.nonzero(is_storm), strict=False))
        candidates = []

        while unvisited:
            start = unvisited.pop()
            component_class = class_idx[start]
            queue = deque([start])
            pixels = [start]

            while queue:
                y, x = queue.popleft()
                for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                    neighbour = (y + dy, x + dx)
                    if (
                        neighbour in unvisited
                        and class_idx[neighbour] == component_class
                    ):
                        unvisited.discard(neighbour)
                        queue.append(neighbour)
                        pixels.append(neighbour)

            ys, xs = zip(*pixels, strict=False)
            ys_arr, xs_arr = np.array(ys), np.array(xs)
            weights = class_prob[ys_arr, xs_arr]
            candidates.append(
                {
                    "lat": float(np.average(self._lats[ys_arr], weights=weights)),
                    "lon": float(np.average(self._lons[xs_arr], weights=weights)),
                    "class_index": float(component_class),
                    "score": float(weights.max()),
                }
            )

        return candidates

    def _merge_nearby_candidates(self, candidates: list[dict]) -> list[dict]:
        """Merge candidates that are too close together into the strongest one.

        Parameters
        ----------
        candidates : list[dict]
            This timestep's candidates, as built by :meth:`_cluster_candidates`.

        Returns
        -------
        list[dict]
            The retained candidates, strongest first.
        """
        if self.parameters.merge_distance_deg <= 0:
            return candidates

        kept: list[dict] = []
        for candidate in sorted(candidates, key=lambda c: c["score"], reverse=True):
            if all(
                np.hypot(
                    candidate["lat"] - other["lat"], candidate["lon"] - other["lon"]
                )
                >= self.parameters.merge_distance_deg
                for other in kept
            ):
                kept.append(candidate)
        return kept

    def _nearest_candidate(
        self,
        track: dict,
        candidates: list[dict],
        unmatched: set[int],
    ) -> int | None:
        """Find the closest unmatched candidate to a track's last point.

        Parameters
        ----------
        track : dict
            Active track state, with a ``last`` dict holding its most recent
            ``lat``/``lon``.
        candidates : list[dict]
            This timestep's candidate detections, as built by
            :meth:`_cluster_candidates`.
        unmatched : set[int]
            Indices into ``candidates`` not yet claimed this timestep.

        Returns
        -------
        int | None
            Index of the nearest candidate within
            :attr:`parameters.stitch_max_distance_deg`, or ``None`` if none
            qualify.
        """
        best_index, best_distance = None, self.parameters.stitch_max_distance_deg #best distance is closest distance found so far
        for i in unmatched:
            candidate = candidates[i]
            distance = np.hypot(
                candidate["lat"] - track["last"]["lat"],
                candidate["lon"] - track["last"]["lon"],
            )
            if distance < best_distance: #condition to check if the distance is less than the best distance found so far
                best_index, best_distance = i, distance
        return best_index
    
    def stitch(self) -> list[Trajectory]:
        """Link the TC locations found by :meth:`detect` into trajectories.

        Consumes the candidates on :attr:`_candidates`, working through the
        timesteps in order and, at each one, extending existing tracks with a
        nearby candidate, retiring tracks that have gone unmatched for too
        long, and starting new tracks from whatever is left over.

        The result is also stored on :attr:`_trajectories`, so that
        :meth:`read_trajectories` - and therefore :meth:`to_netcdf` - can
        reach it.

        Returns
        -------
        list[Trajectory]
            One :class:`~tctrack.core.Trajectory` per track with at least
            :attr:`parameters.stitch_min_length` observations.
        """
        # Group detect()'s flat list of locations by the timestep they belong
        # to, so each timestep's candidates can be looked up by index below.
        by_timestep: dict = {time: [] for time in self._times}
        for candidate in self._candidates:
            by_timestep[candidate["time"]].append(candidate)

        active_tracks: list[dict] = [] #list of dictionaries with trajectories, location of last storm point, missed timesteps
        finished: list[Trajectory] = [] #list of dictionaries with trajectories that have been completed and no longer active
        next_id = 0


        for time in self._times:
            candidates = by_timestep[time]

            #set of candidates that are not yet matched to a track, used to avoid double-counting them
            unmatched = set(range(len(candidates)))
            for track in active_tracks:
                match = self._nearest_candidate(track, candidates, unmatched) #calls nearest candidate check function.
                if match is None:
                    track["missed"] += 1
                    continue
                point = candidates[match]
                track["traj"].add_point(time, _point_variables(point))
                track["last"] = point
                track["missed"] = 0
                unmatched.discard(match)


            still_active = []
            for track in active_tracks:
                if track["missed"] > self.parameters.stitch_max_gap:
                    finished.append(track["traj"])
                else:
                    still_active.append(track)
            active_tracks = still_active


            for i in unmatched:
                point = candidates[i]
                trajectory = Trajectory(next_id, time)
                trajectory.add_point(time, _point_variables(point))
                active_tracks.append(
                    {"traj": trajectory, "last": point, "missed": 0}
                )
                next_id += 1


        finished.extend(track["traj"] for track in active_tracks)
        self._trajectories = [
            trajectory
            for trajectory in finished
            if trajectory.observations >= self.parameters.stitch_min_length
        ]
        return self._trajectories

    def detections_to_netcdf(self, output_file: str) -> None:
        """Write the detected TC locations to a CF-compliant NetCDF point file.

        The detections from :meth:`detect` are a set of independent points -
        each a storm location at one time, with the input variables sampled
        there - rather than connected tracks. They are therefore written with
        ``featureType = "point"`` over a single ``observation`` dimension,
        instead of the trajectory layout :meth:`to_netcdf` produces.

        Parameters
        ----------
        output_file : str
            Filename for the output NetCDF file. Placed in the working
            directory unless a full path is given.

        Warns
        -----
        UserWarning
            If :meth:`detect` found no candidates, in which case no file is
            written.

        References
        ----------
        `CF-Conventions - H.1. Point Data
        <https://cfconventions.org/Data/cf-conventions/cf-conventions-1.11/cf-conventions.html#point-data>`_
        """
        if not self._time_metadata or not self._variable_metadata:
            self.set_metadata()

        if not self._candidates:
            msg = (
                "There are no detections to write, so no output file will be "
                "created. Has detect() been run?"
            )
            warnings.warn(msg, category=UserWarning, stacklevel=2)
            return

        domain_axis = cf.DomainAxis(size=len(self._candidates))
        domain_axis.nc_set_dimension("observation")
        dim_obs = cf.DimensionCoordinate(
            data=cf.Data(range(len(self._candidates))),
            properties={
                "standard_name": "observation",
                "long_name": "detection index",
            },
        )

        time_coord = cf.AuxiliaryCoordinate(
            data=cf.Data(
                date2num(
                    [candidate["time"] for candidate in self._candidates],
                    units=self.time_metadata["units"],
                    calendar=self.time_metadata["calendar"],
                ).tolist()
            ),
            properties={
                "standard_name": "time",
                "long_name": "time of detection",
                "units": cf.Units(
                    self.time_metadata["units"], calendar=self.time_metadata["calendar"]
                ),
            },
        )
        lat_coord = cf.AuxiliaryCoordinate(
            data=cf.Data([candidate["lat"] for candidate in self._candidates]),
            properties={
                "standard_name": "latitude",
                "long_name": "latitude of detection",
                "units": "degrees_north",
            },
        )
        lon_coord = cf.AuxiliaryCoordinate(
            data=cf.Data([candidate["lon"] for candidate in self._candidates]),
            properties={
                "standard_name": "longitude",
                "long_name": "longitude of detection",
                "units": "degrees_east",
            },
        )

        fields = []
        for variable in self._candidates[0]:
            if variable in {"time", "lat", "lon"}:
                continue

            metadata = self.variable_metadata.get(variable, TCTrackerMetadata({}))
            metadata.properties["featureType"] = "point"
            field = cf.Field(properties=metadata.properties)
            if "standard_name" not in metadata.properties:
                field.nc_set_variable(variable)

            axis = field.set_construct(domain_axis)
            field.set_construct(dim_obs, axes=(axis,))
            field.set_construct(time_coord, axes=(axis,))
            field.set_construct(lat_coord, axes=(axis,))
            field.set_construct(lon_coord, axes=(axis,))
            field.set_data(
                cf.Data([candidate[variable] for candidate in self._candidates]),
                axes=(axis,),
            )

            if self.global_metadata:
                field.nc_set_global_attributes(self.global_metadata)

            fields.append(field)

        cf.write(fields, output_file)  # type: ignore[operator]

    def run_tracker(self, output_file: str) -> None:
        """Run the tracker to obtain the tropical cyclone trajectories.

        This runs each of the individual steps of the tracking algorithm.

        Finally, the output is saved as a CF-netCDF file by calling :meth:`to_netcdf`.
        If the tracking algorithm doesn't produce an intermediate output file then the
        trajectories should be passed to it as an argument.

        Arguments
        ---------
        output_file : str
            Filename to which the tropical cyclone track trajectories are saved.

        Examples
        --------
        To set the parameters, instantiate a :class:`MLTracker` instance and run the
        algorithm to generate trajectories:

        >>> params = MLParameters(...)
        >>> tracker = MLTracker(params)
        >>> tracker.run_tracker("trajectories.nc")
        """
        # Run the steps for the tracking algorithm, eg:
        # self.preprocess()
        # self.detect()
        # trajectories = self.stitch()

        # Output the trajectories as a CF-netcdf file.
        # self.to_netcdf(output_file, trajectories)
