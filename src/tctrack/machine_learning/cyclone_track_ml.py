"""Module providing a class to run a machine-learning tracking algorithm.

References
----------
- `Instructions for adding a new tracker.
  <https://tctrack.readthedocs.io/en/latest/developer/adding_algorithms.html>`__
"""

from collections.abc import Sequence
from dataclasses import dataclass

import cf
import numpy as np
import torch

from tctrack.core import (
    TCTrackerMetadata,
    Trajectory,
)
from tctrack.core.ml_tracker import TCMLParameters, TCMLTracker


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
        "temperature",
        "u_component_of_wind",
        "v_component_of_wind",
        "vorticity",
    )
    """CF identities of the pressure-level input variables, in model channel order.

    Matches ``pressure_var_codes`` in the reference training pipeline
    (`cyclone-track-ml <https://github.com/robert-edwin-rouse/cyclone-track-ml>`__).
    """

    pressure_levels: Sequence[float] = (1000, 750, 500)
    """Pressure levels (hPa) to extract for each of :attr:`pressure_variables`."""

    sst_variable: str = "sea_surface_temperature"
    """CF identity of the sea surface temperature variable."""

    t2m_variable: str = "2m_temperature"
    """CF identity of the 2-metre temperature variable, used to gap-fill SST over land."""

    normalisation_stats_path: str | None = None
    """Path to an ``.npz`` file with per-channel ``mean`` and ``range`` arrays
    (each shape ``(channel,)``, ordered as in :meth:`MLTracker.preprocess`),
    used to normalise input as ``(x - mean) / range``.

    These must be the statistics computed over the model's *training* set
    (see ``compiler.py`` in the reference training pipeline); they cannot be
    recomputed from a single inference file. If ``None``, :meth:`MLTracker.preprocess`
    raises rather than silently passing unnormalised data to the model.
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

        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray]
            ``(mean, range)`` arrays, each of shape ``(channel,)``, ordered to
            match the channel stack built by :meth:`preprocess`.

        Raises
        ------
        ValueError
            If :attr:`parameters.normalisation_stats_path` is not set.

        Notes
        -----
        Placeholder for loading the trained model's actual
        training-set statistics.
        """
        if self.parameters.normalisation_stats_path is None:
            msg = (
                "parameters.normalisation_stats_path must be set to normalise "
                "model input - see MLTracker.preprocess."
            )
            raise ValueError(msg)
        stats = np.load(self.parameters.normalisation_stats_path)
        return stats["mean"], stats["range"]

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
                pressure_channels.append(field.subspace(Z=level).squeeze().array)

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

    def detect(self) -> None:
        """Run the U-Net over the grid and save per-pixel class probabilities.

        Raises
        ------
        ValueError
            If the model's output spatial shape doesn't match the input
            file's latitude/longitude grid, it would indicate a bug.
        """
        data = self.preprocess()  # (channel, time, lat, lon)
        _, n_time, _, _ = data.shape

        out_lats, out_lons = self._lats, self._lons

        self._scores = []
        with torch.no_grad():
            for t in range(n_time):
                frame = data[:, t, :, :].unsqueeze(0)  # (1, channel, lat, lon)
                output = self.model(frame)  # (1, 5, lat', lon')
                probs = torch.softmax(output, dim=1)[0].numpy()  # (5, lat', lon')

                if probs.shape[1:] != (len(out_lats), len(out_lons)):
                    msg = (
                        f"Model output shape {probs.shape[1:]} does not match "
                        f"the input grid {(len(out_lats), len(out_lons))}."
                    )
                    raise ValueError(msg)

                for i, lat in enumerate(out_lats):
                    for j, lon in enumerate(out_lons):
                        self._scores.append(
                            {
                                "time": self._times[t],
                                "lat": float(lat),
                                "lon": float(lon),
                                "probs": probs[:, i, j].tolist(),
                            }
                        )

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
