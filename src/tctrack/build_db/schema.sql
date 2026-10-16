-- TCTrack Database Schema
-- Stores collections of cyclone trajectories.
--
-- Hierarchy:
--  collections              Named collections of trajectories.
--  └─ files                 Individual files from TCTrack.
--     └─ trajectories       Single cyclone tracks.
--        └─ observations    Attributes from each trajectory observation.

pragma foreign_keys = on;

-- Ocean Basins
create table oceans (
    id    text primary key,
    name  text not null unique
);

insert into oceans (id, name) values
    ('NAT', 'North Atlantic Ocean'),
    ('SAT', 'South Atlantic Ocean'),
    ('NPA', 'North Pacific Ocean'),
    ('SPA', 'South Pacific Ocean'),
    ('IND', 'Indian Ocean'),
    ('ARC', 'Arctic Ocean'),
    ('SOU', 'Southern Ocean');

-- Optional landfall gates
create table gates (
    id       text primary key not null,
    name     text not null,
    min_lat  real not null,
    max_lat  real not null,
    min_lon  real not null,
    max_lon  real not null
);

-- Collections
-- A named group of track files.
create table collections (
    id           integer primary key,
    name         text    not null unique,

    title        text,
    description  text,

    created      text    not null default current_timestamp
);

-- Files
-- An individual NetCDF file produced by a tracker run.
create table files (
    id                  integer primary key,
    collection_id       integer not null references collections(id) on delete cascade,

    filepath            text    not null unique,
    filename            text    not null,

    tctrack_version     text,
    tracker_name        text    not null,
    tracker_parameters  text    not null,

    trajectories        integer not null,
    observations        integer not null,
    time_units          text    not null,
    time_calendar       text    not null,

    created             text    not null default current_timestamp
);

create index files_collection_idx on files(collection_id);

-- Trajectories
-- A single cyclone trajectory.
create table trajectories (
    id         integer primary key,
    file_id    integer not null references files(id) on delete cascade,

    start_end  text    check (start_end in ('S', 'E', 'SE')),

    ocean_id   text    references oceans(id),
    landfall   integer check (landfall in (0, 1)),
    gate_id    text    references gates(id),
    lat_lon    text    not null default '[]'
);

create index traj_file_idx on trajectories(file_id);
create index traj_ocean_idx on trajectories(ocean_id);
create index traj_gate_idx on trajectories(gate_id);
create index traj_landfall_idx on trajectories(landfall);


-- Observations
-- Individual observations from a trajectory.
create table observations (
    trajectory_id                  integer not null references trajectories(id) on delete cascade,
    sequence                       integer not null,
    date                           text    not null default current_timestamp,

    latitude                       real not null,
    longitude                      real not null,

    air_pressure_at_sea_level      real,
    surface_altitude               real,
    wind_speed                     real,
    atmosphere_relative_vorticity  real,

    sea_surface_temp_anomaly       real,
    roni                           real,

    ocean_id                       text    references oceans(id),
    distance_to_coast_km           real,
    landfall                       integer check (landfall in (0, 1)),
    gate_id                        text    references gates(id),

    primary key (trajectory_id, sequence)
);

create index obs_air_pressure_idx on observations(air_pressure_at_sea_level);
create index obs_surface_altitude_idx on observations(surface_altitude);
create index obs_wind_speed_idx on observations(wind_speed);
create index obs_ssta_idx on observations(sea_surface_temp_anomaly);
create index obs_roni_idx on observations(roni);
create index obs_ocean_idx on observations(ocean_id);
create index obs_distance_to_coast_km_idx on observations(distance_to_coast_km);
create index obs_landfall_idx on observations(landfall);
create index obs_landfall_distance_idx on observations(landfall, distance_to_coast_km);
create index obs_gate_idx on observations(gate_id);


-- Trajectory Statistics
-- Derived from stored observations after all files have been imported.
create table statistics (
    trajectory_id                 integer primary key references trajectories(id) on delete cascade,
    date_start                    text,
    date_end                      text,
    date_median                   text,
    duration_hours                real,
    observation_count             integer not null check (observation_count >= 0),
    wind_speed_min                real,
    wind_speed_max                real,
    wind_speed_median             real,
    wind_speed_mean               real,
    air_pressure_at_sea_level_min real,
    ace                           real,
    sea_surface_temp_anomaly_mean real,
    roni_mean                     real
);
create index stats_duration_hours_idx on statistics(duration_hours);
create index stats_date_start_idx on statistics(date_start);
create index stats_wind_speed_median_idx on statistics(wind_speed_median);
create index stats_wind_speed_max_idx on statistics(wind_speed_max);
create index stats_ace_idx on statistics(ace);


-- Views used for map display

-- Points - Observations without tracks
create view points as
select
	collections.name as collection,
	files.filename,
	trajectory_id,
	cast(substr(date, 1, 4) as integer) as year,
	sequence, date, latitude, longitude,
	air_pressure_at_sea_level, surface_altitude, wind_speed,
	sea_surface_temp_anomaly, roni,
	case
		when roni >= 0.5 then 'El Niño'
		when roni <= -0.5 then 'La Niña'
		else 'Neutral'
	end as enso,
	ob.ocean_id, ob.landfall, ob.gate_id, cast(sequence = 0 as integer) as genesis
from
	observations ob
	join trajectories on trajectories.id = trajectory_id
	join files on files.id = file_id
	join collections on collections.id = collection_id;

-- Points layered by collection
create view points_layer_collection as
select
	collection as layer_collection,
	year, sequence, date, latitude, longitude,
	air_pressure_at_sea_level, surface_altitude, wind_speed,
	sea_surface_temp_anomaly, roni, enso,
	ocean_id, landfall, gate_id, genesis
from points;

-- Points layered by filename
create view points_layer_file as
select
	filename as layer_file,
	year, sequence, date, latitude, longitude,
	air_pressure_at_sea_level, surface_altitude, wind_speed,
	sea_surface_temp_anomaly, roni, enso,
	ocean_id, landfall, gate_id, genesis
from points;

-- Points layered by ENSO
create view points_layer_enso as
select
	enso as layer_enso,
	year, sequence, date, latitude, longitude,
	air_pressure_at_sea_level, surface_altitude, wind_speed,
	sea_surface_temp_anomaly, roni,
	ocean_id, landfall, gate_id, genesis
from points;


-- Tracks - Trajectories with statistical summary of observations
create view tracks as
select
	collections.name as collection,
	files.filename,
	trajectory_id,
	tr.ocean_id, tr.landfall, tr.gate_id,
	cast(substr(stats.date_start, 1, 4) as integer) as year,
	stats.date_start, stats.date_end, stats.date_median, stats.duration_hours, stats.observation_count,
	stats.wind_speed_min, stats.wind_speed_max, stats.wind_speed_median, stats.wind_speed_mean,
	stats.air_pressure_at_sea_level_min, stats.ace,
	stats.sea_surface_temp_anomaly_mean, stats.roni_mean,
	case
		when roni_mean >= 0.5 then 'El Niño'
		when roni_mean <= -0.5 then 'La Niña'
		else 'Neutral'
	end as enso,
	tr.start_end,
	tr.lat_lon
from
	trajectories tr
	join statistics stats on stats.trajectory_id = tr.id
	join files on files.id = file_id
	join collections on collections.id = collection_id;

-- Tracks layered by collection
create view tracks_layer_collection as
select
	collection as layer_collection,
	trajectory_id,
	ocean_id, landfall, gate_id,
	year, date_start, date_end, date_median, duration_hours, observation_count,
	wind_speed_min, wind_speed_max, wind_speed_median, wind_speed_mean,
	air_pressure_at_sea_level_min, ace,
	sea_surface_temp_anomaly_mean, roni_mean, enso,
	start_end,
	lat_lon
from tracks;

-- Tracks layered by filename
create view tracks_layer_file as
select
	filename as layer_file,
	trajectory_id,
	ocean_id, landfall, gate_id,
	year, date_start, date_end, date_median, duration_hours, observation_count,
	wind_speed_min, wind_speed_max, wind_speed_median, wind_speed_mean,
	air_pressure_at_sea_level_min, ace,
	sea_surface_temp_anomaly_mean, roni_mean, enso,
	start_end,
	lat_lon
from tracks;

-- Tracks layered by ENSO phase
create view tracks_layer_enso as
select
	enso as layer_enso,
	trajectory_id,
	ocean_id, landfall, gate_id,
	year, date_start, date_end, date_median, duration_hours, observation_count,
	wind_speed_min, wind_speed_max, wind_speed_median, wind_speed_mean,
	air_pressure_at_sea_level_min, ace,
	sea_surface_temp_anomaly_mean, roni_mean,
	start_end,
	lat_lon
from tracks;
