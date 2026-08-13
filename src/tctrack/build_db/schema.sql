-- TCTrack Database Schema
-- Stores collections of cyclone trajectories.
--
-- Hierarchy:
--  collections              Named collections of trajectories.
--  └─ files                 Individual files from TCTrack.
--     └─ trajectories       Single cyclone tracks with GeoJSON geometry.
--        └─ observations    Attributes from each trajectory observation.

pragma foreign_keys = on;

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
-- A single cyclone trajectory stored as GeoJSON.
--
-- geojson_track  LineString for the full vector path.
-- geojson_points FeatureCollection with a Point for each observation.
create table trajectories (
    id              integer primary key,
    file_id         integer not null references files(id) on delete cascade,

    start_end       text    check (start_end in ('S', 'E', 'SE')),

    geojson_track   text    not null,
    geojson_points  text    not null
);

create index trajectories_file_idx on trajectories(file_id);


-- Observations
-- Individual observations from a trajectory.
create table observations (
    trajectory_id                  integer not null references trajectories(id) on delete cascade,
    sequence                       integer not null,
    date                           text not null default current_timestamp,

    latitude                       real not null,
    longitude                      real not null,

    air_pressure_at_sea_level      real,
    surface_altitude               real,
    wind_speed                     real,
    atmosphere_relative_vorticity  real,

    primary key (trajectory_id, sequence)
);

create index air_pressure_idx on observations(air_pressure_at_sea_level);
create index surface_altitude_idx on observations(surface_altitude);
create index wind_speed_idx on observations(wind_speed);


-- Observation view
create view observation_view as
select file_id, files.filename, trajectory_id,
	ob.sequence, ob.date, ob.latitude, ob.longitude,
	ob.air_pressure_at_sea_level, ob.surface_altitude, ob.wind_speed, ob.atmosphere_relative_vorticity,
	cast(ob.sequence = 0 as integer) as genesis
from observations ob
	join trajectories on trajectories.id = trajectory_id
	join files on files.id = file_id;


-- Trajectory view
create view trajectory_view as
select trajectories.id as trajectory_id, start_end,
	geojson_track as geometry,
	file_id, filename, filepath,
	tctrack_version, tracker_name
from trajectories
	join files on files.id = file_id;
