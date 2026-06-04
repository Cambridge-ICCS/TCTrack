-- TCTrack Database Schema
-- Stores collections of cyclone trajectories.
--
-- Hierarchy:
--  collections         - named collections of trajectories
--    files             - individual files from TCTrack
--      trajectories    - single cyclone tracks with GeoJSON geometry
--        observations  - attributes from each trajectory observation

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

    filename            text    not null unique,

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
-- The geojson column contains a LineString feature for the full vector path
-- plus Point features for each path point along with attributes.
create table trajectories (
    id       integer primary key,
    file_id  integer not null references files(id) on delete cascade,

    geojson  text    not null
);

create index trajectories_file_idx on trajectories(file_id);


-- Observations
-- An individual observation from a trajectory.
--
-- Data is the same as trajectories.geojson but as individual, searchable rows.
create table observations (
    trajectory_id              integer not null references trajectories(id) on delete cascade,
    date                       text not null,

    latitude                   real not null,
    longitude                  real not null,
    grid_i                     real not null,
    grid_j                     real not null,

    air_pressure_at_sea_level  real not null,
    surface_altitude           real not null,
    wind_speed                 real not null
);

create index observations_trajectory_idx on observations(trajectory_id);
create index air_pressure_idx on observations(air_pressure_at_sea_level);
create index surface_altitude_idx on observations(surface_altitude);
create index wind_speed_idx on observations(wind_speed);
