-- TCTrack Database Schema
-- Stores collections of cyclone trajectories.
--
-- Hierarchy:
--  collections              Named collections of trajectories.
--  └─ files                 Individual files from TCTrack.
--     └─ trajectories       Single cyclone tracks.
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
-- A single cyclone trajectory.
create table trajectories (
    id              integer primary key,
    file_id         integer not null references files(id) on delete cascade,

    start_end       text    check (start_end in ('S', 'E', 'SE'))
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


-- Views used for map display

-- Points only (trajectory_id renamed so it is not grouped into tracks as per the metadata group_by setting)
create view points as
select
	files.filename,
	trajectory_id as track_id,
	cast(substr(date, 1, 4) as integer) as year,
	sequence, date, latitude, longitude,
	air_pressure_at_sea_level, surface_altitude, wind_speed,
	cast(sequence = 0 as integer) as genesis
from
	observations ob
	join trajectories on trajectories.id = trajectory_id
	join files on files.id = file_id;

-- Points layered by year
create view points_layer_year as
select *, year as layer_year from points;

-- Points layered by filename
create view points_layer_file as
select *, filename as layer_file from points;


-- Tracks
create view tracks as
select
	files.filename,
	trajectory_id,
	tr.start_end,
	cast(substr(date, 1, 4) as integer) as year,
	sequence, date, latitude, longitude,
	air_pressure_at_sea_level, surface_altitude, wind_speed
from
	observations ob
	join trajectories tr on tr.id = ob.trajectory_id
	join files on files.id = file_id
order by
	trajectory_id, sequence;

-- Tracks layered by year
create view tracks_layer_year as
select *, year as layer_year from tracks;

-- Tracks layered by filename
create view tracks_layer_file as
select *, filename as layer_file from tracks;
