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
create view Points as
select
	collections.name as Collection,
	trajectory_id as Track,
	cast(substr(date, 1, 4) as integer) as Year,
	sequence as Sequence, date as Date, latitude as Latitude, longitude as Longitude,
	air_pressure_at_sea_level as SLP, surface_altitude as 'Surface Altitude', wind_speed as 'Wind Speed',
	sea_surface_temp_anomaly as SSTA, roni as RONI,
	case
		when roni >= 0.5 then 'El Niño'
		when roni <= -0.5 then 'La Niña'
		else 'Neutral'
	end as ENSO,
	ob.ocean_id as Ocean, ob.landfall as Landfall, ob.gate_id as Gate, cast(sequence = 0 as integer) as Genesis,
	files.filename as File
from
	observations ob
	join trajectories on trajectories.id = trajectory_id
	join files on files.id = file_id
	join collections on collections.id = collection_id;

-- Points layered by collection
create view 'Points - Collection Layer' as
select
	Collection as 'Layer Collection',
	Track, Year, Sequence, Date, Latitude, Longitude,
	SLP, "Surface Altitude", "Wind Speed",
	SSTA, RONI, ENSO,
	Ocean, Landfall, Gate, Genesis,
	File
from points;

-- Points layered by filename
create view 'Points - File Layer' as
select
	File as 'Layer File',
	Track, Year, Sequence, Date, Latitude, Longitude,
	SLP, "Surface Altitude", "Wind Speed",
	SSTA, RONI, ENSO,
	Ocean, Landfall, Gate, Genesis
from points;

-- Points layered by ENSO
create view 'Points - ENSO Layer' as
select
	ENSO As 'Layer ENSO',
	Track, Year, Sequence, Date, Latitude, Longitude,
	SLP, "Surface Altitude", "Wind Speed",
	SSTA, RONI,
	Ocean, Landfall, Gate, Genesis,
	Collection, File
from points;


-- Tracks - Trajectories with statistical summary of observations
create view Tracks as
select
	collections.name as Collection,
	trajectory_id as Track,
	tr.ocean_id as Ocean, tr.landfall as Landfall, tr.gate_id as Gate,
	cast(substr(stats.date_start, 1, 4) as integer) as Year,
	stats.date_start as 'Start', stats.date_end as 'End', stats.date_median as 'Median', stats.duration_hours as Hours, stats.observation_count as Points,
	stats.wind_speed_min as 'WS Min', stats.wind_speed_max as 'WS Max', stats.wind_speed_median as 'WS Median', stats.wind_speed_mean as 'WS Mean',
	stats.air_pressure_at_sea_level_min as 'SLP Min', stats.ace as ACE,
	stats.sea_surface_temp_anomaly_mean as 'SSTA Mean', stats.roni_mean as 'RONI Mean',
	case
		when roni_mean >= 0.5 then 'El Niño'
		when roni_mean <= -0.5 then 'La Niña'
		else 'Neutral'
	end as ENSO,
	tr.start_end as 'Start/End',
	files.filename as File,
	tr.lat_lon
from
	trajectories tr
	join statistics stats on stats.trajectory_id = tr.id
	join files on files.id = file_id
	join collections on collections.id = collection_id;

-- Tracks layered by collection
create view 'Tracks - Collection Layer' as
select
	collection as 'Layer Collection',
	Track,
	Ocean, Landfall, Gate,
	Year, "Start", "End", "Median", Hours, Points,
	"WS Min", "WS Max", "WS Median", "WS Mean",
	"SLP Min", ACE,
	"SSTA Mean", "RONI Mean", ENSO,
	"Start/End", File,
	lat_lon
from Tracks;

-- Tracks layered by filename
create view 'Tracks - File Layer' as
select
	File as 'Layer File',
	Track,
	Ocean, Landfall, Gate,
	Year, "Start", "End", "Median", Hours, Points,
	"WS Min", "WS Max", "WS Median", "WS Mean",
	"SLP Min", ACE,
	"SSTA Mean", "RONI Mean", ENSO,
	"Start/End",
	lat_lon
from Tracks;

-- Tracks layered by ENSO phase
create view 'Tracks - ENSO Layer' as
select
	ENSO as 'Layer ENSO',
	Track,
	Ocean, Landfall, Gate,
	Year, "Start", "End", "Median", Hours, Points,
	"WS Min", "WS Max", "WS Median", "WS Mean",
	"SLP Min", ACE,
	"SSTA Mean", "RONI Mean",
	"Start/End", File,
	lat_lon
from Tracks;
