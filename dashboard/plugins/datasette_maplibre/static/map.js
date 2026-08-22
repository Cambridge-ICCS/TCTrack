import * as maplibregl from "https://unpkg.com/maplibre-gl@^6.4.1/dist/maplibre-gl.mjs";

// Pick up config passed from Python __init__ layer
const GEOJSON_COLUMN_NAMES = window.DATASETTE_MAPLIBRE_GEOJSON_COLUMNS || [];
const BASEMAP_STYLE = window.DATASETTE_MAPLIBRE_STYLE || "https://demotiles.maplibre.org/style.json";

/**
	Fetch the current Datasette query as JSON.

	GeoJSON data cannot be taken from the Datasette HTML table because fields
	are truncated there. The row list is also paginated.
	This function retrieves all rows of the current query as JSON.

	A previous version used the query paramter _json=COL to load GeoJSON fields
	as actual JSON. However, the Datasette JSON parsing (load and dump) for every
	field and row was a performance penalty (2 seconds on a dataset of 5000 rows).
	Instead, GeoJSON fields are loaded as strings to avoid JSON processing on the
	server-side. JSON conversion is now performed during buildFeatureCollection.

	The Datasette setting, `max_returned_rows` is a hard limit on the URL
	_size parameter. It should be set to accommodate the total number
	of trajectories in the dataset.
*/
async function fetchRows() {

	// Get the current dataset as an array of pure JSON row objects.
	// Only the specified GeoJSON columns are retrieved (_col=).
	const url = location.pathname + ".json" + location.search
		+ (location.search ? "&" : "?") + "_size=max&_shape=array"
		+ GEOJSON_COLUMN_NAMES.map((column) => "&_col=" + column).join("");

	const res = await fetch(url);
	if (!res.ok) throw new Error("Fetch failed: " + res.status);

	return await res.json();
}


/**
	Build a single FeatureCollection from GeoJSON fields across all rows.

	Row fields are passed as string-encoded JSON for performance.
	They are parsed into actual JSON during this build.

	@returns A GeoJSON FeatureCollection for all data.
*/
function buildFeatureCollection(rows, geojson_fieldnames) {
	const features = rows.flatMap((row) =>
		geojson_fieldnames.flatMap((field) => {
			const geojson = JSON.parse(row[field]);
			return geojson.type == "Feature" ? geojson : geojson.features;
		})
	);
	console.log(`Features: ${features.length}`);

	return { type: "FeatureCollection", features };
}

/**
	Load data as GeoJSON for the current Datasette query.

	@returns A GeoJSON FeatureCollection for all data.
*/
async function loadGeoJSONData() {
	let rows;
	try {
		rows = await fetchRows();
	} catch (e) {
		console.error("[datasette-maplibre] Unable to load data.", e);
		return;
	}
	console.log(`Rows loaded: ${rows.length}`);

	return buildFeatureCollection(rows, GEOJSON_COLUMN_NAMES);
}


/**
	Calculate the bounding box for a GeoJSON dataset.
*/
function calcGeoJSONBounds(geojson) {
	const bounds = new maplibregl.LngLatBounds();

	function extend(coords) {
		if (typeof coords[0] === "number") bounds.extend(coords);
		else coords.forEach(extend);
	}
	geojson.features.forEach((f) => f.geometry && extend(f.geometry.coordinates));

	return bounds;
}


/**
	Make an HTML element for Feature properties.
*/
function propertiesHtml(properties) {
	const entries = Object.entries(properties || {});
	if (entries.length == 0) return "";

	const items = entries
		.map(([k, v]) => `<dt>${k}</dt><dd>${String(v)}</dd>`)
		.join("");

	return `<dl class="properties">${items}</dl>`;
}


/**
	Initialise map and load data source.
*/
function init() {

	if (!GEOJSON_COLUMN_NAMES.length) {
		console.log("GeoJSON column names missing. No data to display.");
		return;
	}

	// Load data concurrently with map loading
	const t0 = performance.now();
	const dataPromise = loadGeoJSONData();

	// Add map container and element
	const parent = document.querySelector("section.content");
	if (!parent) return;

	const container = document.createElement("div");
	container.id = "datasette-maplibre";
	parent.prepend(container);

	const map = new maplibregl.Map({
		container: container,
		style: BASEMAP_STYLE,
		center: [0, 0],
		zoom: 1,
	});
	map.addControl(new maplibregl.NavigationControl({ showCompass: false }));
	map.addControl(new maplibregl.GlobeControl(), 'top-right');

	map.on("load", async () => {
		// Wait for data load then set as the datasource
		const geojson = await dataPromise;
		const sourceId = "datasette-geojson";
		map.addSource(sourceId, { type: "geojson", data: geojson });

		console.log(`Data load: ${Math.round(performance.now() - t0)}ms`);

		// Add layers for tracks and points
		// The single data source is filtered by geometry type
		map.addLayer({
			id: "geojson_lines",
			source: "datasette-geojson",
			type: "line",
			paint: {
				"line-color": "#fafafa",
				"line-width": 3,
			},
			"filter": ["==", "$type", "LineString"]
		});
		map.addLayer({
			id: "geojson_points",
			source: "datasette-geojson",
			type: "circle",
			paint: {
				"circle-radius": 4,
				"circle-color": "#d32f2f",
				"circle-stroke-color": "#fafafa",
				"circle-stroke-width": 1,
			},
			"filter": ["==", "$type", "Point"]
		});

		// Set on-click popups for both layer types
		for (const layer of ["geojson_lines", "geojson_points"]) {
			map.on("click", layer, (ev) => {
				const feature = ev.features[0];
				const html = propertiesHtml(feature.properties);
				new maplibregl.Popup()
					.setLngLat(ev.lngLat)
					.setHTML(html || "(no properties)")
					.addTo(map);
			});
		}

		// Zoom to data bounds
		const bounds = calcGeoJSONBounds(geojson);
		map.fitBounds(bounds, { padding: 40, maxZoom: 15 });

		// Expose the map API for debugging
		window.datasette_maplibre_map = map;
	});
}


// Trigger initialisation
if (document.readyState === "loading") {
	document.addEventListener("DOMContentLoaded", init);
} else {
	init();
}
