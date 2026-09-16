import * as maplibregl from "https://unpkg.com/maplibre-gl@6.8.0/dist/maplibre-gl.mjs";
import { LayerControl } from "https://unpkg.com/maplibre-gl-layer-control@0.17.4/dist/index.mjs";

// Pick up config passed from Python __init__ layer
const BASEMAP_STYLE = window.DATASETTE_MAPLIBRE_STYLE || "https://demotiles.maplibre.org/style.json";
const GROUP_BY = window.DATASETTE_MAPLIBRE_GROUP_BY || null;
const LAYER_COLUMN = window.DATASETTE_MAPLIBRE_LAYER_COLUMN || null;
const LAYER_PALETTE = window.DATASETTE_MAPLIBRE_LAYER_PALETTE || null;

/**
	Fetch the current Datasette query as JSON.

	Data cannot be taken from the Datasette HTML table because fields
	are truncated there. The row list is also paginated.

	The Datasette setting, `max_returned_rows` is a hard limit on the URL
	_size parameter. It should be set to accommodate the total number
	of observations in the dataset.
*/
async function fetchRows() {

	// Get the current dataset as row arrays (the most efficient and compact form)
	// Facet and suggestion features are turned off for performance.
	const url = location.pathname + ".json" + location.search
		+ (location.search ? "&" : "?")
		+ "_size=max&_shape=arrays&_nocount=on&_nofacet=on&_nosuggest=on";

	const res = await fetch(url);
	if (!res.ok) throw new Error("Fetch failed: " + res.status);

	return await res.json();
}


/**
	Build a GeoJSON FeatureCollection from a Datasette JSON response (_shape=arrays).

	Expects columns: latitude, longitude

	@returns A GeoJSON FeatureCollection for all data.
	         Includes a root-level key, groups, as the number of groups created from the dataset
	         if a group-by is defined and matched.
*/
function buildGeoJSON({ columns, rows }) {
	// Locate latitude and longitude column indices
	const lat_idx = columns.findIndex((c) => c.toLowerCase() === "latitude");
	const lon_idx = columns.findIndex((c) => c.toLowerCase() === "longitude");
	if (lat_idx == -1 || lon_idx == -1)
		throw new Error("No latitude/longitude in dataset");

	// Pair property column names with their row indexes
	const prop_pairs = [];
	columns.forEach((name, i) => {
		if (i != lat_idx && i != lon_idx) prop_pairs.push([name, i]);
	});

	// Build lat/lon Point features with all other fields as properties
	const features = rows.map((row) => {
		const properties = {};
		for (let i = 0; i < prop_pairs.length; i++)
			properties[prop_pairs[i][0]] = row[prop_pairs[i][1]];

		return {
			type: "Feature",
			geometry: { type: "Point", coordinates: [row[lon_idx], row[lat_idx]] },
			properties,
		};
	});

	// Group-by processing
	const groups = new Map();
	const group_idx = GROUP_BY?.column ? columns.indexOf(GROUP_BY.column) : -1;

	if (group_idx != -1) {
		// Pair the group property columns with their indexes
		// Group-by and layer columns are always included as properties
		const group_prop_pairs = prop_pairs.filter(([name]) =>
			name == GROUP_BY.column ||
			name == LAYER_COLUMN ||
			GROUP_BY.properties?.includes(name)
		);

		// Build group map
		// Iterate rows to find groups and accumulate coordinates
		// Points are chained together in row order so group rows must be contiguous
		for (const row of rows) {
			const key = row[group_idx];

			// Get group for this row or create a new one
			let group = groups.get(key);
			if (!group) {
				// Get properties for this new group
				// Properties are considered to be the same for all points in a group so the first row values are taken
				const properties = {};
				for (let p = 0; p < group_prop_pairs.length; p++)
					properties[group_prop_pairs[p][0]] = row[group_prop_pairs[p][1]];

				// Create new group
				group = { coordinates: [], properties };
				groups.set(key, group);
			}

			group.coordinates.push([row[lon_idx], row[lat_idx]]);
		}

		// Create a LineString for each group entry
		groups.forEach((group) => {
			// A LineString needs at least two positions
			if (group.coordinates.length < 2) return;

			features.push({
				type: "Feature",
				geometry: { type: "LineString", coordinates: group.coordinates },
				properties: group.properties,
			});
		});
	}

	console.log(`Features: ${features.length}`);
	console.log(`Groups: ${groups.size}`);

	return { type: "FeatureCollection", features, groups: groups.size };
}

/**
	Load data as GeoJSON for the current Datasette query.

	A top-level "layers" key is added to expose distinct layer values.

	@returns A GeoJSON FeatureCollection for all data.
*/
async function loadGeoJSON() {
	const t0 = performance.now();

	const data = await fetchRows();
	console.log(`Rows loaded: ${data.rows.length}`);

	if (data.rows.length == 0)
		throw new Error("Query returned no rows");

	const collection = buildGeoJSON(data);

	// Collect distinct values from the layers column if set
	const layer_idx = data.columns.indexOf(LAYER_COLUMN);
	if (layer_idx != -1) {
		collection.layers = [...new Set(data.rows.map((row) => row[layer_idx]))];
	} else {
		// Set a single base layer when no dynamic layers are specified
		collection.layers = ["geojson"];
	}

	console.log(`Data load: ${Math.round(performance.now() - t0)}ms`);

	return collection;
}


/**
	Calculate the bounding box for a GeoJSON dataset.
*/
function calcGeoJSONBounds(geojson) {
	const bounds = new maplibregl.LngLatBounds();

	if (!geojson.features?.length)
		return bounds.extend([0, 0, 360, 0]);

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

	const t0 = performance.now();

	// Add map container and element
	const parent = document.querySelector("section.content");
	if (!parent) {
		console.log("Unable to find section.content in HTML template");
		return;
	}

	// Load data concurrently with map loading
	const data_promise = loadGeoJSON();

	// Create map element
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
	map.addControl(new LayerControl({ panelWidth: 500 }), 'top-right');

	map.on("load", async () => {
		// Wait for data load then set as data source
		const geojson = await data_promise.catch((e) => {
			console.error("[datasette-maplibre] Unable to load data.", e);
			return null;
		});
		if (!geojson) return;

		const source_id = "datasette-geojson";
		map.addSource(source_id, { type: "geojson", data: geojson });

		// Add new layers for each distinct layer value
		let layer_id;
		let layer_filter;
		let layers_added = [];
		let colour_idx = 0;
		for (const layer of geojson.layers) {

			// Set filter expression to restrict data to this layer only
			if (geojson.layers.length > 1)
				layer_filter = ["==", ["get", LAYER_COLUMN], layer];

			// Line layer - lines added first so that points can be rendered on top
			layer_id = layer + "_lines";
			map.addLayer({
				id: layer_id,
				source: "datasette-geojson",
				type: "line",
				paint: {
					"line-color": LAYER_PALETTE[colour_idx],
					"line-width": 3,
				},
				filter: [
					"all",
					["==", ["geometry-type"], "LineString"],
					...(layer_filter ? [layer_filter] : [])
				],
			});
			layers_added.push(layer_id);

			// Points layer
			layer_id = layer + "_points";
			map.addLayer({
				id: layer_id,
				source: "datasette-geojson",
				type: "circle",
				paint: {
					"circle-radius": 4,
					"circle-color": LAYER_PALETTE[colour_idx],
					"circle-stroke-color": "#00000080",
					"circle-stroke-width": 1,
				},
				filter: [
					"all",
					["==", ["geometry-type"], "Point"],
					...(layer_filter ? [layer_filter] : [])
				],
			});
			layers_added.push(layer_id);

			// Advance layer colour index - wraps at the end of palette (not ideal)
			colour_idx = (colour_idx + 1) % LAYER_PALETTE.length;
		}

		// Set on-click popups for all added layers
		for (const layer of layers_added) {
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

		console.log(`Map load: ${Math.round(performance.now() - t0)}ms`);
		console.log("Layers: ", geojson.layers);

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
