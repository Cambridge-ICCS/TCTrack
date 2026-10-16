import * as maplibregl from "https://unpkg.com/maplibre-gl@6.12.0/dist/maplibre-gl.mjs";
import { LayerControl } from "https://unpkg.com/maplibre-gl-layer-control@0.17.5/dist/index.mjs";

// Pick up config passed from Python __init__ layer
const TABLE_NAME = window.DATASETTE_MAPLIBRE_TABLE_NAME || "Data";
const BASEMAP = window.DATASETTE_MAPLIBRE_BASEMAP || { style: "https://demotiles.maplibre.org/style.json" };
const POINT = window.DATASETTE_MAPLIBRE_POINT || { radius: 4, opacity: 1 };
const LINE = window.DATASETTE_MAPLIBRE_LINE || { thickness: 2, opacity: 1 };
const HEATMAP = window.DATASETTE_MAPLIBRE_HEATMAP || null;
const LAYER_COLUMN = window.DATASETTE_MAPLIBRE_LAYER_COLUMN || null;
const PALETTE = window.DATASETTE_MAPLIBRE_PALETTE || { single: "#000000", layers: ["#ffffff"] };
const MAX_LAYERS = window.DATASETTE_MAPLIBRE_MAX_LAYERS || 100;

// Local storage key for persisting map view settings
const STORAGE_KEY = "datasette-maplibre-view";

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
		+ "_shape=arrays&_extra=columns&_size=max&_nocount=on&_nofacet=on&_nosuggest=on";

	const res = await fetch(url);
	if (!res.ok) throw new Error("Fetch failed: " + res.status);

	return await res.json();
}


/**
	Build a GeoJSON FeatureCollection from a Datasette JSON response (_shape=arrays).

	Expects columns: lat_lon JSON pairs, or latitude and longitude.

	@returns A GeoJSON FeatureCollection for all data.
	         Includes two root-level keys:
	           points = number of Point features added
	           lines  = number of LineString features added

*/
function buildGeoJSON(data) {
	if (data.columns.some((c) => c.toLowerCase() === "lat_lon"))
		return buildGeoJSONFromPairs(data);

	return buildGeoJSONFromRows(data);
}

/**
	Build trajectory lines from a JSON array of [latitude, longitude] pairs
*/
function buildGeoJSONFromPairs({ columns, rows }) {
	const pairs_idx = columns.findIndex((c) => c.toLowerCase() === "lat_lon");
	if (pairs_idx == -1)
		throw new Error("No lat_lon in dataset");

	// Pair property column names with their row indexes to aid attachment to line features below
	const prop_pairs = [];
	columns.forEach((name, i) => {
		if (i != pairs_idx) prop_pairs.push([name, i]);
	});

	const features = [];
	rows.forEach((row) => {
		const pairs = JSON.parse(row[pairs_idx]);
		if (pairs.length < 2) return;

		// Unwrap longitude when a track crosses the antimeridian (±180°).
		// For example, 179° followed by -179° becomes 179° followed by 181°.
		// Adding or subtracting 360° preserves the location but continues the
		// line to the nearby map edge instead of stretching across the whole map.
		let previous_lon = pairs[0][1];
		const coordinates = pairs.map(([latitude, lon]) => {
			let longitude = lon;
			if (longitude < previous_lon - 180) longitude += 360;
			else if (longitude > previous_lon + 180) longitude -= 360;
			previous_lon = longitude;

			return [longitude, latitude];
		});

		// Add row property values to this line feature
		const properties = {};
		for (let i = 0; i < prop_pairs.length; i++)
			properties[prop_pairs[i][0]] = row[prop_pairs[i][1]];

		features.push({
			type: "Feature",
			geometry: { type: "LineString", coordinates },
			properties,
		});
	});

	console.log(`Features: ${features.length}`);

	return { type: "FeatureCollection", features, points: 0, lines: features.length };
}

function buildGeoJSONFromRows({ columns, rows }) {
	// Locate latitude and longitude column indices
	const lat_idx = columns.findIndex((c) => c.toLowerCase() === "latitude");
	const lon_idx = columns.findIndex((c) => c.toLowerCase() === "longitude");
	if (lat_idx == -1 || lon_idx == -1)
		throw new Error("No latitude/longitude in dataset");

	// Pair property column names with their row indexes to aid attachment to features below
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

	console.log(`Features: ${features.length}`);

	return { type: "FeatureCollection", features, points: features.length, lines: 0 };
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
		if (collection.layers.length > MAX_LAYERS) {
			collection.layers = collection.layers.slice(0, MAX_LAYERS);
			console.log(`WARNING Layers truncated at max: ${MAX_LAYERS}`);
		}
	} else {
		// Set a single base layer when no dynamic layers are specified
		collection.layers = [TABLE_NAME];
	}

	console.log(`Data load: ${Math.round(performance.now() - t0)}ms`);
	console.log(`Layers: ${collection.layers}`);

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
	Load a JSON value from browser local storage.

	@returns The stored value or null when unset or unavailable.
*/
function loadStored(key) {
	try {
		return JSON.parse(localStorage.getItem(key));
	} catch {
		return null;
	}
}

/**
	Save a JSON value to browser local storage.
*/
function saveStored(key, value) {
	try {
		localStorage.setItem(key, JSON.stringify(value));
	} catch {
		// Quietly ignore
	}
}


/**
	Initialise map and load data source.
*/
function init() {

	const t0 = performance.now();

	// Add map container just before the data table
	const parent = document.querySelector("div.table-wrapper");
	if (!parent) {
		console.log("Unable to find `div.table-wrapper` in HTML template");
		return;
	}

	// Create map element early to minimise visual pop-in
	const container = document.createElement("div");
	container.id = "datasette-maplibre";
	parent.before(container);

	// Load data concurrently with map loading
	const data_promise = loadGeoJSON();

	const map = new maplibregl.Map({
		container: container,
		style: BASEMAP.style,
		center: [0, 0],
		zoom: 1,
		cooperativeGestures: true,  // ctrl+scroll for zoom to allow page scrolling
		dragRotate: false,
		attributionControl: false
	});
	map.addControl(new maplibregl.NavigationControl({ showCompass: true }));
	map.addControl(new maplibregl.GlobeControl(), 'top-right');
	map.addControl(new LayerControl({ basemapStyleUrl: BASEMAP.style }), 'top-right');

	// Restore the persisted view settings (globe and camera position)
	const view = loadStored(STORAGE_KEY);
	map.on("style.load", () => {
		if (view) map.setProjection({ type: view.globe ? "globe" : "mercator" });
		if (view?.camera) map.jumpTo(view.camera);

		// Persist view changes
		const saveView = () => {
			const center = map.getCenter();

			saveStored(STORAGE_KEY, {
				globe: map.getProjection()?.type === "globe",
				camera: {
					center: [center.lng, center.lat],
					zoom: map.getZoom(),
					bearing: map.getBearing(),
					pitch: map.getPitch(),
				},
			});
		};
		map.on("projectiontransition", saveView);
		map.on("moveend", saveView);

		// Map outline - added as soon as the basemap has loaded to avoid it popping in later
		if (BASEMAP.outline) {
			map.addLayer({
				id: BASEMAP.outline.name,
				type: "line",
				source: BASEMAP.outline.tile_source,
				"source-layer": BASEMAP.outline.tile_layer,
				...(BASEMAP.outline.tile_class ? {
					filter: ["==", ["get", "class"], BASEMAP.outline.tile_class]
				} : {}),
				paint: {
					"line-color": BASEMAP.outline.colour,
					"line-opacity": BASEMAP.outline.opacity,
					"line-width": 1,
				}
			});
		}
	});

	map.on("load", async () => {
		// Wait for data load then set as data source
		const geojson = await data_promise.catch((e) => {
			console.error("[datasette-maplibre] Unable to load data.", e);
			return null;
		});
		if (!geojson) return;

		const source_id = "datasette-geojson";
		map.addSource(source_id, { type: "geojson", data: geojson });

		// Link heatmap weight to property value if defined and column exists
		let heatmap_weight = null;
		if (HEATMAP?.weight_property?.column && HEATMAP.weight_property.column in geojson.features[0].properties)
			heatmap_weight = [
				"interpolate",
				["exponential", HEATMAP.weight_property.exponential],
				["coalesce", ["get", HEATMAP.weight_property.column], 0.0],  // default to zero when no value present
				HEATMAP.weight_property.min[0], HEATMAP.weight_property.min[1],
				HEATMAP.weight_property.max[0], HEATMAP.weight_property.max[1],
			];

		// Add heatmap if configured
		if (heatmap_weight)
			map.addLayer({
				id: HEATMAP.weight_property.column,
				type: "heatmap",
				source: source_id,
				maxzoom: HEATMAP.max_zoom,
				paint: {
					"heatmap-weight": heatmap_weight,
					"heatmap-color": [
						"interpolate",
						["linear"],
						["heatmap-density"],
						0.0, HEATMAP.palette.low + "00",  // add transparent alpha channel for zero values
						0.1, HEATMAP.palette.low,
						0.5, HEATMAP.palette.mid,
						1.0, HEATMAP.palette.high,
					],
					// Adjust the heatmap radius by zoom level
					"heatmap-radius": [
						"interpolate",
						["linear"],
						["zoom"],
						0, 1,
						HEATMAP.max_zoom, HEATMAP.influence_radius,
					],
					// Fade out heatmap as it reaches maximum zoom
					"heatmap-opacity": [
						"interpolate",
						["linear"],
						["zoom"],
						HEATMAP.max_zoom - 2, HEATMAP.opacity,
						HEATMAP.max_zoom, 0,
					],
				},
			});

		// Add new layers for each distinct layer value
		let layer_id;
		let layer_filter;
		let layers_added = [];

		// Use feature colour for single layers and the layer palette for multiple
		let colour_idx = 0;
		let colour = geojson.layers.length > 1 ? PALETTE.layers[colour_idx] : PALETTE.single;

		// Create opacity expressions for point and line layers
		let circle_opacity = POINT.opacity || 1.0;
		if (POINT.opacity_property && POINT.opacity_property.column in geojson.features[0].properties)
			circle_opacity = [
				"interpolate",
				["linear"],
				["get", POINT.opacity_property.column],
				POINT.opacity_property.min[0], POINT.opacity_property.min[1],
				POINT.opacity_property.max[0], POINT.opacity_property.max[1]
			];
		let line_opacity = LINE.opacity || 1.0;
		if (LINE.opacity_property && LINE.opacity_property.column in geojson.features[0].properties)
			line_opacity = [
				"interpolate",
				["linear"],
				["get", LINE.opacity_property.column],
				LINE.opacity_property.min[0], LINE.opacity_property.min[1],
				LINE.opacity_property.max[0], LINE.opacity_property.max[1]
			];

		// Create circle-radius expression
		let circle_radius = POINT.radius;
		// Link to property value if defined and column exists
		if (POINT.radius_property && POINT.radius_property.column in geojson.features[0].properties)
			circle_radius = [
				"interpolate",
				["exponential", POINT.radius_property.exponential],
				["get", POINT.radius_property.column],
				POINT.radius_property.min[0], POINT.radius_property.min[1],
				POINT.radius_property.max[0], POINT.radius_property.max[1]
			];

		for (const layer of geojson.layers) {

			// Set filter expression to restrict data to this layer only
			if (geojson.layers.length > 1)
				layer_filter = ["==", ["get", LAYER_COLUMN], layer];

			// Line layer - lines added first so that points can be rendered on top
			// Only added if lines were made
			if (geojson.lines > 0) {
				layer_id = String(layer) + " (lines)";
				map.addLayer({
					id: layer_id,
					source: "datasette-geojson",
					type: "line",
					paint: {
						"line-color": colour,
						"line-width": LINE.thickness,
						"line-opacity": line_opacity,
					},
					filter: [
						"all",
						["==", ["geometry-type"], "LineString"],
						...(layer_filter ? [layer_filter] : [])
					],
				});
				layers_added.push(layer_id);
			}

			// Points layer
			if (geojson.points > 0) {
				layer_id = String(layer);
				map.addLayer({
					id: layer_id,
					source: "datasette-geojson",
					type: "circle",
					paint: {
						"circle-radius": circle_radius,
						"circle-color": colour,
						"circle-opacity": circle_opacity,
						"circle-stroke-color": "#00000020",
						"circle-stroke-width": 1,
					},
					filter: [
						"all",
						["==", ["geometry-type"], "Point"],
						...(layer_filter ? [layer_filter] : []),
					],
				});
				layers_added.push(layer_id);
			}

			// Advance layer colour index, reusing once exhausted
			colour_idx = (colour_idx + 1) % PALETTE.layers.length;
			colour = PALETTE.layers[colour_idx];
		}

		// Move map outline to top of the stack so it appears above all other layers
		if (BASEMAP.outline) map.moveLayer(BASEMAP.outline.name);

		// Set on-click popups and pointer style for all added layers
		for (const layer of layers_added) {
			map.on("mouseenter", layer, () => {
				map.getCanvas().style.cursor = "pointer";
			});
			map.on("mouseleave", layer, () => {
				map.getCanvas().style.cursor = "";
			});

			map.on("click", layer, (ev) => {
				const feature = ev.features[0];
				const html = propertiesHtml(feature.properties);
				new maplibregl.Popup()
					.setLngLat(ev.lngLat)
					.setHTML(html || "(no properties)")
					.addTo(map);
			});
		}

		console.log(`Map load: ${Math.round(performance.now() - t0)}ms`);

		// Zoom to data bounds when no view saved or no points are visible in the view
		const data_bounds = calcGeoJSONBounds(geojson);
		if (!view || !map.getBounds().intersects(data_bounds))
			map.fitBounds(data_bounds, { padding: 40, maxZoom: 8 });

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
