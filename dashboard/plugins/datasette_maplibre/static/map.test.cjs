const assert = require("node:assert/strict");
const { readFileSync } = require("node:fs");
const { join } = require("node:path");
const { test } = require("node:test");
const vm = require("node:vm");

// Skip CDN imports and defer browser initialisation; exercise the actual script.
const source = readFileSync(join(__dirname, "map.js"), "utf8")
	.replace(/^import .*;\r?\n/gm, "");

function loadScript(config = {}, data = {}) {
	const context = vm.createContext({
		window: config,
		document: { readyState: "loading", addEventListener() {} },
		console: { log() {} },
		performance: { now: () => 0 },
		location: { pathname: "/tracks/trajectories", search: "" },
		fetch: async () => ({ ok: true, json: async () => data }),
	});
	vm.runInContext(source, context, { filename: "map.js" });
	return vm.runInContext("({ buildGeoJSON, buildGeoJSONFromPairs, loadGeoJSON, featureAtPoint })", context);
}

// Normalize VM-realm objects for strict comparisons with test expectations.
function plain(value) {
	return JSON.parse(JSON.stringify(value));
}

function line(coordinates, properties) {
	return { type: "Feature", geometry: { type: "LineString", coordinates }, properties };
}

function point(coordinates, properties) {
	return { type: "Feature", geometry: { type: "Point", coordinates }, properties };
}

test("missing coordinates and malformed JSON report clear errors", () => {
	const { buildGeoJSON, buildGeoJSONFromPairs } = loadScript();
	for (const columns of [[], ["latitude"], ["longitude"]]) {
		assert.throws(() => buildGeoJSON({ columns, rows: [] }), /No latitude\/longitude in dataset/);
	}
	assert.throws(() => buildGeoJSONFromPairs({ columns: [], rows: [] }), /No lat_lon in dataset/);

	const invalid = ["", "not JSON", "[[NaN,2]]", "[[1,2]"];
	for (const value of invalid) {
		assert.throws(() => buildGeoJSON({
			columns: ["LAT_LON", "latitude", "longitude"],
			rows: [["[]", 10, 20], [value, 30, 40]],
		}), { name: "SyntaxError" });
	}
});

test("hit testing prefers exact hits and expands only line queries by five pixels", () => {
	const { featureAtPoint } = loadScript();
	const position = { x: 100, y: 200 };
	const layers = ["points", "tracks"];
	const line_layers = ["tracks"];
	const exact = { properties: { name: "point" } };
	const nearby = { properties: { name: "line" } };
	const calls = [];
	let results;
	const map = {
		queryRenderedFeatures(geometry, options) {
			calls.push(plain({ geometry, options }));
			return results.shift();
		},
	};

	results = [[exact]];
	assert.equal(featureAtPoint(map, position, layers, line_layers), exact);
	assert.deepEqual(calls.splice(0), [{ geometry: position, options: { layers } }]);

	results = [[], [nearby]];
	assert.equal(featureAtPoint(map, position, layers, line_layers), nearby);
	assert.deepEqual(calls.splice(0), [
		{ geometry: position, options: { layers } },
		{ geometry: [[95, 195], [105, 205]], options: { layers: line_layers } },
	]);

	results = [[], []];
	assert.equal(featureAtPoint(map, position, layers, line_layers), undefined);
	assert.equal(calls.splice(0).length, 2);

	results = [[]];
	assert.equal(featureAtPoint(map, position, ["points"], []), undefined);
	assert.equal(calls.length, 1);
});

test("observation rows retain points, properties and longitude unwrapping", () => {
	const data = {
		columns: ["Latitude", "LONGITUDE", "track", "layer", "score", "extra"],
		rows: [
			[10, 170, "a", "north", 5, "first"],
			[11, -170, "a", "north", 6, "second"],
			[12, -175, "b", "south", 7, "third"],
		],
	};
	const points = [
		point([170, 10], { track: "a", layer: "north", score: 5, extra: "first" }),
		point([-170, 11], { track: "a", layer: "north", score: 6, extra: "second" }),
		point([-175, 12], { track: "b", layer: "south", score: 7, extra: "third" }),
	];
	assert.deepEqual(plain(loadScript().buildGeoJSON(data)), {
		type: "FeatureCollection", features: points, points: 3, lines: 0,
	});
});

test("pairs take precedence, preserve order and all other properties", () => {
	const config = {
		DATASETTE_MAPLIBRE_LAYER_COLUMN: "layer",
	};
	const data = {
		columns: ["track", "LaT_LoN", "latitude", "longitude", "layer", "score", "extra"],
		rows: [
			["same", "[[10,170],[11,-170],[12,-160]]", 80, 90, "north", 5, "first"],
			["same", "[[20,-170],[21,170],[22,160]]", 81, 91, "south", 6, "second"],
		],
	};
	const expected = {
		type: "FeatureCollection",
		features: [
			line([[170, 10], [190, 11], [200, 12]], {
				track: "same", latitude: 80, longitude: 90, layer: "north", score: 5, extra: "first",
			}),
			line([[-170, 20], [-190, 21], [-200, 22]], {
				track: "same", latitude: 81, longitude: 91, layer: "south", score: 6, extra: "second",
			}),
		],
		points: 0,
		lines: 2,
	};
	assert.deepEqual(plain(loadScript(config).buildGeoJSON(data)), expected);
	assert.deepEqual(plain(loadScript().buildGeoJSON(data)), expected);
});

test("empty and single-position trajectories emit no geometry; exact 180-degree steps stay unchanged", () => {
	const { buildGeoJSON } = loadScript();
	for (const rows of [[], [["[]"]], [["[[10,20]]"]], [["[]"], ["[[10,20]]"]]]) {
		assert.deepEqual(plain(buildGeoJSON({ columns: ["lat_lon"], rows })), {
			type: "FeatureCollection", features: [], points: 0, lines: 0,
		});
	}
	assert.deepEqual(plain(buildGeoJSON({
		columns: ["lat_lon"],
		rows: [["[]"], ["[[10,20]]"], ["[[1,0],[2,180],[3,0],[4,-180]]"]],
	})), {
		type: "FeatureCollection",
		features: [line([[0, 1], [180, 2], [0, 3], [-180, 4]], {})],
		points: 0,
		lines: 1,
	});
});

test("loading rejects non-drawable collections and retains layer discovery for trajectory lines", async () => {
	await assert.rejects(loadScript({}, { columns: ["lat_lon"], rows: [] }).loadGeoJSON(),
		/Query returned no rows/);

	const data = {
		columns: ["lat_lon", "layer", "score"],
		rows: [["[[10,20],[11,21]]", "north", 5], ["[[12,22],[13,23]]", "south", 6]],
	};
	const collection = plain(await loadScript({ DATASETTE_MAPLIBRE_LAYER_COLUMN: "layer" }, data).loadGeoJSON());
	assert.deepEqual(collection, {
		type: "FeatureCollection",
		features: [
			line([[20, 10], [21, 11]], { layer: "north", score: 5 }),
			line([[22, 12], [23, 13]], { layer: "south", score: 6 }),
		],
		points: 0,
		lines: 2,
		layers: ["north", "south"],
	});
	assert.deepEqual(plain(await loadScript({}, data).loadGeoJSON()).layers, ["Data"]);
});
