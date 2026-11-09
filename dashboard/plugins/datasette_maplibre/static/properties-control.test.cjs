const assert = require("node:assert/strict");
const { readFileSync } = require("node:fs");
const { join } = require("node:path");
const { test } = require("node:test");
const vm = require("node:vm");

function createElement(tagName) {
	return {
		tagName,
		children: [],
		attributes: {},
		textContent: "",
		setAttribute(name, value) { this.attributes[name] = value; },
		append(...children) { this.children.push(...children); },
		replaceChildren() { this.children = []; this.textContent = ""; },
		remove() { this.removed = true; },
	};
}

const source = readFileSync(join(__dirname, "properties-control.js"), "utf8")
	.replace("export class", "class");
const context = vm.createContext({ document: { createElement } });
vm.runInContext(source, context);
const PropertiesControl = vm.runInContext("PropertiesControl", context);

function propertyPairs(panel) {
	return Array.from(panel.children[1].children[0].children, node => node.textContent);
}

test("missing and empty properties display an explicit empty state", () => {
	const control = new PropertiesControl();
	const panel = control.onAdd();
	for (const properties of [undefined, null, {}]) {
		control.setProperties({ previous: "value" });
		control.setProperties(properties);
		assert.equal(panel.children[1].textContent, "(no properties)");
		assert.equal(panel.children[1].children.length, 0);
	}
});

test("control displays and replaces properties as text, then removes its panel", () => {
	const control = new PropertiesControl();
	const panel = control.onAdd();
	assert.equal(panel.className, "maplibregl-ctrl properties-control");
	assert.equal(panel.attributes["aria-label"], "Feature properties");
	assert.equal(panel.tagName, "details");
		assert.equal(panel.children[0].tagName, "summary");
		assert.equal(panel.open, false);
		assert.equal(panel.children[1].textContent, "");

	control.setProperties({ "<b>name</b>": "<script>alert(1)</script>", zero: 0, flag: false, empty: null });
	assert.deepEqual(propertyPairs(panel), [
		"<b>name</b>", "<script>alert(1)</script>", "zero", "0", "flag", "false", "empty", "null",
	]);
	assert.equal(panel.scrollTop, 0);
	assert.equal(panel.open, true);

	panel.open = false;
	control.setProperties({ next: "selection" });
	assert.equal(panel.open, true);
	assert.deepEqual(propertyPairs(panel), ["next", "selection"]);

	control.onRemove();
	assert.equal(panel.removed, true);
});
