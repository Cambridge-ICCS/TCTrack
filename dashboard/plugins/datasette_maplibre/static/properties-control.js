export class PropertiesControl {
	onAdd() {
		this.container = document.createElement("details");
				this.container.open = false;
		this.container.className = "maplibregl-ctrl properties-control";
		this.container.setAttribute("role", "region");
		this.container.setAttribute("aria-label", "Feature properties");

		const heading = document.createElement("summary");
		heading.textContent = "Properties";

		this.content = document.createElement("div");
		this.content.setAttribute("aria-live", "polite");

		this.container.append(heading, this.content);

		return this.container;
	}

	onRemove() {
		this.container.remove();
	}

	setProperties(properties) {
		this.container.open = true;
		this.container.scrollTop = 0;
		this.content.replaceChildren();

		const entries = Object.entries(properties || {});
		if (entries.length === 0) {
			this.content.textContent = "(no properties)";
			return;
		}

		const list = document.createElement("dl");
		for (const [name, value] of entries) {
			const term = document.createElement("dt");
			term.textContent = name;

			const description = document.createElement("dd");
			description.textContent = String(value);
			list.append(term, description);
		}

		this.content.append(list);
	}
}
