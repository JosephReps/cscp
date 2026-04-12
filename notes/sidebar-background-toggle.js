document.addEventListener("DOMContentLoaded", function () {
  const sidebar =
    document.querySelector("#quarto-sidebar .sidebar-menu-container") ||
    document.querySelector("#quarto-sidebar .sidebar-item-container") ||
    document.querySelector("#quarto-sidebar nav") ||
    document.querySelector("#quarto-sidebar");

  if (!sidebar) return;

  const backgroundTitles = new Set([
    "Point Process Foundations",
    "First-Order Structure",
    "Second-Order Structure",
    "Structural Assumptions",
    "Random Fields",
    "Cox Processes",
    "Log-Gaussian Cox Processes"
  ]);

  const candidates = Array.from(sidebar.querySelectorAll("*"));

  const backgroundNodes = [];
  const matchedTitles = new Set();

  for (const el of candidates) {
    const txt = (el.textContent || "").trim();

    for (const title of backgroundTitles) {
      if (!matchedTitles.has(title) && txt === title) {
        const container =
          el.closest("li") ||
          el.closest(".sidebar-item") ||
          el.closest(".chapter") ||
          el;

        if (!backgroundNodes.includes(container)) {
          backgroundNodes.push(container);
          matchedTitles.add(title);
        }
        break;
      }
    }
  }

  if (backgroundNodes.length === 0) return;

  backgroundNodes.sort((a, b) => {
    if (a === b) return 0;
    const pos = a.compareDocumentPosition(b);
    return pos & Node.DOCUMENT_POSITION_FOLLOWING ? -1 : 1;
  });

  const firstNode = backgroundNodes[0];
  const parent = firstNode.parentNode;
  if (!parent) return;

  const wrapper = document.createElement("li");
  wrapper.className = "background-group";

  const button = document.createElement("button");
  button.className = "background-toggle";
  button.type = "button";
  button.setAttribute("aria-expanded", "false");
  button.innerHTML = `
    <span class="background-toggle-label">Background Material</span>
    <span class="background-toggle-chevron" aria-hidden="true">▸</span>
  `;

  const panel = document.createElement("ul");
  panel.className = "background-panel";
  panel.hidden = true;

  wrapper.appendChild(button);
  wrapper.appendChild(panel);

  parent.insertBefore(wrapper, firstNode);

  backgroundNodes.forEach((node) => {
    panel.appendChild(node);
  });

  button.addEventListener("click", () => {
    const expanded = button.getAttribute("aria-expanded") === "true";
    button.setAttribute("aria-expanded", String(!expanded));
    panel.hidden = expanded;
    wrapper.classList.toggle("is-open", !expanded);
  });

  const currentLink =
    panel.querySelector(".active") ||
    panel.querySelector('[aria-current="page"]');

  if (currentLink) {
    button.setAttribute("aria-expanded", "true");
    panel.hidden = false;
    wrapper.classList.add("is-open");
  }
});
