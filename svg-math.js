/* Copyright (c) 2026 Neil Yuanting Li. MIT license; see LICENSE.
 * Load before vendor/mathjax/tex-svg.js; await TransportSvgMath.ready before
 * rendering. MathJax handles SVG labels only; KaTeX handles document math.
 */
(function (root) {
  "use strict";

  const namespace = "http://www.w3.org/2000/svg";
  const cache = new Map();
  const cacheLimit = 256;
  let available = false, resolveReady, rejectReady;
  const ready = new Promise(function (resolve, reject) { resolveReady = resolve; rejectReady = reject; });
  // A failed local script load must reject readiness rather than strand the UI.
  function loadError(event) {
    const source = event.target;
    if (source && source.tagName === "SCRIPT" && /\/mathjax\/tex-svg\.js(?:[?#]|$)/.test(source.src)) {
      rejectReady(new Error("The local MathJax SVG runtime could not be loaded."));
      root.removeEventListener("error", loadError, true);
    }
  }
  root.addEventListener("error", loadError, true);

  function finite(value, name) {
    if (typeof value !== "number" || !Number.isFinite(value)) throw new TypeError(name + " must be finite");
    return value;
  }

  function template(tex) {
    if (cache.has(tex)) {
      const hit = cache.get(tex);
      cache.delete(tex);
      cache.set(tex, hit);
      return hit;
    }
    const container = root.MathJax.tex2svg(tex, {display: false, em: 16, ex: 8, containerWidth: 10000});
    const svg = container.querySelector("svg");
    if (!svg || svg.querySelector('[data-mml-node="merror"]')) throw new Error("Unable to render SVG label: " + tex);
    // This helper promises self-contained vector geometry. Missing glyphs,
    // linked resources, foreign HTML, and font-dependent SVG text are errors.
    const allowed = new Set(["svg", "g", "path", "rect", "line", "polyline", "polygon", "circle", "ellipse"]);
    for (const element of [svg].concat(Array.from(svg.querySelectorAll("*")))) {
      if (!allowed.has(element.localName)) throw new Error("SVG label requires unsupported nonvector content: " + tex);
      for (const attribute of Array.from(element.attributes)) {
        if (/href$/i.test(attribute.name) || /^on/i.test(attribute.name) || /url\s*\(/i.test(attribute.value)) {
          throw new Error("SVG label contains a resource reference: " + tex);
        }
      }
    }
    const box = (svg.getAttribute("viewBox") || "").trim().split(/\s+/).map(Number);
    if (box.length !== 4 || box.some(function (value) { return !Number.isFinite(value); }) || box[2] < 0 || box[3] < 0) {
      throw new Error("Invalid SVG label dimensions: " + tex);
    }
    const result = {box: box, children: Array.from(svg.childNodes).map(function (node) { return node.cloneNode(true); })};
    cache.set(tex, result);
    if (cache.size > cacheLimit) cache.delete(cache.keys().next().value);
    return result;
  }

  /** Return an SVG <g> containing inline glyph paths. (x, baselineY) is the
   * text anchor in parent SVG units; fontSize is the em size in those units.
   * anchor: start|middle|end; rotate: degrees around that anchor and baseline.
   * The result carries its original TeX, an accessible title, and its own color.
   */
  function render(tex, x, baselineY, options) {
    if (!available) throw new Error("Await TransportSvgMath.ready before rendering SVG math.");
    if (typeof tex !== "string" || tex.length === 0 || tex.length > 4096) throw new TypeError("tex must contain 1 to 4096 characters");
    finite(x, "x"); finite(baselineY, "baselineY");
    const p = Object.assign({fontSize: 16, anchor: "start", color: "#111", rotate: 0}, options || {});
    finite(p.fontSize, "fontSize"); finite(p.rotate, "rotate");
    if (p.fontSize <= 0) throw new RangeError("fontSize must be positive");
    if (!["start", "middle", "end"].includes(p.anchor)) throw new RangeError("anchor must be start, middle, or end");
    if (typeof p.color !== "string" || !/^(?:#[\da-f]{3,4}|#[\da-f]{6}|#[\da-f]{8}|[a-z]+)$/i.test(p.color)
        || p.color.toLowerCase() === "currentcolor" || !root.CSS.supports("color", p.color)) {
      throw new TypeError("color must be an explicit CSS color name or hexadecimal value");
    }
    const source = template(tex);
    // MathJax's TeX SVG font uses 1000 internal units per em. Its children
    // already place the baseline at y=0; viewBox y is the top, not the baseline.
    const scale = p.fontSize / 1000;
    const width = source.box[2] * scale;
    const offset = (p.anchor === "end" ? -width : p.anchor === "middle" ? -width / 2 : 0) - source.box[0] * scale;
    const group = root.document.createElementNS(namespace, "g");
    group.setAttribute("transform", "translate(" + x + " " + baselineY + ") rotate(" + p.rotate + ")");
    group.setAttribute("color", p.color);
    group.setAttribute("fill", p.color);
    group.setAttribute("role", "img");
    group.setAttribute("aria-label", tex);
    group.setAttribute("data-tex", tex);
    const title = root.document.createElementNS(namespace, "title");
    title.textContent = tex;
    group.appendChild(title);
    const geometry = root.document.createElementNS(namespace, "g");
    geometry.setAttribute("transform", "translate(" + offset + " 0) scale(" + scale + ")");
    for (const node of source.children) geometry.appendChild(node.cloneNode(true));
    group.appendChild(geometry);
    return new root.XMLSerializer().serializeToString(group);
  }

  root.TransportSvgMath = Object.freeze({ready: ready, render: render});
  root.MathJax = {
    loader: {
      load: [],
      failed: function (error) { rejectReady(new Error("MathJax SVG startup failed: " + error.message)); }
    },
    // All three packages are inside the pinned bundle. In particular, neither
    // autoload nor require is enabled, so labels cannot fetch TeX extensions.
    tex: {packages: ["base", "ams", "newcommand"]},
    svg: {fontCache: "none", mtextInheritFont: false},
    options: {enableMenu: false, enableAssistiveMml: false},
    startup: {
      typeset: false,
      ready: function () {
        try {
          root.MathJax.startup.defaultReady();
          root.MathJax.startup.promise.then(function () {
            available = true;
            root.removeEventListener("error", loadError, true);
            resolveReady(root.TransportSvgMath);
          }, rejectReady);
        } catch (error) { rejectReady(error); }
      }
    }
  };
})(window);
