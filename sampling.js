/* Copyright (c) 2026 Neil Yuanting Li. MIT license; see LICENSE. */
(function (root, factory) {
  "use strict";
  if (typeof module === "object" && module.exports) module.exports = factory(require("./physics.js"));
  else root.TransportSampling = factory(root.TransportPhysics);
})(typeof globalThis !== "undefined" ? globalThis : this, function (P) {
  "use strict";

  // Polynomials use ascending powers. Only degree <=10 is needed here.
  function normalized(coefficients) {
    const scale = Math.max.apply(null, coefficients.map(Math.abs));
    if (scale === 0) return [0];
    const result = coefficients.map(function (value) { return value / scale; });
    while (result.length > 1 && Math.abs(result[result.length - 1]) < 1e-14) result.pop();
    return result;
  }
  function evaluate(coefficients, E) {
    let value = 0;
    for (let i = coefficients.length - 1; i >= 0; i--) value = value * E + coefficients[i];
    return value;
  }
  function derivative(coefficients) {
    return coefficients.slice(1).map(function (value, i) { return (i + 1) * value; });
  }
  function sortedUnique(values, tolerance) {
    return values.sort(function (a, b) { return a - b; }).filter(function (value, i, all) {
      return i === 0 || value - all[i - 1] > tolerance;
    });
  }

  // Derivative roots partition a polynomial into monotonic intervals. A sign
  // change then brackets each simple root; a zero at a derivative root also
  // identifies an even-multiplicity root. No energy mesh locates resonances.
  function rootsInBand(coefficients) {
    const polynomial = normalized(coefficients);
    if (polynomial.length === 1) return [];
    if (polynomial.length === 2) {
      const root = -polynomial[0] / polynomial[1];
      return root > -2 && root < 2 ? [root] : [];
    }
    const critical = rootsInBand(derivative(polynomial));
    const boundaries = [-2].concat(critical, [2]);
    const roots = [];
    for (const E of critical) {
      if (Math.abs(evaluate(polynomial, E)) < 2e-12) roots.push(E);
    }
    for (let i = 1; i < boundaries.length; i++) {
      let lower = boundaries[i - 1], upper = boundaries[i];
      let lowerValue = evaluate(polynomial, lower);
      const upperValue = evaluate(polynomial, upper);
      if (lowerValue === 0 && lower > -2) roots.push(lower);
      if (upperValue === 0 && upper < 2) roots.push(upper);
      if (lowerValue * upperValue >= 0) continue;
      for (let step = 0; step < 60 && upper - lower > 2e-13; step++) {
        const middle = (lower + upper) / 2;
        const middleValue = evaluate(polynomial, middle);
        if (middleValue === 0) { lower = middle; upper = middle; break; }
        if ((middleValue > 0) === (lowerValue > 0)) {
          lower = middle;
          lowerValue = middleValue;
        } else upper = middle;
      }
      roots.push((lower + upper) / 2);
    }
    return sortedUnique(roots, 2e-12);
  }

  function reflectionPolynomial(sites, hopping, contact, onsite) {
    if (sites === 1) return [-onsite, 1 - contact * contact];
    const determinants = [[1], [-onsite, 1]];
    for (let n = 2; n <= sites; n++) {
      const next = Array(n + 1).fill(0);
      determinants[n - 1].forEach(function (value, i) {
        next[i] -= onsite * value;
        next[i + 1] += value;
      });
      determinants[n - 2].forEach(function (value, i) { next[i] -= hopping * hopping * value; });
      determinants.push(next);
    }
    const result = determinants[sites].slice();
    determinants[sites - 1].forEach(function (value, i) { result[i + 1] -= contact * contact * value; });
    determinants[sites - 2].forEach(function (value, i) { result[i] += Math.pow(contact, 4) * value; });
    return result;
  }

  function chainCurve(parameters) {
    const input = parameters || {};
    const p = Object.assign({sites: 5, hopping: 1, contact: 1, onsite: 0, eta: 0,
      leadHopping: 1, leadOnsite: 0}, input);
    // This sampler is for the lossless, identical-lead experiment. A different
    // model needs its own extrema analysis rather than an unverified reuse.
    if (p.eta !== 0 || p.leadHopping !== 1 || p.leadOnsite !== 0) {
      throw new RangeError("chainCurve requires eta=0 and identical leads with hopping=1, onsite=0");
    }
    if (!Number.isInteger(p.sites) || p.sites < 1 || p.sites > 9 ||
        !Number.isFinite(p.hopping) || p.hopping < 0.4 || p.hopping > 1.4 ||
        !Number.isFinite(p.contact) || p.contact < 0.15 || p.contact > 1.4 ||
        !Number.isFinite(p.onsite) || Math.abs(p.onsite) > 1.5) {
      throw new RangeError("chainCurve supports N=1..9, hopping=0.4..1.4, contact=0.15..1.4, and |onsite|<=1.5");
    }
    const Q = reflectionPolynomial(p.sites, p.hopping, p.contact, p.onsite);
    const stationary = Array(Q.length + 1).fill(0);
    Q.forEach(function (value, i) {
      stationary[i + 1] += (1 - i) * value;
      if (i > 0) stationary[i - 1] += 4 * i * value;
    });

    // With P_n the isolated-chain determinant, the real reflection numerator
    // is Q=P_N-c^2 E P_(N-1)+c^4 P_(N-2). Inside the lead band,
    // T=K(4-E^2)/(Q^2+K(4-E^2)), K=c^4 t^(2N-2).
    // Thus all interior extrema are roots of Q or (4-E^2)Q'+E Q.
    // These formulas choose nodes; every plotted T is computed by physics.js.
    const edgeOffset = 2e-12;
    const lower = -2 + edgeOffset, upper = 2 - edgeOffset;
    const seeds = Array.from({length: 97}, function (_, i) { return lower + (upper - lower) * i / 96; });
    seeds.push.apply(seeds, rootsInBand(Q));
    seeds.push.apply(seeds, rootsInBand(stationary));
    // Resolve possible threshold cusps independently of the smooth-band mesh.
    for (let power = 2; power <= 11; power++) {
      seeds.push(-2 + Math.pow(10, -power), 2 - Math.pow(10, -power));
    }
    if (Number.isFinite(p.probe) && p.probe > lower && p.probe < upper) seeds.push(p.probe);
    const energies = sortedUnique(seeds.filter(function (E) { return E >= lower && E <= upper; }), 1e-14);
    const cache = new Map();
    function transmission(E) {
      if (!cache.has(E)) cache.set(E, P.chain(E, p).T);
      return cache.get(E);
    }
    const inside = [[energies[0], transmission(energies[0])]];
    const tolerance = 0.0005;
    function refine(a, ta, b, tb, depth) {
      const middle = (a + b) / 2;
      const quarter = (a + middle) / 2, threeQuarter = (middle + b) / 2;
      const tq = transmission(quarter), tm = transmission(middle), tt = transmission(threeQuarter);
      const error = Math.max(Math.abs(tq - (3 * ta + tb) / 4),
        Math.abs(tm - (ta + tb) / 2), Math.abs(tt - (ta + 3 * tb) / 4));
      if (error <= tolerance || depth >= 40 || b - a < 1e-13) {
        inside.push([b, tb]);
      } else {
        refine(a, ta, middle, tm, depth + 1);
        refine(middle, tm, b, tb, depth + 1);
      }
    }
    for (let i = 1; i < energies.length; i++) {
      refine(energies[i - 1], transmission(energies[i - 1]), energies[i], transmission(energies[i]), 0);
    }
    // Null entries break the SVG path. At E=+/-2 the channel count is zero;
    // the one-sided transmission limit need not equal that endpoint value.
    return [[-2.4, 0], [-2, P.chain(-2, p).T], [-2, null]].concat(inside,
      [[2, null], [2, P.chain(2, p).T], [2.4, 0]]);
  }

  return Object.freeze({chainCurve: chainCurve});
});
