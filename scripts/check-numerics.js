#!/usr/bin/env node
'use strict';
// Copyright (c) 2026 Neil Yuanting Li. MIT license; see LICENSE.

// Independent reference calculations use scattering wave amplitudes and exact
// Lorentzian integrals. All energies are in eV; currents are per spin in microamperes.
const fs = require('node:fs');
const path = require('node:path');
const P = require('../physics.js');

// Exact SI e and h, with a second factor e converting an eV integration to J.
const currentMicroampPerEV = (1.602176634e-19) ** 2 / 6.62607015e-34 * 1e6;
const add = (a, b) => [a[0] + b[0], a[1] + b[1]];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1]];
const scale = (a, x) => [a[0] * x, a[1] * x];
const multiply = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const normSquared = a => a[0] * a[0] + a[1] * a[1];
const divide = (a, b) => scale(multiply(a, [b[0], -b[1]]), 1 / normSquared(b));

/** Solve the discrete Schrödinger equation by propagating an outgoing wave
 * backwards from the right lead, then resolving incident/reflected left waves.
 * No Green function, self-energy, or matrix inversion is used here. The two
 * leads are identical, so transmitted/incident probability currents have the
 * same velocity factor. Device and contact Hamiltonian elements are -hopping.
 */
function scatteringReference(E, { sites, hopping, contact, onsite = 0, leadHopping = 1, leadOnsite = 0 }) {
  if (!(Math.abs(E - leadOnsite) < 2 * leadHopping) || contact <= 0) {
    throw new RangeError('Reference needs an open lead channel and nonzero contacts.');
  }
  const k = Math.acos(-(E - leadOnsite) / (2 * leadHopping));
  const forward = [Math.cos(k), Math.sin(k)];
  const backward = [forward[0], -forward[1]];
  let right = [1, 0]; // \psi_{N+1}, fixed normalization of the outgoing wave
  let here = scale(backward, leadHopping / contact); // right lead equation
  for (let site = sites; site >= 1; site--) {
    const rightHopping = site === sites ? contact : hopping;
    const leftHopping = site === 1 ? contact : hopping;
    const left = scale(add(scale(here, E - onsite), scale(right, rightHopping)), -1 / leftHopping);
    right = here;
    here = left;
  }
  const psiZero = here;
  const psiMinusOne = scale(add(scale(psiZero, E - leadOnsite), scale(right, contact)), -1 / leadHopping);
  const incident = divide(sub(psiMinusOne, multiply(psiZero, forward)), sub(backward, forward));
  const reflected = sub(psiZero, incident);
  return { T: 1 / normSquared(incident), R: normSquared(reflected) / normSquared(incident) };
}

function lorentzianArea(epsilon, width, lower, upper) {
  const halfWidth = width / 2;
  return (Math.atan((upper - epsilon) / halfWidth) - Math.atan((lower - epsilon) / halfWidth)) / Math.PI;
}

/** Exact zero-temperature, finite-window WBL results, including an empty
 * absorbing channel of width 2\eta when \eta>0. Occupied ranges are clipped to the
 * actual (possibly shifted) integration window, not necessarily to either chemical potential.
 */
function zeroTemperatureReference(p) {
  const gamma = p.gammaL + p.gammaR;
  const width = gamma + 2 * (p.eta || 0);
  const a = width / 2;
  const lower = -p.limit + (p.offset || 0);
  const upper = p.limit + (p.offset || 0);
  const muL = Math.max(lower, Math.min(upper, p.bias / 2));
  const muR = Math.max(lower, Math.min(upper, -p.bias / 2));
  if (!(a > 0 && upper > lower)) {
    throw new RangeError('Reference requires positive spectral width and integration-window width.');
  }
  const primitive = x => Math.atan((x - p.epsilon) / a) / a;
  const filledL = primitive(muL) - primitive(lower);
  const filledR = primitive(muR) - primitive(lower);
  const pair = p.gammaL * p.gammaR * (primitive(muL) - primitive(muR));
  const sinkL = 2 * (p.eta || 0) * p.gammaL * filledL;
  const sinkR = 2 * (p.eta || 0) * p.gammaR * filledR;
  return {
    I: currentMicroampPerEV * pair,
    IL: currentMicroampPerEV * (pair + sinkL),
    IR: currentMicroampPerEV * (-pair + sinkR),
    Isink: currentMicroampPerEV * (sinkL + sinkR),
    n: (p.gammaL * filledL + p.gammaR * filledR) / (2 * Math.PI),
    spectralArea: lorentzianArea(p.epsilon, width, lower, upper)
  };
}

function maximum(values) { return Math.max(...values.map(Math.abs)); }

// Independent quadrature for the finite-temperature interface checks. Its
// adaptive Simpson rule and analytic Lorentzian integrands share neither the
// Green-function evaluations nor the uniform trapezoid grid of physics.js.
function adaptiveSimpson(fn, lower, upper, tolerance) {
  function refine(a, fa, m, fm, b, fb, whole, budget, depth) {
    const lm = (a + m) / 2, rm = (m + b) / 2;
    const flm = fn(lm), frm = fn(rm);
    const left = (m - a) * (fa + 4 * flm + fm) / 6;
    const right = (b - m) * (fm + 4 * frm + fb) / 6;
    const correction = left + right - whole;
    if (Math.abs(correction) <= 15 * budget) return left + right + correction / 15;
    if (depth === 0) throw new RangeError('Independent Simpson reference did not converge.');
    return refine(a, fa, lm, flm, m, fm, left, budget / 2, depth - 1)
      + refine(m, fm, rm, frm, b, fb, right, budget / 2, depth - 1);
  }
  const middle = (lower + upper) / 2;
  const a = fn(lower), m = fn(middle), b = fn(upper);
  return refine(lower, a, middle, m, upper, b, (upper - lower) * (a + 4 * m + b) / 6, tolerance, 24);
}

function finiteTemperatureReference(p) {
  if (!(p.kT > 0)) throw new RangeError('Finite-temperature reference requires kT > 0.');
  const lower = -p.limit + (p.offset || 0), upper = p.limit + (p.offset || 0);
  const halfWidth = (p.gammaL + p.gammaR) / 2 + (p.eta || 0);
  const muL = p.bias / 2, muR = -p.bias / 2;
  const anchors = [lower, upper, p.epsilon - halfWidth, p.epsilon, p.epsilon + halfWidth];
  for (const mu of [muL, muR]) {
    for (const offset of [-8, -2, 0, 2, 8]) anchors.push(mu + offset * p.kT);
  }
  const nodes = [...new Set(anchors.filter(E => E >= lower && E <= upper))].sort((a, b) => a - b);
  const f = (E, mu) => (1 - Math.tanh((E - mu) / (2 * p.kT))) / 2;
  function integral(fn) {
    let result = 0;
    for (let i = 1; i < nodes.length; i++) result += adaptiveSimpson(fn, nodes[i - 1], nodes[i], 2e-12 / nodes.length);
    return result;
  }
  const denominator = E => (E - p.epsilon) ** 2 + halfWidth ** 2;
  return {
    I: currentMicroampPerEV * integral(E => p.gammaL * p.gammaR * (f(E, muL) - f(E, muR)) / denominator(E)),
    n: integral(E => (p.gammaL * f(E, muL) + p.gammaR * f(E, muR)) / denominator(E)) / (2 * Math.PI),
    spectralArea: lorentzianArea(p.epsilon, 2 * halfWidth, lower, upper)
  };
}

function interpolateCurve(curve, E) {
  let left = 0, right = curve.length - 1;
  if (E < curve[left][0] || E > curve[right][0]) throw new RangeError('Energy outside plotted curve.');
  while (right - left > 1) {
    const middle = (left + right) >> 1;
    if (curve[middle][0] <= E) left = middle;
    else right = middle;
  }
  const a = curve[left], b = curve[right];
  if (a[0] === E && Number.isFinite(a[1])) return a[1];
  if (b[0] === E && Number.isFinite(b[1])) return b[1];
  if (!Number.isFinite(a[1]) || !Number.isFinite(b[1])) throw new RangeError('Cannot interpolate across a plot path break.');
  return a[1] + (b[1] - a[1]) * (E - a[0]) / (b[0] - a[0]);
}

/** Independently scan the open lead band by wave recurrence. Refine detected
 * maxima on a finer uniform local mesh, locate their half-height crossings,
 * then check plotted interpolation both at peaks and along their flanks.
 * These are sampled error measurements, not a uniform mathematical bound.
 */
function checkChainCurve(curve, parameters, intervals = 80000) {
  const lower = -1.99999, upper = 1.99999;
  const step = (upper - lower) / intervals;
  let previous = null, beforePrevious = null;
  let maxError = 0, worstEnergy = lower, maxLibraryDifference = 0;
  const candidates = [];
  const reference = E => scatteringReference(E, parameters).T;
  function checkPoint(E, truth) {
    const error = Math.abs(interpolateCurve(curve, E) - truth);
    if (error > maxError) { maxError = error; worstEnergy = E; }
    return error;
  }
  for (let i = 0; i <= intervals; i++) {
    const E = lower + i * step;
    const T = reference(E);
    checkPoint(E, T);
    if (beforePrevious && previous[1] > beforePrevious[1] && previous[1] >= T && previous[1] > 0.5) candidates.push(previous[0]);
    beforePrevious = previous;
    previous = [E, T];
  }
  function halfHeight(energy, peakT, direction) {
    let inside = energy, distance = step;
    let outside = Math.max(lower, Math.min(upper, energy + direction * distance));
    while (reference(outside) > peakT / 2) {
      if (outside === lower || outside === upper) return outside;
      inside = outside;
      distance *= 1.5;
      outside = Math.max(lower, Math.min(upper, energy + direction * distance));
    }
    for (let i = 0; i < 40; i++) {
      const midpoint = (inside + outside) / 2;
      if (reference(midpoint) > peakT / 2) inside = midpoint;
      else outside = midpoint;
    }
    return (inside + outside) / 2;
  }
  const peaks = candidates.map(candidate => {
    let energy = candidate, transmission = reference(candidate);
    for (let i = 0; i <= 800; i++) {
      const E = Math.max(lower, Math.min(upper, candidate - step + 2 * step * i / 800));
      const T = reference(E);
      if (T > transmission) { energy = E; transmission = T; }
    }
    const leftHalf = halfHeight(energy, transmission, -1);
    const rightHalf = halfHeight(energy, transmission, 1);
    let localError = 0;
    for (const fraction of [-1.5, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.5]) {
      const scale = fraction < 0 ? energy - leftHalf : rightHalf - energy;
      const E = Math.max(lower, Math.min(upper, energy + fraction * scale));
      const truth = reference(E);
      localError = Math.max(localError, checkPoint(E, truth));
      maxLibraryDifference = Math.max(maxLibraryDifference, Math.abs(P.chain(E, { ...parameters, eta: 0 }).T - truth));
    }
    return { energy, transmission, halfMaximumWidth: rightHalf - leftHalf, plottedPeak: interpolateCurve(curve, energy), maximumLocalInterpolationError: localError };
  });
  return { parameters, plotSamples: curve.length, referenceScan: { lower, upper, intervals, step, localPeakStep: 2 * step / 800 },
    sampledMaximumInterpolationError: maxError, worstEnergy, maximumLibraryReferenceDifference: maxLibraryDifference, peaks };
}

function fanoReference(E, { sideEnergy, coupling }) {
  if (coupling === 0) return 0.25 / (E * E + 0.25);
  const detuning = E - sideEnergy;
  return 0.25 * detuning ** 2 / ((E * detuning - coupling ** 2) ** 2 + 0.25 * detuning ** 2);
}

function checkFanoCurve(curve, parameters, intervals = 60000) {
  const lower = -3, upper = 3, step = (upper - lower) / intervals;
  let maxError = 0, worstEnergy = lower;
  function check(E) {
    const error = Math.abs(interpolateCurve(curve, E) - fanoReference(E, parameters));
    if (error > maxError) { maxError = error; worstEnergy = E; }
  }
  for (let i = 0; i <= intervals; i++) check(lower + i * step);
  // Also inspect each rendered segment. A uniform reference grid alone can
  // skip narrow peaks and flanks at the minimum allowed nonzero coupling.
  for (let i = 1; i < curve.length; i++) {
    for (const fraction of [0.1, 0.2, 0.4, 0.6, 0.8, 0.9]) {
      check(curve[i - 1][0] + fraction * (curve[i][0] - curve[i - 1][0]));
    }
  }
  const { sideEnergy, coupling } = parameters;
  const split = Math.hypot(sideEnergy / 2, coupling);
  const peakEnergies = coupling === 0 ? [0] : [sideEnergy / 2 - split, sideEnergy / 2 + split];
  return { parameters, plotSamples: curve.length, referenceScan: { lower, upper, intervals, step },
    sampledMaximumInterpolationError: maxError, worstEnergy,
    atSideEnergy: interpolateCurve(curve, sideEnergy),
    peaks: peakEnergies.map(energy => ({ energy, exactTransmission: fanoReference(energy, parameters), plottedTransmission: interpolateCurve(curve, energy) })) };
}

function rounded(value) { return Number(value.toPrecision(12)); }
function roundTree(value) {
  if (typeof value === 'number') return rounded(value);
  if (Array.isArray(value)) return value.map(roundTree);
  if (value && typeof value === 'object') return Object.fromEntries(Object.entries(value).map(([key, item]) => [key, roundTree(item)]));
  return value;
}

function buildReport() {
  const S = require('../sampling.js');
  const openEnergies = [-1.8, -1.1, -0.3, 0, 0.47, 1.2, 1.8];
  const matched = [];
  const mismatch = [];
  for (const sites of [1, 2, 5, 9]) {
    const params = { sites, hopping: 1, contact: 1, onsite: 0, eta: 0 };
    const values = openEnergies.map(E => P.chain(E, params));
    matched.push({ sites, maximumTransmissionError: maximum(values.map(q => q.T - 1)), maximumConservationError: maximum(values.map(q => q.T + q.R - 1)) });
  }
  for (const params of [
    { sites: 3, hopping: 0.83, contact: 0.57, onsite: 0.21, leadHopping: 1, leadOnsite: 0 },
    { sites: 7, hopping: 1.17, contact: 0.76, onsite: -0.14, leadHopping: 1, leadOnsite: 0 }
  ]) {
    const residuals = openEnergies.map(E => {
      const green = P.chain(E, { ...params, eta: 0 });
      const wave = scatteringReference(E, params);
      return { transmission: green.T - wave.T, reflection: green.R - wave.R };
    });
    mismatch.push({ parameters: params, maximumTransmissionError: maximum(residuals.map(r => r.transmission)), maximumReflectionError: maximum(residuals.map(r => r.reflection)) });
  }

  const biasParameters = { epsilon: 0.13, gammaL: 0.18, gammaR: 0.32, bias: 0.4, kT: 0, eta: 0, limit: 10, intervals: 40000 };
  const bias = [0, 0.03].map(eta => {
    const parameters = { ...biasParameters, eta };
    const numeric = P.bias(parameters);
    const exact = zeroTemperatureReference(parameters);
    return { parameters,
      integration: { scheme: numeric.integrationScheme, nominalStep: numeric.gridStep, actualIntervals: numeric.actualIntervals, insertedBreakpoints: numeric.insertedBreakpoints, samples: numeric.samples },
      numerical: { I: numeric.I, IL: numeric.IL, IR: numeric.IR, Isink: numeric.Isink, n: numeric.n, spectralArea: numeric.spectralArea }, exact,
      currentErrorMicroamp: numeric.I - exact.I, occupationError: numeric.n - exact.n, conservationResidualMicroamp: numeric.IL + numeric.IR - numeric.Isink };
  });

  const etaConvergence = [0.1, 0.03, 0.01, 0.003, 0.001, 0].map(eta => {
    const parameters = { epsilon: 0.13, gammaL: 0.2, gammaR: 0.2, eta };
    const q = P.resonantLevel(parameters.epsilon, parameters);
    return { eta, T: q.T, R: q.R, loss: q.absorption, totalSpectralFunction: q.A,
      spectralFunctionFromPhysicalLeads: q.AL + q.AR, missingSpectralFunctionWithoutEtaTerm: q.A - q.AL - q.AR,
      probabilityLifetimeFemtoseconds: 0.6582119569 / (0.4 + 2 * eta) };
  });

  // An isolated level has a delta-function spectrum at \eta\to0. Its sampled peak
  // can land on a grid point or halfway between points: a smaller eta alone does
  // not produce a converged quadrature.
  const spectralGrid = [];
  for (const intervals of [100, 400, 1600, 6400]) {
    for (const fractionOfStep of [0, 0.5]) {
      const limit = 2;
      const eta = 0.003;
      const gridStep = 2 * limit / intervals;
      const offset = fractionOfStep * gridStep;
      const q = P.etaDiagnostics({ epsilon: 0, gammaL: 0, gammaR: 0, eta, limit, intervals, offset });
      const exact = lorentzianArea(0, 2 * eta, -limit + offset, limit + offset);
      spectralGrid.push({ eta, intervals, gridStep, fractionOfStep, lower: -limit + offset, upper: limit + offset,
        sampledArea: q.area, exactFiniteWindowArea: exact, quadratureError: q.area - exact, missingTail: 1 - exact });
    }
  }

  const spectralWindow = [0.5, 2, 10, 100].map(limit => {
    const parameters = { epsilon: 0.07, gammaL: 0.1, gammaR: 0.1, eta: 0.03, limit, intervals: Math.round(2 * limit / 0.002), offset: 0 };
    const q = P.etaDiagnostics(parameters);
    const exact = lorentzianArea(parameters.epsilon, 0.26, -limit, limit);
    return { limit, intervals: parameters.intervals, sampledArea: q.area, exactFiniteWindowArea: exact, quadratureError: q.area - exact, missingTail: 1 - exact };
  });

  const chainPlotSampling = [
    { sites: 9, hopping: 1, contact: 0.15, onsite: 0 },
    { sites: 9, hopping: 0.4, contact: 0.15, onsite: 1.5 },
    { sites: 9, hopping: 1.4, contact: 1.4, onsite: -1.5 }
  ].map(parameters => checkChainCurve(S.chainCurve(parameters), parameters));

  const fanoPlotSampling = [
    { sideEnergy: 1.25, coupling: 0.025 },
    { sideEnergy: -1.5, coupling: 0.025 },
    { sideEnergy: 0, coupling: 0.025 },
    { sideEnergy: -0.1, coupling: 0.05 },
    { sideEnergy: 0.7, coupling: 0.35 },
    { sideEnergy: 0, coupling: 1 },
    { sideEnergy: 1.5, coupling: 0 }
  ].map(parameters => checkFanoCurve(S.fanoCurve(parameters), parameters));

  const finiteTemperatureBias = [
    { epsilon: 0.05, bias: 0.1, kT: 0.005, fraction: 0.05 },
    { epsilon: -2, bias: -4, kT: 0.005, fraction: 0.95 },
    { epsilon: 0.35, bias: 0.7, kT: 0.2, fraction: 0.3 },
    { epsilon: 0, bias: 0, kT: 0.005, fraction: 0.5 }
  ].map(settings => {
    const { epsilon, bias, kT, fraction } = settings;
    const parameters = { epsilon, bias, kT, gammaL: 0.4 * fraction, gammaR: 0.4 * (1 - fraction), eta: 0, limit: 40, intervals: 32000 };
    const numerical = P.bias(parameters);
    const reference = finiteTemperatureReference(parameters);
    return { parameters, numerical: { I: numerical.I, n: numerical.n, spectralArea: numerical.spectralArea }, reference,
      currentErrorMicroamp: numerical.I - reference.I, occupationError: numerical.n - reference.n };
  });

  return roundTree({
    formatVersion: 1,
    runtime: { node: process.version },
    conventions: { energies: 'eV', currents: 'microampere per spin', lifetime: 'femtosecond; single wide-band pole, one-electron probability survival time \\(\\hbar/(\\Gamma+2\\eta)\\); maintained occupation determined separately', hopping: 'H[n,n+1] = -t', bias: 'muL = bias/2; muR = -bias/2', current: 'IL and IR enter device; Isink leaves device', eta: 'Finite eta is implemented as an empty absorbing channel of width 2 eta for current balance.', integration: 'Composite trapezoid. Zero-temperature bias integration inserts exact chemical-potential breakpoints and uses one-sided occupations on each segment. Spectral-grid diagnostics retain the specified uniform grid.', randomness: 'none' },
    independentReferences: { scattering: 'Outgoing-wave recurrence of the discrete Schrodinger equation, resolved into incident and reflected lead waves.', integration: 'Closed-form arctangent integrals of a Lorentzian, with the same finite energy window as the numerical integral.', finiteTemperature: 'Adaptive Simpson integration of analytic Lorentzian and tanh Fermi functions; integral tolerance 2e-12 per scalar integrand.', fano: 'T(E) = 0.25 (E-epsilon_s)^2 / ([E(E-epsilon_s)-v^2]^2 + 0.25(E-epsilon_s)^2); v=0 uses the central Lorentzian directly.', currentMicroampPerEV },
    matchedChain: { energies: openEnergies, hopping: 1, contact: 1, onsite: 0, eta: 0, rows: matched },
    mismatchedChain: { energies: openEnergies, eta: 0, rows: mismatch },
    zeroTemperatureBias: bias,
    finiteTemperatureBias,
    etaConvergence: { epsilon: 0.13, gammaL: 0.2, gammaR: 0.2, evaluationEnergy: 0.13, rows: etaConvergence },
    spectralGrid: { epsilon: 0, gammaL: 0, gammaR: 0, interpretation: 'Compare sampled area with exact finite-window area, not with one alone. Grid error and omitted tails are distinct.', rows: spectralGrid },
    spectralWindow: { epsilon: 0.07, gammaL: 0.1, gammaR: 0.1, eta: 0.03, gridStep: 0.002, rows: spectralWindow },
    chainPlotSampling: { interpolation: 'Piecewise linear between samples; do not cross band-edge path breaks.', reference: 'Independent outgoing-wave recurrence on a dense uniform mesh, with separately refined local maxima and half-height flanks.', rows: chainPlotSampling },
    fanoPlotSampling: { interpolation: 'Piecewise linear between samples, with seeds at the exact zero and both unit-transmission energies when v>0.', reference: 'Independent rational transmission formula on a dense uniform mesh and six interior points per rendered segment. Reported maxima are sampled errors, not uniform mathematical bounds.', rows: fanoPlotSampling }
  });
}

if (require.main === module) {
  const report = buildReport();
  console.log('Quantum transport numerical checks (eV, microamperes per spin; no randomness)');
  console.log(`Node ${process.version}; composite trapezoid integration`);
  console.log('\nMatched chain');
  console.table(report.matchedChain.rows);
  console.log('Independent wave matching versus Green functions');
  console.table(report.mismatchedChain.rows.map(row => ({ sites: row.parameters.sites, maximumTransmissionError: row.maximumTransmissionError, maximumReflectionError: row.maximumReflectionError })));
  console.log('Zero-temperature finite-window bias integrals');
  console.table(report.zeroTemperatureBias.map(row => ({ eta: row.parameters.eta, currentMicroamp: row.numerical.I, currentErrorMicroamp: row.currentErrorMicroamp, occupation: row.numerical.n, occupationError: row.occupationError, conservationResidualMicroamp: row.conservationResidualMicroamp })));
  console.log('Finite-temperature bias versus independent adaptive Simpson integration');
  console.table(report.finiteTemperatureBias.map(row => ({ epsilon: row.parameters.epsilon, bias: row.parameters.bias, kT: row.parameters.kT, currentErrorMicroamp: row.currentErrorMicroamp, occupationError: row.occupationError })));
  console.log('Finite-eta convergence on an open resonance, Gamma_L = Gamma_R = 0.2 eV');
  console.table(report.etaConvergence.rows.map(row => ({ eta: row.eta, T: row.T, R: row.R, loss: row.loss, spectralWeightMissingWithoutEta: row.missingSpectralFunctionWithoutEtaTerm })));
  console.log('Spectral grid convergence for an isolated level, eta = 0.003 eV');
  console.table(report.spectralGrid.rows);
  console.log('Finite-window spectral weight, delta E = 0.002 eV');
  console.table(report.spectralWindow.rows);
  console.log('Chain plot versus independent dense wave-matching scans');
  console.table(report.chainPlotSampling.rows.map(row => ({ hopping: row.parameters.hopping, contact: row.parameters.contact, onsite: row.parameters.onsite, plotSamples: row.plotSamples, resolvedPeaks: row.peaks.length, sampledMaximumError: row.sampledMaximumInterpolationError, libraryReferenceDifference: row.maximumLibraryReferenceDifference })));
  console.log('Fano plot versus independent rational transmission formula');
  console.table(report.fanoPlotSampling.rows.map(row => ({ sideEnergy: row.parameters.sideEnergy, coupling: row.parameters.coupling, plotSamples: row.plotSamples, sampledMaximumError: row.sampledMaximumInterpolationError })));
  if (process.argv.includes('--write')) {
    const destination = path.join(__dirname, '..', 'docs', 'numerical-checks.json');
    fs.mkdirSync(path.dirname(destination), { recursive: true });
    fs.writeFileSync(destination, JSON.stringify(report, null, 2) + '\n');
    console.log(`Saved ${path.relative(process.cwd(), destination)}`);
  }
}

module.exports = { scatteringReference, lorentzianArea, zeroTemperatureReference, finiteTemperatureReference, currentMicroampPerEV, interpolateCurve, checkChainCurve, checkFanoCurve, buildReport };
