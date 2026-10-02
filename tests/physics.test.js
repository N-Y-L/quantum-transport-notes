'use strict';
// Copyright (c) 2026 Neil Yuanting Li. MIT license; see LICENSE.

const test = require('node:test');
const assert = require('node:assert/strict');
const P = require('../physics.js');
const { scatteringReference, lorentzianArea, zeroTemperatureReference, finiteTemperatureReference, currentMicroampPerEV, checkChainCurve, checkFanoCurve, interpolateCurve } = require('../scripts/check-numerics.js');

function close(actual, expected, tolerance = 1e-11, message = '') {
  assert.ok(Number.isFinite(actual), `${message}: nonfinite result ${actual}`);
  assert.ok(Math.abs(actual - expected) <= tolerance, `${message}: got ${actual}, expected ${expected}, tolerance ${tolerance}`);
}
const cMultiply = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const cConjugate = a => [a[0], -a[1]];

// Compare every matrix entry of i(G-G†) with G W G†. Here W is a diagonal
// broadening supplied from the model specification, including 2η on every site.
function assertSpectralIdentity(G, widths, tolerance = 2e-10) {
  for (let i = 0; i < G.length; i++) {
    for (let j = 0; j < G.length; j++) {
      const spectral = [-(G[i][j][1] + G[j][i][1]), G[i][j][0] - G[j][i][0]];
      const product = [0, 0];
      for (let k = 0; k < G.length; k++) {
        const term = cMultiply(G[i][k], cConjugate(G[j][k]));
        product[0] += widths[k] * term[0];
        product[1] += widths[k] * term[1];
      }
      close(spectral[0], product[0], tolerance, `spectral real (${i},${j})`);
      close(spectral[1], product[1], tolerance, `spectral imaginary (${i},${j})`);
    }
  }
}

test('resonant level agrees with the independently evaluated Lorentzian', () => {
  for (const epsilon of [-0.37, 0, 0.28]) {
    for (const [gammaL, gammaR] of [[0.5, 0.5], [0.12, 0.71], [0, 0.4]]) {
      for (const eta of [0, 0.017, 0.2]) {
        for (const E of [-1.3, -0.17, 0, 0.43, 1.1]) {
          const q = P.resonantLevel(E, { epsilon, gammaL, gammaR, eta });
          const a = (gammaL + gammaR) / 2 + eta;
          const denominator = (E - epsilon) ** 2 + a * a;
          close(q.T, gammaL * gammaR / denominator);
          close(q.A, 2 * a / denominator);
          close(q.G[0], (E - epsilon) / denominator);
          close(q.G[1], -a / denominator);
          close(q.absorption, 2 * eta * gammaL / denominator);
        }
      }
    }
  }
});

test('resonance height, FWHM and lifetime factors use the same Gamma convention', () => {
  const p = { epsilon: 0.17, gammaL: 0.2, gammaR: 0.6, eta: 0.03 };
  const width = p.gammaL + p.gammaR + 2 * p.eta;
  const peak = P.resonantLevel(p.epsilon, p);
  close(peak.T, 4 * p.gammaL * p.gammaR / width ** 2);
  close(P.resonantLevel(p.epsilon + width / 2, p).T, peak.T / 2);
  close(P.resonantLevel(p.epsilon - width / 2, p).T, peak.T / 2);
  const hbar = 0.6582119569;
  const probabilityLifetime = hbar / width;
  const amplitudeLifetime = 2 * hbar / width;
  close(Math.exp(-width * probabilityLifetime / hbar), 1 / Math.E);
  close(Math.exp(-width * amplitudeLifetime / (2 * hbar)), 1 / Math.E);
  close(amplitudeLifetime, 2 * probabilityLifetime);
});

test('one-channel resonant scattering is unitary at eta=0 and subunitary for an empty sink', () => {
  for (const eta of [0, 0.001, 0.07, 0.4]) {
    for (const E of [-2, -0.31, 0.12, 0.9]) {
      const q = P.resonantLevel(E, { epsilon: 0.12, gammaL: 0.17, gammaR: 0.63, eta });
      assert.ok(q.T >= 0 && q.T <= 1 + 1e-12);
      assert.ok(q.R >= 0 && q.R <= 1 + 1e-12);
      close(q.T + q.R + q.absorption, 1);
      close(q.T + q.Rright + q.absorptionRight, 1);
      close(q.A, q.AL + q.AR + q.Aeta);
    }
  }
});

test('omitting the finite-eta spectral term demonstrably breaks the identity', () => {
  const q = P.resonantLevel(0, { epsilon: 0, gammaL: 0.2, gammaR: 0.2, eta: 0.1 });
  close(q.A - q.AL - q.AR, q.Aeta);
  assert.ok(q.Aeta / q.A > 0.3, 'missing sink spectral weight must be visible');
  assert.ok(1 - q.T - q.R > 0.4, 'finite eta must not be presented as conserved two-terminal transport');
});

test('semi-infinite lead branch is causal, continuous at both edges and solves its quadratic', () => {
  for (const hopping of [0.6, 1, 1.7]) {
    const contact = 0.43;
    const onsite = 0.27;
    for (const scaledEnergy of [-6, -2.001, -2, -1.5, 0, 1.5, 2, 2.001, 6]) {
      const E = onsite + scaledEnergy * hopping;
      const q = P.surfaceSelfEnergy(E, { hopping, contact, onsite });
      const g = [q.re / contact ** 2, q.im / contact ** 2];
      const square = cMultiply(g, g);
      close(hopping ** 2 * square[0] - (E - onsite) * g[0] + 1, 0, 2e-12, 'surface quadratic real part');
      close(hopping ** 2 * square[1] - (E - onsite) * g[1], 0, 2e-12, 'surface quadratic imaginary part');
      assert.ok(q.im <= 1e-14, 'retarded self-energy must have nonpositive imaginary part');
      close(q.gamma, -2 * q.im);
      if (Math.abs(scaledEnergy) > 2) {
        close(q.im, 0);
        assert.ok(Math.abs(g[0]) < 1 / hopping, 'outside-band solution must select the decaying branch');
        assert.equal(Math.sign(g[0]), Math.sign(E - onsite));
      }
    }
    for (const edge of [-2, 2]) {
      const boundary = P.surfaceSelfEnergy(onsite + edge * hopping, { hopping, contact, onsite });
      const outside = P.surfaceSelfEnergy(onsite + edge * hopping * (1 + 1e-9), { hopping, contact, onsite });
      assert.ok(Math.abs(boundary.re - outside.re) < 2e-5);
    }
  }
});

test('surface self-energy has the correct far-band 1/E limit and contact-square scaling', () => {
  for (const sign of [-1, 1]) {
    const E = sign * 10000;
    const q = P.surfaceSelfEnergy(E, { hopping: 1.3, contact: 0.7, onsite: 0.2 });
    close(q.re * (E - 0.2) / 0.7 ** 2, 1, 5e-8);
  }
  const one = P.surfaceSelfEnergy(0.3, { hopping: 1.1, contact: 0.4, onsite: 0 });
  const two = P.surfaceSelfEnergy(0.3, { hopping: 1.1, contact: 0.8, onsite: 0 });
  close(two.re, 4 * one.re);
  close(two.im, 4 * one.im);
});

test('a matched finite chain transmits exactly one open channel independently of length', () => {
  for (const sites of [1, 2, 3, 8, 16]) {
    for (const E of [-1.99, -1.73, -0.8, 0, 0.31, 1.61, 1.99]) {
      const q = P.chain(E, { sites, hopping: 1, contact: 1, onsite: 0, eta: 0 });
      close(q.T, 1, 2e-10, `matched chain sites=${sites}, E=${E}`);
      close(q.R, 0, 2e-10);
      close(q.absorption, 0);
    }
  }
});

test('mismatched chain matches independent discrete wave scattering', () => {
  const parameterSets = [
    { sites: 1, hopping: 1, contact: 0.4, onsite: 0.2, leadHopping: 1, leadOnsite: 0 },
    { sites: 3, hopping: 0.83, contact: 0.57, onsite: 0.21, leadHopping: 1, leadOnsite: 0 },
    { sites: 7, hopping: 1.17, contact: 0.76, onsite: -0.14, leadHopping: 1, leadOnsite: 0 },
    { sites: 4, hopping: 0.7, contact: 0.6, onsite: -0.17, leadHopping: 0.9, leadOnsite: 0.13 }
  ];
  for (const p of parameterSets) {
    for (const E of [-1.4, -0.72, -0.17, 0, 0.46, 1.3]) {
      const wave = scatteringReference(E, p);
      const green = P.chain(E, { ...p, eta: 0 });
      close(wave.T + wave.R, 1, 2e-11, 'reference wave flux');
      close(green.T, wave.T, 2e-10, `T, N=${p.sites}, E=${E}`);
      close(green.R, wave.R, 2e-10, `R, N=${p.sites}, E=${E}`);
    }
  }
});

test('chain spectral identity holds as a complex matrix, including absorbing channels', () => {
  for (const sites of [1, 2, 5]) {
    for (const eta of [0, 0.023]) {
      const p = { sites, hopping: 0.84, contact: 0.69, onsite: 0.16, eta };
      for (const E of [-1.6, -0.42, 0.37, 1.5]) {
        const q = P.chain(E, p);
        const leadGamma = 2 * p.contact ** 2 * Math.sqrt(1 - (E / 2) ** 2);
        const widths = Array(sites).fill(2 * eta);
        widths[0] += leadGamma;
        widths[sites - 1] += leadGamma;
        assertSpectralIdentity(q.G, widths);
        close(q.A, q.AL + q.AR + q.Aeta, 2e-10);
        close(q.T + q.R + q.absorption, 1, 2e-10);
        assert.ok(q.T >= 0 && q.T <= 1 + 2e-10);
        assert.ok(q.absorption >= -2e-10);
      }
    }
  }
});

test('closed lead channels have zero transmission and no incident-wave reflection probability', () => {
  for (const E of [-3.2, -2, 2, 3.2]) {
    const q = P.chain(E, { sites: 4, hopping: 0.7, contact: 0.5, onsite: 0.1, eta: 0.01 });
    close(q.T, 0);
    close(q.absorption, 0);
    assert.equal(q.channels, 0);
    assert.equal(q.R, null);
  }
});

test('Fano transmission has the exact interference zero and agrees with its rational formula', () => {
  const p = { epsilon: -0.11, sideEnergy: 0.23, coupling: 0.37, gammaL: 0.21, gammaR: 0.49, eta: 0 };
  for (const E of [-1.2, -0.31, 0, p.sideEnergy, 0.49, 1.1]) {
    const delta = E - p.sideEnergy;
    const numerator = p.gammaL * p.gammaR * delta ** 2;
    const denominator = ((E - p.epsilon) * delta - p.coupling ** 2) ** 2 + ((p.gammaL + p.gammaR) * delta / 2) ** 2;
    const q = P.fano(E, p);
    close(q.T, numerator / denominator, 2e-11);
    close(q.T + q.R, 1, 2e-11);
  }
  const zero = P.fano(p.sideEnergy, p);
  close(zero.T, 0, 1e-13);
  close(zero.R, 1, 1e-13);
  assert.ok(P.fano(p.sideEnergy, { ...p, eta: 0.03 }).T > 0, 'an absorber on the side orbital fills the exact zero');
});

test('decoupling the side orbital recovers central resonant transmission, even at the side energy', () => {
  const p = { epsilon: -0.1, sideEnergy: 0.23, coupling: 0, gammaL: 0.2, gammaR: 0.5, eta: 0 };
  for (const E of [-1, 0, p.sideEnergy, 0.7]) {
    close(P.fano(E, p).T, P.resonantLevel(E, p).T, 2e-11);
  }
  const pole = P.fano(p.sideEnergy, p);
  assert.equal(pole.A, null, 'an isolated delta function must not be reported as finite sampled DOS');
  assert.deepEqual(pole.boundStates, [{ energy: p.sideEnergy, weight: 1, site: 1 }]);
  assert.equal(pole.spectralIncludesBoundStates, false);
});

test('Fano spectral and flux identities include eta on both orbitals', () => {
  const p = { epsilon: 0.1, sideEnergy: -0.23, coupling: 0.4, gammaL: 0.2, gammaR: 0.5, eta: 0.02 };
  for (const E of [-0.7, -0.23, 0, 0.5]) {
    const q = P.fano(E, p);
    assertSpectralIdentity(q.G, [p.gammaL + p.gammaR + 2 * p.eta, 2 * p.eta]);
    close(q.A, q.AL + q.AR + q.Aeta, 2e-10);
    close(q.T + q.R + q.absorption, 1, 2e-10);
    assert.ok(q.T >= 0 && q.T <= 1 + 1e-12);
  }
});

test('fully isolated Fano orbitals report both delta-function states and distributional poles', () => {
  // H=[[0.3,-0.4],[-0.4,-0.3]] has eigenvalues ±0.5 and normalized
  // eigenvectors (1,2)/sqrt(5), (2,-1)/sqrt(5), respectively.
  const cases = [
    { epsilon: 0.3, sideEnergy: -0.3, coupling: 0.4, energies: [-0.5, 0.5], weights: [[0.2, 0.8], [0.8, 0.2]] },
    { epsilon: 0, sideEnergy: 0, coupling: 0.5, energies: [-0.5, 0.5], weights: [[0.5, 0.5], [0.5, 0.5]] },
    { epsilon: 0.3, sideEnergy: -0.3, coupling: 0, energies: [-0.3, 0.3], weights: [[0, 1], [1, 0]] },
    { epsilon: 0.3, sideEnergy: 0.3, coupling: 0, energies: [0.3, 0.3], weights: [[1, 0], [0, 1]] }
  ];
  for (const example of cases) {
    const p = { epsilon: example.epsilon, sideEnergy: example.sideEnergy, coupling: example.coupling, gammaL: 0, gammaR: 0, eta: 0 };
    const q = P.fano(0.1, p);
    close(q.A, 0, 0, 'regular spectrum away from the isolated poles');
    assert.equal(q.spectralIncludesBoundStates, false);
    assert.equal(q.boundStates.length, 2);
    for (let i = 0; i < 2; i++) {
      const state = q.boundStates[i];
      close(state.energy, example.energies[i]);
      close(state.weight, 1);
      for (let site = 0; site < 2; site++) close(state.siteWeights[site], example.weights[i][site]);
      const pole = P.fano(state.energy, p);
      assert.equal(pole.singular, true);
      assert.equal(pole.G, null);
      assert.equal(pole.A, null, 'a delta-function pole has no ordinary finite spectral value');
      assert.equal(pole.spectralIncludesBoundStates, false);
      close(pole.T, 0);
      close(pole.R, 1);
      close(pole.absorption, 0);
    }
    for (let site = 0; site < 2; site++) {
      close(q.boundStates.reduce((total, state) => total + state.siteWeights[site], 0), 1, 1e-12, 'one integrated state per orbital');
    }
  }
});

test('broadening an isolated Fano device recovers its eigenstate spectral projectors', () => {
  const p = { epsilon: 0.3, sideEnergy: -0.3, coupling: 0.4, gammaL: 0, gammaR: 0, eta: 0.02 };
  for (const E of [-1, -0.5, 0, 0.5, 1]) {
    const lower = 2 * p.eta / ((E + 0.5) ** 2 + p.eta ** 2);
    const upper = 2 * p.eta / ((E - 0.5) ** 2 + p.eta ** 2);
    const q = P.fano(E, p);
    close(q.spectral[0], 0.2 * lower + 0.8 * upper, 2e-10);
    close(q.spectral[1], 0.8 * lower + 0.2 * upper, 2e-10);
    close(q.A, lower + upper, 2e-10);
    assert.equal(q.spectralIncludesBoundStates, true);
    assert.deepEqual(q.boundStates, []);
  }
});

test('coupled two-orbital spectral trace has two units of total weight, with exact finite-window tails', () => {
  // For equal onsite energies and v>Gamma/4, the two effective-Hamiltonian
  // poles are ±sqrt(v²−Gamma²/16)−i(Gamma/4+eta). Tr G has residue one
  // at each pole, even though the individual site spectra are not this sum.
  const p = { epsilon: 0, sideEnergy: 0, coupling: 0.5, gammaL: 0.3, gammaR: 0.5, eta: 0.02 };
  const gamma = p.gammaL + p.gammaR;
  const pole = Math.sqrt(p.coupling ** 2 - gamma ** 2 / 16);
  const halfWidth = gamma / 4 + p.eta;
  for (const E of [-1, -pole, 0, pole, 1]) {
    const exactA = 2 * halfWidth / ((E - pole) ** 2 + halfWidth ** 2)
      + 2 * halfWidth / ((E + pole) ** 2 + halfWidth ** 2);
    close(P.fano(E, p).A, exactA, 2e-11);
  }
  const lower = -4, upper = 4;
  const area = P.integrate(E => P.fano(E, p).A / (2 * Math.PI), lower, upper, 16000);
  const exact = lorentzianArea(pole, 2 * halfWidth, lower, upper) + lorentzianArea(-pole, 2 * halfWidth, lower, upper);
  close(area, exact, 3e-10);
  assert.ok(exact > 1.92 && exact < 2, 'both broadened orbitals contribute, with finite-window tails excluded');
});

test('Fermi function has correct limits, symmetry and half-filled zero-temperature step', () => {
  close(P.fermi(0.3, 0.3, 0), 0.5);
  close(P.fermi(0.2, 0.3, 0), 1);
  close(P.fermi(0.4, 0.3, 0), 0);
  for (const x of [0, 0.02, 0.2, 20]) close(P.fermi(0.3 + x, 0.3, 0.025) + P.fermi(0.3 - x, 0.3, 0.025), 1);
  close(P.fermi(-100, 0, 0.001), 1);
  close(P.fermi(100, 0, 0.001), 0);
});

test('zero-temperature current and occupations match exact finite-window arctangent integrals', () => {
  for (const bias of [-0.5, 0, 0.5]) {
    for (const eta of [0, 0.03]) {
      const p = { epsilon: 0.125, gammaL: 0.18, gammaR: 0.32, eta, bias, kT: 0, limit: 8, intervals: 65536 };
      const exact = zeroTemperatureReference(p);
      const q = P.bias(p);
      for (const key of ['I', 'IL', 'IR', 'Isink']) close(q[key], exact[key], 3e-6, key);
      close(q.n, exact.n, 5e-8, 'occupation');
      close(q.spectralArea, exact.spectralArea, 2e-10, 'spectral area');
    }
  }
});

test('finite-temperature currents conserve flux only after the modeled sink is included', () => {
  for (const eta of [0, 0.04]) {
    const p = { epsilon: 0.19, gammaL: 0.13, gammaR: 0.47, bias: 0.38, kT: 0.027, eta, limit: 10, intervals: 12000 };
    const q = P.bias(p);
    close(q.IL + q.IR, q.Isink, 2e-10);
    assert.ok(q.n >= 0 && q.n <= 1, 'spinless occupation must lie in [0,1]');
    if (eta === 0) {
      close(q.IL, q.I);
      close(q.IR, -q.I);
      close(q.Isink, 0);
    } else {
      assert.ok(q.Isink > 1, 'physical-lead current imbalance from the empty sink must be resolved');
      assert.ok(Math.abs(q.IL + q.IR) > 1);
    }
  }
});

test('zero-temperature step integration handles shifted windows, endpoint steps, and steps outside the window', () => {
  for (const bias of [-2.2, -1.9, -0.5, 0, 0.5, 1.9, 2.2]) {
    const p = { epsilon: 0.137, gammaL: 0.18, gammaR: 0.32, eta: 0.03, bias, kT: 0, limit: 0.6, offset: 0.35, intervals: 10000 };
    const q = P.bias(p), exact = zeroTemperatureReference(p);
    close(q.I, exact.I, 1e-6, `current, bias=${bias}`);
    close(q.n, exact.n, 2e-8, `occupation, bias=${bias}`);
    close(q.spectralArea, exact.spectralArea, 2e-8);
    assert.ok(q.actualIntervals >= p.intervals && q.actualIntervals <= p.intervals + 2);
  }
});

test('finite-temperature interface grids agree with independent adaptive Simpson reference', () => {
  for (const [epsilon, bias, kT, fraction] of [[0.05, 0.1, 0.005, 0.05], [-2, -4, 0.005, 0.95], [0.35, 0.7, 0.2, 0.3]]) {
    const p = { epsilon, bias, kT, gammaL: 0.4 * fraction, gammaR: 0.4 * (1 - fraction), eta: 0, limit: 40, intervals: 32000 };
    const q = P.bias(p), reference = finiteTemperatureReference(p);
    close(q.I, reference.I, 1e-7, 'finite-temperature current');
    close(q.n, reference.n, 1e-8, 'finite-temperature occupation');
  }
  // At particle-hole symmetry and equilibrium, f(E)+f(-E)=1 fixes n to
  // half the finite-window spectral weight at every positive temperature.
  const p = { epsilon: 0, gammaL: 0.1, gammaR: 0.3, eta: 0, bias: 0, kT: 0.005, limit: 40, intervals: 32000 };
  const reference = finiteTemperatureReference(p);
  close(reference.I, 0);
  close(reference.n, reference.spectralArea / 2, 2e-12, 'independent quadrature symmetry');
  close(P.bias(p).n, reference.n, 1e-9);
});

test('small-bias conductance approaches (e²/h)T(EF) for a spinless channel', () => {
  for (const bias of [0.0002, 0.0002037, -0.0002037]) {
    const p = { epsilon: 0.17, gammaL: 0.21, gammaR: 0.39, eta: 0, bias, kT: 0, limit: 1, intervals: 200000 };
    const q = P.bias(p);
    const conductance = q.I / p.bias;
    const expected = currentMicroampPerEV * P.resonantLevel(0, p).T;
    close(conductance, expected, 2e-5, `bias=${bias}`);
  }
});

test('trapezoid quadrature is exact for linear data and converges for a smooth polynomial', () => {
  close(P.integrate(x => 2 * x + 3, -2, 4, 17), 30, 1e-12);
  close(P.trapz([0, 0.2, 0.9, 1], [1, 1.4, 2.8, 3]), 2);
  const coarseError = Math.abs(P.integrate(x => x * x, 0, 1, 20) - 1 / 3);
  const fineError = Math.abs(P.integrate(x => x * x, 0, 1, 40) - 1 / 3);
  close(coarseError / fineError, 4, 2e-10);
});

test('spectral finite-window normalization agrees with the exact arctangent, not automatically one', () => {
  const p = { epsilon: 0.07, gammaL: 0.1, gammaR: 0.1, eta: 0.03, limit: 2, intervals: 16000, offset: 0 };
  const q = P.etaDiagnostics(p);
  const exact = lorentzianArea(p.epsilon, p.gammaL + p.gammaR + 2 * p.eta, -p.limit, p.limit);
  close(q.analyticArea, exact);
  close(q.area, exact, 2e-10);
  assert.ok(1 - exact > 0.04, 'finite-window tails must be non-negligible in this check');
  const largerWindow = P.etaDiagnostics({ ...p, limit: 20, intervals: 160000 });
  assert.ok(largerWindow.area > q.area);
  assert.ok(largerWindow.area < 1);
});

test('a narrow isolated pole exposes grid-offset error; resolving its width restores the sum rule', () => {
  const base = { epsilon: 0, gammaL: 0, gammaR: 0, eta: 0.003, limit: 2 };
  const coarse = P.etaDiagnostics({ ...base, intervals: 100, offset: 0 });
  const coarseShifted = P.etaDiagnostics({ ...base, intervals: 100, offset: 0.02 });
  assert.ok(coarse.area > 4, 'grid point through a narrow pole must substantially overestimate its integral');
  assert.ok(coarseShifted.area < 0.3, 'half-step shift must miss much of the same narrow pole');
  for (const offset of [0, 0.0003125]) {
    const fine = P.etaDiagnostics({ ...base, intervals: 6400, offset });
    close(fine.area, lorentzianArea(0, 0.006, -2 + offset, 2 + offset), 5e-11);
  }
});

test('eta→0 converges for an already open resonance but fails as an unresolved delta-peak quadrature', () => {
  const p = { epsilon: 0, gammaL: 0.2, gammaR: 0.2 };
  let previousError = Infinity;
  for (const eta of [0.1, 0.03, 0.01, 0.003, 0.001, 0]) {
    const error = 1 - P.resonantLevel(0, { ...p, eta }).T;
    assert.ok(error >= -1e-13 && error < previousError);
    previousError = error;
  }
  const narrow = P.etaDiagnostics({ epsilon: 0, gammaL: 0, gammaR: 0, eta: 0.001, limit: 2, intervals: 400, offset: 0 });
  const broad = P.etaDiagnostics({ epsilon: 0, gammaL: 0, gammaR: 0, eta: 0.03, limit: 2, intervals: 400, offset: 0 });
  assert.ok(Math.abs(narrow.area - narrow.analyticArea) > 1);
  assert.ok(Math.abs(broad.area - broad.analyticArea) < 1e-6);
});

test('invalid broadening, nonfinite energies and malformed discretizations are rejected', () => {
  assert.throws(() => P.resonantLevel(0, { gammaL: -0.1 }), RangeError);
  assert.throws(() => P.resonantLevel(0, { eta: -0.1 }), RangeError);
  assert.throws(() => P.resonantLevel(NaN, {}));
  assert.throws(() => P.chain(0, { sites: 2.5 }), RangeError);
  assert.throws(() => P.chain(0, { hopping: 0 }), RangeError);
  assert.throws(() => P.surfaceSelfEnergy(0, { contact: -1 }), RangeError);
  assert.throws(() => P.etaDiagnostics({ limit: 0 }), RangeError);
  assert.throws(() => P.etaDiagnostics({ intervals: 0 }), RangeError);
  assert.throws(() => P.resonantLevel(0, { epsilon: 0, gammaL: 0, gammaR: 0, eta: 0 }), RangeError);
});

for (const parameters of [
  { sites: 9, hopping: 1, contact: 0.15, onsite: 0 },
  { sites: 9, hopping: 0.4, contact: 0.15, onsite: 1.5 },
  { sites: 9, hopping: 1.4, contact: 1.4, onsite: -1.5 }
]) {
  test(`chain plot resolves peaks and flanks against independent wave matching: t=${parameters.hopping}, contact=${parameters.contact}, onsite=${parameters.onsite}`, t => {
    const S = require('../sampling.js');
    const curve = S.chainCurve(parameters);
    for (let i = 0; i < curve.length; i++) {
      assert.ok(Number.isFinite(curve[i][0]), 'energy samples must be finite');
      if (i > 0) assert.ok(curve[i][0] >= curve[i - 1][0], 'energy samples must be ordered');
      assert.ok(curve[i][1] === null || (Number.isFinite(curve[i][1]) && curve[i][1] >= -1e-12 && curve[i][1] <= 1 + 1e-10));
    }
    const check = checkChainCurve(curve, parameters);
    assert.ok(check.sampledMaximumInterpolationError < 0.002, `sampled interpolation error ${check.sampledMaximumInterpolationError} at E=${check.worstEnergy}`);
    assert.ok(check.maximumLibraryReferenceDifference < 2e-9, 'wave reference must agree independently with the matrix transport library');
    assert.ok(check.peaks.length >= 4, 'challenging spectra must contain their resolved resonances');
    for (const peak of check.peaks) {
      assert.ok(peak.maximumLocalInterpolationError < 0.002, `peak/flank error near E=${peak.energy}`);
      assert.ok(peak.plottedPeak > 0.998, `unit-transmission resonance at E=${peak.energy} must be visible`);
    }
    if (parameters.hopping === 1 && parameters.contact === 0.15 && parameters.onsite === 0) {
      assert.equal(check.peaks.length, 9, 'all nine weak-contact resonances must be resolved');
      assert.ok(Math.min(...check.peaks.map(peak => peak.halfMaximumWidth)) < 0.0006, 'reference must exercise the narrow outer resonances');
    }
    t.diagnostic(`${curve.length} plotted samples; ${check.peaks.length} peaks; sampled maximum |T_plot−T_wave|=${check.sampledMaximumInterpolationError.toExponential(3)}`);
  });
}

test('matched-chain plots preserve unit transmission inside the lead band and break at thresholds', () => {
  const S = require('../sampling.js');
  const curve = S.chainCurve({ sites: 9, hopping: 1, contact: 1, onsite: 0 });
  close(curve[0][0], -2.4);
  close(curve.at(-1)[0], 2.4);
  for (const [E, T] of curve) {
    if (T === null) continue;
    close(T, Math.abs(E) < 2 ? 1 : 0, 2e-9);
  }
  assert.ok(curve.some(([E, T]) => E === -2 && T === null), 'left channel threshold requires a path break');
  assert.ok(curve.some(([E, T]) => E === 2 && T === null), 'right channel threshold requires a path break');
  for (const E of [-1.99, -1, 0, 0.7, 1.99]) close(interpolateCurve(curve, E), 1, 2e-9);
});

test('Fano curves resolve displaced narrow resonances, their tails, and the exact zero', t => {
  const S = require('../sampling.js');
  for (const parameters of [
    { sideEnergy: 1.25, coupling: 0.025 },
    { sideEnergy: -1.5, coupling: 0.025 },
    { sideEnergy: 0, coupling: 0.025 },
    { sideEnergy: -0.1, coupling: 0.05 },
    { sideEnergy: 0.7, coupling: 0.35 },
    { sideEnergy: 0, coupling: 1 }
  ]) {
    const curve = S.fanoCurve(parameters);
    const check = checkFanoCurve(curve, parameters);
    close(curve[0][0], -3);
    close(curve.at(-1)[0], 3);
    for (let i = 1; i < curve.length; i++) assert.ok(curve[i][0] > curve[i - 1][0]);
    close(check.atSideEnergy, 0, 1e-13);
    for (const peak of check.peaks) close(peak.plottedTransmission, 1, 2e-10);
    assert.ok(check.sampledMaximumInterpolationError < 0.0007,
      `Fano interpolation error ${check.sampledMaximumInterpolationError} at E=${check.worstEnergy}`);
    t.diagnostic(`sideEnergy=${parameters.sideEnergy}, coupling=${parameters.coupling}: ${curve.length} samples, max sampled error=${check.sampledMaximumInterpolationError.toExponential(3)}`);
  }
});

test('disconnected Fano curves retain the central Lorentzian without a false transmission notch', () => {
  const S = require('../sampling.js');
  for (const sideEnergy of [-1.5, 0, 1.25]) {
    const p = { sideEnergy, coupling: 0 };
    const curve = S.fanoCurve(p), check = checkFanoCurve(curve, p);
    close(check.atSideEnergy, 0.25 / (sideEnergy ** 2 + 0.25), 0.0007);
    close(interpolateCurve(curve, 0), 1);
    assert.ok(check.sampledMaximumInterpolationError < 0.0007);
  }
});
