/* Copyright (c) 2026 Neil Yuanting Li. MIT license; see LICENSE.
 * One spin channel; energies in eV, spectral functions in 1/eV.
 * Complex numbers are [real, imaginary]. No numerical eta is put in lead Sigma.
 */
(function (root, factory) {
  "use strict";
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.TransportPhysics = factory();
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";

  const constants = Object.freeze({
    hbarEVfs: 0.6582119569,
    currentMicroampPerEV: 38.74045865,
    kBoltzmannEV: 8.617333262145e-5
  });

  function finite(value, name) {
    if (typeof value !== "number" || !Number.isFinite(value)) {
      throw new TypeError(name + " must be a finite number");
    }
    return value;
  }
  function nonnegative(value, name) {
    finite(value, name);
    if (value < 0) throw new RangeError(name + " must be nonnegative");
    return value;
  }
  function positive(value, name) {
    finite(value, name);
    if (value <= 0) throw new RangeError(name + " must be positive");
    return value;
  }
  function integer(value, name, min, max) {
    finite(value, name);
    if (!Number.isInteger(value) || value < min || value > max) {
      throw new RangeError(name + " must be an integer from " + min + " to " + max);
    }
    return value;
  }
  function abs2(z) { return z[0] * z[0] + z[1] * z[1]; }
  function sum(values) { return values.reduce(function (a, b) { return a + b; }, 0); }

  // Stable complex division avoids squaring a very large denominator.
  function divide(a, b) {
    if (Math.abs(b[0]) >= Math.abs(b[1])) {
      const ratio = b[1] / b[0];
      const denominator = b[0] + b[1] * ratio;
      return [(a[0] + a[1] * ratio) / denominator,
        (a[1] - a[0] * ratio) / denominator];
    }
    const ratio = b[0] / b[1];
    const denominator = b[1] + b[0] * ratio;
    return [(a[0] * ratio + a[1]) / denominator,
      (a[1] * ratio - a[0]) / denominator];
  }

  function levelParameters(p) {
    p = p || {};
    const result = {
      epsilon: p.epsilon === undefined ? 0 : finite(p.epsilon, "epsilon"),
      gammaL: p.gammaL === undefined ? 0.5 : nonnegative(p.gammaL, "gammaL"),
      gammaR: p.gammaR === undefined ? 0.5 : nonnegative(p.gammaR, "gammaR"),
      eta: p.eta === undefined ? 0 : nonnegative(p.eta, "eta")
    };
    result.gamma = result.gammaL + result.gammaR;
    result.gammaEta = 2 * result.eta;
    result.width = result.gamma + result.gammaEta;
    return result;
  }

  // The eta term has exactly the same retarded part as an empty absorbing
  // reservoir with Gamma_eta=2 eta. This is a diagnostic interpretation, not
  // an assertion that a numerical regulator models dephasing or a lifetime.
  function evaluateLevel(E, p) {
    const detuning = E - p.epsilon;
    const halfWidth = p.width / 2;
    if (detuning === 0 && halfWidth === 0) {
      throw new RangeError("An isolated level at its pole has delta-function spectral weight; choose eta > 0 to sample it");
    }
    const G = divide([1, 0], [detuning, halfWidth]);
    const probability = abs2(G);
    const AL = p.gammaL * probability;
    const AR = p.gammaR * probability;
    const Aeta = p.gammaEta * probability;
    return {
      G: G,
      T: p.gammaL * AR,
      R: Math.pow(1 + p.gammaL * G[1], 2) + Math.pow(p.gammaL * G[0], 2),
      Rright: Math.pow(1 + p.gammaR * G[1], 2) + Math.pow(p.gammaR * G[0], 2),
      absorption: p.gammaEta * AL,
      absorptionRight: p.gammaEta * AR,
      A: -2 * G[1], AL: AL, AR: AR, Aeta: Aeta,
      gamma: p.gamma, gammaEta: p.gammaEta, width: p.width,
      gammaL: p.gammaL, gammaR: p.gammaR,
      channels: 1
    };
  }

  function resonantLevel(E, parameters) {
    return evaluateLevel(finite(E, "E"), levelParameters(parameters));
  }

  // Semi-infinite lead H_nn=onsite, H_n,n+1=-hopping, joined by -contact.
  // g_s=[x-sqrt(x^2-4t^2)]/(2t^2), with the retarded square-root branch.
  // Outside the band the rationalized form avoids catastrophic cancellation.
  function surfaceSelfEnergy(E, parameters) {
    finite(E, "E");
    const p = parameters || {};
    const hopping = p.hopping === undefined ? 1 : positive(p.hopping, "hopping");
    const contact = p.contact === undefined ? 1 : nonnegative(p.contact, "contact");
    const onsite = p.onsite === undefined ? 0 : finite(p.onsite, "onsite");
    const x = E - onsite;
    const bandEdge = 2 * hopping;
    const propagating = Math.abs(x) < bandEdge;
    let surface;
    if (propagating) {
      const imaginary = -Math.sqrt((bandEdge - x) * (bandEdge + x));
      surface = [x / (2 * hopping * hopping), imaginary / (2 * hopping * hopping)];
    } else {
      const radical = Math.sqrt((Math.abs(x) - bandEdge) * (Math.abs(x) + bandEdge));
      surface = [2 / (x + (x < 0 ? -radical : radical)), 0];
    }
    const re = contact * contact * surface[0];
    const im = contact * contact * surface[1];
    return {re: re, im: im, gamma: -2 * im, propagating: propagating,
      surfaceG: surface};
  }

  // Gauss-Jordan inversion with partial pivoting. These examples use small
  // dense matrices; no hidden imaginary regularizer is introduced on failure.
  function inverse(matrix) {
    const size = matrix.length;
    const augmented = matrix.map(function (row, i) {
      return row.map(function (z) { return z.slice(); }).concat(
        Array.from({length: size}, function (_, j) { return [i === j ? 1 : 0, 0]; }));
    });
    for (let column = 0; column < size; column++) {
      let pivotRow = column;
      for (let row = column + 1; row < size; row++) {
        if (abs2(augmented[row][column]) > abs2(augmented[pivotRow][column])) pivotRow = row;
      }
      if (abs2(augmented[pivotRow][column]) === 0) return null;
      const oldRow = augmented[column];
      augmented[column] = augmented[pivotRow];
      augmented[pivotRow] = oldRow;
      const pivot = augmented[column][column].slice();
      for (let j = column; j < 2 * size; j++) {
        augmented[column][j] = divide(augmented[column][j], pivot);
      }
      for (let row = 0; row < size; row++) {
        if (row === column) continue;
        const coefficient = augmented[row][column].slice();
        for (let j = column; j < 2 * size; j++) {
          const z = augmented[column][j];
          augmented[row][j][0] -= coefficient[0] * z[0] - coefficient[1] * z[1];
          augmented[row][j][1] -= coefficient[0] * z[1] + coefficient[1] * z[0];
        }
      }
    }
    return augmented.map(function (row) { return row.slice(size); });
  }

  // Rank-one leads on sites left and right. Loss is the outgoing flux into
  // Gamma_eta=2 eta I, evaluated independently from 1-R-T.
  function transportFromGreen(G, gammaL, gammaR, eta, left, right, channels) {
    const size = G.length;
    const spectral = G.map(function (row, i) { return -2 * row[i][1]; });
    const partialL = G.map(function (row) { return gammaL * abs2(row[left]); });
    const partialR = G.map(function (row) { return gammaR * abs2(row[right]); });
    const partialEta = G.map(function (row) { return 2 * eta * sum(row.map(abs2)); });
    const reflection = function (gamma, site) {
      return Math.pow(1 + gamma * G[site][site][1], 2) + Math.pow(gamma * G[site][site][0], 2);
    };
    const result = {
      G: G, T: gammaL * gammaR * abs2(G[right][left]),
      R: channels ? reflection(gammaL, left) : null,
      Rright: channels ? reflection(gammaR, right) : null,
      absorption: 2 * eta * sum(partialL),
      absorptionRight: 2 * eta * sum(partialR),
      A: sum(spectral), AL: sum(partialL), AR: sum(partialR), Aeta: sum(partialEta),
      spectral: spectral, partialL: partialL, partialR: partialR, partialEta: partialEta,
      gammaL: gammaL, gammaR: gammaR, gammaEta: 2 * eta,
      channels: channels, singular: false, sites: size
    };
    result.spectralResidual = result.A - result.AL - result.AR - result.Aeta;
    result.fluxResidual = channels ? 1 - result.R - result.T - result.absorption : null;
    return result;
  }

  function chain(E, parameters) {
    finite(E, "E");
    const p = parameters || {};
    const sites = p.sites === undefined ? 5 : integer(p.sites, "sites", 1, 64);
    const hopping = p.hopping === undefined ? 1 : positive(p.hopping, "hopping");
    const onsite = p.onsite === undefined ? 0 : finite(p.onsite, "onsite");
    const contact = p.contact === undefined ? 1 : nonnegative(p.contact, "contact");
    const eta = p.eta === undefined ? 0 : nonnegative(p.eta, "eta");
    const leadHopping = p.leadHopping === undefined ? 1 : positive(p.leadHopping, "leadHopping");
    const leadOnsite = p.leadOnsite === undefined ? 0 : finite(p.leadOnsite, "leadOnsite");
    const selfEnergy = surfaceSelfEnergy(E, {hopping: leadHopping, onsite: leadOnsite, contact: contact});
    const matrix = Array.from({length: sites}, function (_, i) {
      return Array.from({length: sites}, function (_, j) {
        return [i === j ? E - onsite : (Math.abs(i - j) === 1 ? hopping : 0), i === j ? eta : 0];
      });
    });
    // For a one-site device both self-energies act on the same diagonal.
    [0, sites - 1].forEach(function (site) {
      matrix[site][site][0] -= selfEnergy.re;
      matrix[site][site][1] -= selfEnergy.im;
    });
    const G = inverse(matrix);
    const channels = selfEnergy.propagating ? 1 : 0;
    if (!G) {
      if (eta > 0 || selfEnergy.gamma > 0) {
        throw new RangeError("Numerically singular open-device resolvent");
      }
      // An isolated pole or a threshold singularity is a distributional limit.
      // Scattering is still reflection-only for a disconnected open lead.
      return {T: 0, R: channels ? 1 : null, Rright: channels ? 1 : null,
        absorption: 0, absorptionRight: 0, A: null, AL: 0, AR: 0, Aeta: 0,
        spectral: Array(sites).fill(null), partialL: Array(sites).fill(0),
        partialR: Array(sites).fill(0), partialEta: Array(sites).fill(0), G: null,
        gammaL: 0, gammaR: 0, gammaEta: 0, channels: channels, sites: sites,
        singular: true, selfEnergy: selfEnergy,
        spectralResidual: null, fluxResidual: channels ? 0 : null};
    }
    const result = transportFromGreen(G, selfEnergy.gamma, selfEnergy.gamma, eta, 0, sites - 1, channels);
    result.selfEnergy = selfEnergy;
    return result;
  }

  function fano(E, parameters) {
    finite(E, "E");
    const input = parameters || {};
    const p = levelParameters(input);
    const sideEnergy = input.sideEnergy === undefined ? 0 : finite(input.sideEnergy, "sideEnergy");
    const coupling = input.coupling === undefined ? 0.5 : nonnegative(input.coupling, "coupling");
    if (coupling === 0) {
      const central = evaluateLevel(E, p);
      const atBoundPole = p.eta === 0 && E === sideEnergy;
      const sideG = atBoundPole ? null : divide([1, 0], [E - sideEnergy, p.eta]);
      const sideA = atBoundPole ? null : -2 * sideG[1];
      const result = Object.assign({}, central, {
        G: [[central.G, [0, 0]], [[0, 0], sideG]],
        A: sideA === null ? null : central.A + sideA,
        Aeta: sideA === null ? null : central.Aeta + sideA,
        spectral: [central.A, sideA], partialL: [central.AL, 0],
        partialR: [central.AR, 0], partialEta: [central.Aeta, sideA],
        sites: 2, singular: atBoundPole,
        boundStates: p.eta === 0 ? [{energy: sideEnergy, weight: 1, site: 1}] : [],
        spectralIncludesBoundStates: p.eta > 0
      });
      result.fluxResidual = 1 - result.R - result.T - result.absorption;
      result.spectralResidual = sideA === null ? null : result.A - result.AL - result.AR - result.Aeta;
      return result;
    }
    const G = inverse([
      [[E - p.epsilon, p.width / 2], [coupling, 0]],
      [[coupling, 0], [E - sideEnergy, p.eta]]
    ]);
    if (!G) throw new RangeError("An isolated two-level pole requires eta > 0 to sample its spectral weight");
    const result = transportFromGreen(G, p.gammaL, p.gammaR, p.eta, 0, 0, 1);
    result.boundStates = [];
    result.spectralIncludesBoundStates = true;
    return result;
  }

  function fermi(E, mu, kT) {
    finite(E, "E"); finite(mu, "mu"); nonnegative(kT, "kT");
    if (kT === 0) return E < mu ? 1 : (E > mu ? 0 : 0.5);
    const x = (E - mu) / kT;
    // These forms remain finite even at |x| much larger than 700.
    if (x >= 0) { const tail = Math.exp(-x); return tail / (1 + tail); }
    return 1 / (1 + Math.exp(x));
  }

  function integrate(fn, lower, upper, intervals) {
    if (typeof fn !== "function") throw new TypeError("fn must be a function");
    finite(lower, "lower"); finite(upper, "upper");
    integer(intervals, "intervals", 1, 10000000);
    if (upper < lower) throw new RangeError("upper must be at least lower");
    if (upper === lower) return 0;
    const step = (upper - lower) / intervals;
    let total = 0;
    for (let i = 0; i <= intervals; i++) {
      const E = i === intervals ? upper : lower + i * step;
      const value = finite(fn(E), "integrand");
      total += (i === 0 || i === intervals ? 0.5 : 1) * value;
    }
    return total * step;
  }

  function trapz(x, y) {
    if (!x || !y || x.length !== y.length || x.length < 2) {
      throw new RangeError("trapz requires x and y arrays of equal length, with at least two samples");
    }
    let total = 0;
    for (let i = 0; i < x.length; i++) {
      finite(x[i], "x[" + i + "]"); finite(y[i], "y[" + i + "]");
      if (i > 0) {
        if (x[i] <= x[i - 1]) throw new RangeError("x samples must be strictly increasing");
        total += (x[i] - x[i - 1]) * (y[i] + y[i - 1]) / 2;
      }
    }
    return total;
  }

  function integrationParameters(input, defaultLimit, defaultIntervals) {
    const p = input || {};
    const limit = p.limit === undefined ? defaultLimit : positive(p.limit, "limit");
    const intervals = p.intervals === undefined ? defaultIntervals : integer(p.intervals, "intervals", 1, 10000000);
    const offset = p.offset === undefined ? 0 : finite(p.offset, "offset");
    return {limit: limit, intervals: intervals, lower: -limit + offset,
      upper: limit + offset, gridStep: 2 * limit / intervals};
  }

  function etaDiagnostics(parameters) {
    const p = levelParameters(parameters);
    positive(p.width, "total spectral width");
    const grid = integrationParameters(parameters, 8, 1000);
    const area = integrate(function (E) { return evaluateLevel(E, p).A / (2 * Math.PI); },
      grid.lower, grid.upper, grid.intervals);
    const halfWidth = p.width / 2;
    const analyticArea = (Math.atan((grid.upper - p.epsilon) / halfWidth)
      - Math.atan((grid.lower - p.epsilon) / halfWidth)) / Math.PI;
    const peak = evaluateLevel(p.epsilon, p);
    return Object.assign({}, grid, {
      area: area, analyticArea: analyticArea, infiniteArea: 1,
      quadratureError: area - analyticArea, omittedTail: 1 - analyticArea,
      peakT: peak.T, peakR: peak.R, peakLoss: peak.absorption,
      width: p.width, peak: peak, samples: grid.intervals + 1
    });
  }

  // Fixed level and fixed couplings under bias. For eta>0 the deliberately
  // EMPTY sink has no injection. I is the L/R pair current; it is not then
  // the full left contact current. IL and IR point into the device, Isink out.
  // A finite window and unresolved peaks may spoil n<=1; values are not clipped.
  function bias(parameters) {
    const input = parameters || {};
    const p = levelParameters(input);
    positive(p.width, "total spectral width");
    const voltage = input.bias === undefined ? 1 : finite(input.bias, "bias");
    const kT = input.kT === undefined ? 0.025 : nonnegative(input.kT, "kT");
    const grid = integrationParameters(input, 12, 2400);
    const muL = voltage / 2, muR = -voltage / 2;
    let pair = 0, fromL = 0, fromR = 0, spectral = 0;
    // At zero temperature the Fermi steps are discontinuities, not narrow
    // smooth curves. Insert their exact locations and integrate each side
    // using its one-sided occupation. Assigning f(mu)=1/2 at uniform-grid
    // samples alone gives an O(dE) error for off-grid or rounded endpoints.
    const breakpoints = kT === 0 ? Array.from(new Set([muL, muR])).filter(function (mu) {
      return mu > grid.lower && mu < grid.upper;
    }).sort(function (a, b) { return a - b; }) : [];
    let breakpointIndex = 0, actualIntervals = 0;
    let leftEnergy = grid.lower;
    let leftLevel = evaluateLevel(leftEnergy, p);
    function addInterval(rightEnergy) {
      const rightLevel = evaluateLevel(rightEnergy, p);
      const halfStep = (rightEnergy - leftEnergy) / 2;
      if (kT === 0) {
        // No chemical potential lies strictly inside this interval.
        const fL = leftEnergy < muL ? 1 : 0;
        const fR = leftEnergy < muR ? 1 : 0;
        pair += halfStep * (leftLevel.T + rightLevel.T) * (fL - fR);
        fromL += halfStep * (leftLevel.AL + rightLevel.AL) * fL;
        fromR += halfStep * (leftLevel.AR + rightLevel.AR) * fR;
      } else {
        const fLeftL = fermi(leftEnergy, muL, kT), fLeftR = fermi(leftEnergy, muR, kT);
        const fRightL = fermi(rightEnergy, muL, kT), fRightR = fermi(rightEnergy, muR, kT);
        pair += halfStep * (leftLevel.T * (fLeftL - fLeftR) + rightLevel.T * (fRightL - fRightR));
        fromL += halfStep * (leftLevel.AL * fLeftL + rightLevel.AL * fRightL);
        fromR += halfStep * (leftLevel.AR * fLeftR + rightLevel.AR * fRightR);
      }
      spectral += halfStep * (leftLevel.A + rightLevel.A);
      leftEnergy = rightEnergy;
      leftLevel = rightLevel;
      actualIntervals++;
    }
    for (let i = 1; i <= grid.intervals; i++) {
      const rightEnergy = i === grid.intervals ? grid.upper : grid.lower + i * grid.gridStep;
      while (breakpointIndex < breakpoints.length && breakpoints[breakpointIndex] < rightEnergy) {
        const edge = breakpoints[breakpointIndex++];
        if (edge > leftEnergy) addInterval(edge);
      }
      addInterval(rightEnergy);
    }
    const factor = constants.currentMicroampPerEV;
    const energyIntegral = pair;
    const injectionL = fromL;
    const injectionR = fromR;
    const I = factor * energyIntegral;
    const sinkL = factor * p.gammaEta * injectionL;
    const sinkR = factor * p.gammaEta * injectionR;
    const IL = I + sinkL, IR = -I + sinkR, Isink = sinkL + sinkR;
    return Object.assign({}, grid, {
      I: I, IL: IL, IR: IR, Isink: Isink,
      n: (injectionL + injectionR) / (2 * Math.PI),
      nL: injectionL / (2 * Math.PI), nR: injectionR / (2 * Math.PI),
      spectralArea: spectral / (2 * Math.PI),
      energyIntegral: energyIntegral, muL: muL, muR: muR, kT: kT,
      width: p.width, gammaEta: p.gammaEta, conservationResidual: IL + IR - Isink,
      currentUnit: "microampere", occupationConvention: "empty eta sink",
      samples: actualIntervals + 1, actualIntervals: actualIntervals,
      insertedBreakpoints: actualIntervals - grid.intervals,
      integrationScheme: kT === 0 ? "trapezoids split at chemical potentials" : "uniform trapezoids"
    });
  }

  return Object.freeze({constants: constants, resonantLevel: resonantLevel,
    surfaceSelfEnergy: surfaceSelfEnergy, chain: chain, fano: fano,
    fermi: fermi, integrate: integrate, trapz: trapz, bias: bias,
    etaDiagnostics: etaDiagnostics});
});
