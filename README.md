# Quantum transport: visual notes

Five compact experiments follow an electron from an incoming scattering wave to
a measured current. The reader needs basic graduate quantum mechanics and some
condensed matter, but no quantum field theory or Green-function transport
background. The pages use ordinary serif text, locally rendered LaTeX, and
adjustable SVG figures.

## Open the notes

**Double-click `index.html`.** No installation, server, account, build step, or
internet connection is needed for the calculations. External references need a
connection. Equations use the bundled KaTeX distribution.

| Experiment | Try this |
|---|---|
| One level, two doors | Weaken both contacts equally: a resonance narrows without losing its unit peak. |
| Let the leads continue | Change the length of a matched chain, then weaken its contacts. |
| A detour that stops traffic | Move a side orbital and follow the transmission zero. |
| Who fills the level? | Compare zero current at equilibrium with unequal reservoir filling under bias. |
| The imaginary-number audit | Reduce η, shift the integration grid, and separate quadrature error from missing spectral tails. |

Read the [connected derivation](docs/notes.html) for the retarded resolvent,
lead elimination, Γ, the Caroli trace, occupations, and lifetime conventions.
[Annotated references](docs/references.html) identify primary papers and Datta's
teaching. Each experiment has reset and example buttons. “Save figure” exports
the plotted SVG; “Save calculation” exports model parameters, units, results,
and curve samples as JSON. A persistent download link permits another attempt.
The export source is also available in an expandable text panel for inspection or copying.
These exports are records, not importable sessions.

## Models and conventions

All calculations use a static, noninteracting Hamiltonian in an orthonormal
orbital basis, one spin species, energies in **eV**, currents in **μA**, and
time in **fs**. A second identical independent spin species doubles the current.
Positive displayed current counts positive elementary charge times electron
transfer from left to right; signed electric charge flow in that direction has
the opposite sign.

- The resonant level and side orbital use constant wide-band lead self-energies.
- The finite chain is embedded in two exact semi-infinite nearest-neighbor
  chains. At a band threshold there is no nonzero-velocity incident channel;
  the displayed threshold transmission is zero by convention.
- The bias experiment fixes the Hamiltonian and contact spectra while changing
  reservoir chemical potentials. It does not solve electrostatics, Coulomb
  interactions, or inelastic scattering. Its occupation is a finite-window
  integral; the omitted spectral tail is displayed as a bound on omitted
  occupation.
- The adjustable η occurs only in the **device** resolvent. Physical lead
  self-energies use the retarded real-axis limit. Retaining finite η without
  reinjection is mathematically an empty absorbing reservoir with width 2η.
  The current audit includes this loss explicitly. It is not a conserving
  dephasing model.

The probability survival time ℏ/Γ and amplitude e-folding time 2ℏ/Γ apply to a
single wide-band pole with physical width Γ. A numerical retarded infinitesimal
does not establish a material lifetime. Bound-state spectral weight is
distributional at η = 0; isolated occupations require additional preparation
information.

## Reproduce the numerical checks

Node.js 22 or newer is sufficient; **no package installation is needed**:

```sh
node --test tests/physics.test.js
node scripts/check-numerics.js
```

The second command prints the report. To regenerate the included
[numerical report](docs/numerical-checks.json), run:

```sh
node scripts/check-numerics.js --write
```

The report records model parameters, units, integration windows, grid
alignment, runtime version, and exact reference formulas. The calculations are
deterministic and use no random sampling. Independent checks cover:

- Lorentzian line shapes, peak heights, widths, and lifetime factors.
- Causal self-energy signs, the surface quadratic, and the far-band limit.
- Chain scattering against a separate discrete Schrödinger wave-matching
  calculation; one-channel bounds and flux balance.
- Fano zeros, the disconnected side orbital, and the two-orbital spectral sum.
- Zero-temperature current and occupation against exact arctangent integrals;
  small-bias conductance and finite-temperature current balance.
- Finite η spectral identities, explicit absorbed flux, finite-window tails,
  grid-offset aliasing, and resolution of narrow plotted resonances.

The spectral integration is composite trapezoidal quadrature. At zero
temperature, current and occupation integrals also split at the chemical
potentials, using the one-sided occupations on each segment. This avoids
grid-dependent Fermi-step errors. The bias interface compares two grid
spacings. The η lab intentionally permits unresolved integrals and reports
their errors instead of clipping the answers.

Plot samples and integration samples serve different purposes. The chain
curve resolves resonances with adaptive sampling; the Fano curve adds samples
on the antiresonance width scale. The η curve resolves the analytic peak for
display while its dots and table expose the actual uniform integration grid.
The finite-η sink-current readout uses exact zero-temperature integrals over
the infinite energy axis, separately from that intentionally fallible grid.

## File guide

| File | Purpose |
|---|---|
| `index.html` | Reading map |
| `experiments.html` | Physical questions, equations, and controls |
| `physics.js` | Numerical models, independent of the interface |
| `sampling.js` | Resonance-aware curve sampling |
| `app.js` | Controls, SVG drawings, readouts, and exports |
| `docs/notes.html`, `docs/references.html` | Connected derivations and sources |
| `tests/physics.test.js` | Independent physical and numerical checks |
| `scripts/check-numerics.js` | Reproducible numerical report |

Edit the HTML directly; there is no documentation build pipeline.
For an optional local preview with Python 3:

```sh
python3 -m http.server 8000 --bind 127.0.0.1
```

Open [localhost:8000](http://127.0.0.1:8000). On Windows use `py -3` in place
of `python3`. Stop the server with Ctrl+C.

## License

Original code, text, and figures: [MIT](LICENSE).
Copyright (c) 2026 Neil Yuanting Li.
Bundled KaTeX 0.19.0 retains its [third-party MIT license](vendor/katex/LICENSE).
