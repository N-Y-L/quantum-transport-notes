# MathJax SVG runtime

This folder contains the unmodified `es5/tex-svg.js` bundle from **MathJax 3.2.2**,
distributed under the [Apache License 2.0](LICENSE). The upstream README is
preserved as [README.upstream.md](README.upstream.md). This third-party license
is separate from the notes' MIT license.

Source: [official npm package archive](https://registry.npmjs.org/mathjax/-/mathjax-3.2.2.tgz).
The archive was checked against the `dist.integrity` value in the
[official version metadata](https://registry.npmjs.org/mathjax/3.2.2):

```text
sha512-Bt+SSVU8eBG27zChVewOicYs7Xsdt40qm4+UpHyX7k0/O9NliPc+x77k1/FEsPsjKPZGJvtRZM1vO+geW0OhGw==
```

The bundled `tex-svg.js` has SHA-256:

```text
d4295dc33744836935c1399feece5159577b34c5c8ffb9f1c6324cd82e03a882
```

The notes load `svg-math.js` before this runtime and await
`TransportSvgMath.ready` before drawing figures. MathJax converts individual TeX
labels to SVG paths. Document equations continue to use KaTeX. Automatic page
typesetting is disabled; only the bundled `base`, `ams`, and `newcommand` TeX
packages are enabled. Labels do not load extensions or external fonts.

`TransportSvgMath.render(tex, x, baselineY, options)` returns an SVG group for
insertion into an existing SVG figure. Options are `fontSize` (default 16),
`anchor` (`start`, `middle`, or `end`), `color` (an explicit CSS color name or hex
value), and `rotate` (degrees about the anchor and baseline). All coordinates
use the enclosing SVG's units. The group retains the source as `data-tex` and
an accessible title. Paths are inline, with no shared font cache, image, font,
or document references. The helper caches at most 256 TeX templates.

See the official [conversion API](https://docs.mathjax.org/en/v3.2/web/typeset.html#converting-a-math-string-to-other-formats)
and [SVG output options](https://docs.mathjax.org/en/v3.2/options/output/svg.html).
