# Metrics

## Kerr-Schild metric

The Kerr-Schild metric takes the form:

```math
    g_{\mu \nu} = \eta_{\mu \nu} + f \, k_\mu \, k_\nu 
```

where 
```math
  k_\mu = \left( \frac{r \, x + a \, y}{r^2+a^2}
                , \frac{r \, y - a \, x}{r^2+a^2} \frac{z}{r}\right)
```
```math
  f = \frac{2 \, G \, M \, r^3}{r^4+a^2 \, z^2},
```

and ``r`` is implicitly defined by:

```math
\frac{x^2+y^2}{r^2+a^2} + \frac{z^2}{r^2} = 1 .
```

## Weak field metric

For solar system and terrestrial positioning, the weak field metric suffices. The weak field metric has the form:

```math
g_{\mu \nu} = \eta_{\mu \nu} - 2 \, V \, \delta_{\mu \nu},
```

where ``V`` is the gravitational potential.

## Gordon metric

To incorporate atmospheric and ionospheric effects, one uses the analogue Gordon metric, which takes the form (with ``g_{\mu \nu}`` being the gravitational metric):

```math
    \bar{g}_{\mu \nu} = g_{\mu \nu} + \left(1-\frac{1}{n^2}\right) u_\mu u_\nu ,
```
