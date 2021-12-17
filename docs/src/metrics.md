# Metrics

## Kerr-Schild metric

The Kerr-Schild metric takes the form:

```math
    g_{\mu \nu} = \eta_{\mu \nu} + f \, k_\mu \, k_\nu 
```

where 
```math
  k_\mu = \left( 1 , \frac{r \, x + a \, y}{r^2+a^2}
                , \frac{r \, y - a \, x}{r^2+a^2} , \frac{z}{r}\right)
```
```math
  f = \frac{2 \, G \, M \, r^3}{r^4+a^2 \, z^2},
```

and ``r`` is implicitly defined by:

```math
\frac{x^2+y^2}{r^2+a^2} + \frac{z^2}{r^2} = 1 .
```

The components of the Kerr-Schild metric can be calculated using
`gks` function in general

```@docs	
squirrel.metric.gks
```

or using the `ge` function for Earth-like parameters

```@docs	
squirrel.metric.ge
```

## Weak field metric

For solar system and terrestrial positioning, the weak field metric
suffices. The weak field metric has the form:
```math
g_{\mu \nu} = \eta_{\mu \nu} - 2 \, V \, \delta_{\mu \nu},
```

where ``V`` is the gravitational potential.

## Gordon metric

To incorporate atmospheric and ionospheric effects, one uses the
analogue Gordon metric, which takes the form (with ``g_{\mu \nu}`` being
the gravitational metric):
```math
    \bar{g}_{\mu \nu} = g_{\mu \nu} + 
    \left(1-\frac{1}{n^2}\right) u_\mu u_\nu ,
```

The general Gordon metric is implemented in the following function:
```@docs	
squirrel.metric.gGordon
```

In principle, the function implemented above can be used with any user
defined index of refraction profile. However, for evaluation purposes,
we implement some simplified index of refraction profiles in the
following sections.

### Index of refraction models

The index of refraction is written in the following form:
```math
n = 1 + \Delta n_{\rm atm} + \Delta n_{\rm ion}
```

#### Atmospheric component

The atmospheric component for the index of refraction is given here as a
function of the height ``h`` from the surface of the Earth:

```math
\Delta n_{atm} = \frac{A_1}{B_1 + C_1 \left(h-H_1)\right)}
  + \frac{A_2}{B_2 + C_2 \left(h-H_2)\right)}
```
with the parameter choices:
```math
\begin{aligned}
    A_1 &= -222.666    &\qquad   A_2 &= -253.499  \\
    B_1 &= 99.0621     &\qquad   B_2 &= 112.757  \\
    C_1 &= 0.157823 ~\text{km}^{-1}
        &\qquad C_2 &= 0.179669 ~\text{km}^{-1}\\
    H_1 &= -7.1541  ~\text{km}
        &\qquad   H_2 &= -7.15654 ~\text{km}
\end{aligned}
```

The function ``\Delta n_{atm}`` is implemented in the following:

```@docs	
squirrel.metric.ΔnatmStd
```

#### Ionospheric component

The ionospheric component index of refraction is given by the following:
```math
\begin{aligned}
    \Delta n_{ion} \approx 
              \left(4.024 \times 10^{-17}\right) [N_e/\text{m}^{-3}]
\end{aligned}
```

where the electron density is given by:
```math
\begin{aligned}
    N_e :=& \biggl[\alpha_D \, \tilde Ep(h,h_D,b_D)
          + \alpha_E \, \tilde Ep(h,h_E,b_E) 
          + \alpha_F \, \tilde Ep(h,h_F,b_F)\biggr],
\end{aligned}
```

with the parameter choices:
```math
\begin{aligned}
    \alpha_D &= 10^{12}~\text{m}^{-3}
                   &\quad   h_D &=  75~\text{km}
                   &\quad   b_D &=   5~\text{km} \\
    \alpha_E &= 2.5 \times 10^{11}~\text{m}^{-3}
                   &\quad   h_E &= 130~\text{km}
                   &\quad   b_E &=  30~\text{km} \\
    \alpha_F &= 10^{11}~\text{m}^{-3}
                   &\quad   h_F &= 300~\text{km}
                   &\quad   b_F &=  50~\text{km} .
\end{aligned}
```

The function `` \tilde Ep(h,h_E,b_E)`` is the pseudo-Epstein function:
```math
\begin{aligned}
    \tilde Ep(h,h_c,B) := &\frac{1}{16}
    \biggl\{
      \left[1+\left(\tfrac{h-h_c}{2 B}\right)^2\right]^{-1}
      +\left[1+\left(\tfrac{h-h_c}{4 B}\right)^2\right]^{-2}
      +\left[1+\left(\tfrac{h-h_c}{6 B}\right)^2\right]^{-3}
      +\left[1+\left(\tfrac{h-h_c}{7 B}\right)^2\right]^{-4}
      \biggr\}^2,
\end{aligned}
```

defined as:

```@docs	
squirrel.metric.pEp
```

The ionospheric index of refraction, also written in terms of the height
``h`` from the surface of the Earth, is implemented in the following: 
```@docs	
squirrel.metric.Δnios
```

The unperturbed total index of refraction ``\Delta n_{\rm tot} = \Delta
n_{\rm atm} + \Delta n_{\rm ion}`` is given in the following (in
surface-adapted coordinates):
```@docs	
squirrel.metric.Δntot
```

#### Perturbation models

To model uncertainties in the index of refraction profile, perturbations
of the following form are introduced (also in surface-adapted
coordinates):
```math
\begin{aligned}
    \Delta n_{\rm pert} &= \Delta n_{\rm atm} \left(1
                 + \delta_1 \, \tilde p_1(h) \right) \\
    & \qquad 
    + \Delta n_{\rm ion} \left(1 + \delta_2 \, \tilde  p_2(h)\right),
\end{aligned}
```
and implemented in the following.
```@docs	
squirrel.metric.Δntp
```

The respective perturbation profile function ``p_A`` and line shape
function ``Ls`` are defined:
```math
\tilde p_A(h) = \sum_{i} \alpha_i \tilde Ls(h,h_{0,i},\sigma_i)
```
```@docs	
squirrel.metric.P
```
```math
Ls(h,h_{0},\sigma) =  \frac{\sigma^2}{\sigma^2 
                      + (h-h_0)^2} \frac{\sigma^4}
                        {\sigma^4 + (h-h_0)^4}
```
```@docs	
squirrel.metric.LSF
```

The parameters ``h_i`` and  ``\sigma_i`` are fed into the perturbation
profile functions as vectors.

#### Conversion functions

Index of refraction functions are converted from surface-adapted
coordinates to spherical coordinates with the following function:
```@docs
squirrel.metric.nIRs
```

Finally, the index of refraction functions are converted to Cartesian
coordinates with the following function:
```@docs
squirrel.metric.nIR
```
