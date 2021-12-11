# Squirrel algorithm

## Description

The squirrel algorithm takes the coordinates of `ne ≥ 4` emission
points, arranged as column vectors in a `4×ne` matrix `X`,  and computes
the intersection of future pointing null geodesics defined with respect
to a (slightly curved) spacetime metric `g`. 

Null geodesics may be described as solutions to the geodesic equation
(implemented here in the form of Hamilton's equations), which may be
formally written as the functions:

```math
    x^μ_I(λ,{X}_I,{\bf v}_I),
```

where the indices ``I∈\{1,2,...,n_e\}`` distinguish the emission
points ``X_I`` and their associated null geodesics, and ``{\bf v}_I``
denotes the spatial components of the initial four-velocity vector for
the null geodesic (the time component is determined by the requirement
that the four-velocity is null). 

Given a collection of four geodesic functions ``\{x_1,x_2,x_3,x_4\}``
(each having the form ``x^μ_I=x^μ_I(λ,{x}_0,{\bf v})``), the condition
that they intersect is the vanishing of the following vector-valued
function:

```math
f := \left( x_1 - x_2 , x_1 - x_3 , x_1 - x_4 \right)
```

which may be formally written as ``f=f(X_1,X_2,X_3,X_4,{\bf v}_1,{\bf
v}_2,{\bf v}_3,{\bf v}_4)``. A variable counting exercise reveals that
``f`` has ``12`` components; since there are a total of 12 undetermined
quantities in the four initial velocities ``{\bf v}_I`` (each of which
have three components), the condition ``f=0`` can thought of as a set of
``12`` equations for the ``12`` unknowns ``v_I``. The squirrel algorithm
seeks to find the roots of the function ``f``.

The squirrel algorithm is then summarized:

1.  Apply the flat spacetime algorithm to the emission points `X` to
    obtain an initial guess for the zeros of the function `f`.

2.  Apply a root finding algorithm to the function `f` to obtain the
    initial velocities ``v_I``.

3.  Integrate the geodesics with the resulting initial velocities
    ``v_I`` and emission points ``X`` to find the intersection point.

A quasi-Newton Broyden algorithm (which will be described in detail
below) is employed to do the root-finding; in such a method, the
Jacobian for `f` is computed once in the first iteration of the root
finding algorithm, and is updated in the subsequent iterations. The
function `f` is computed by way of numerical integration of geodesics;
if the numerical integration is performed using native Julia libraries,
one can compute the Jacobian by way of automatic differentiation.

## Geodesic endpoint function

To compute the function `f`, the geodesic endpoint function
``x^μ_I=x^μ_I(λ,X_I,{\bf v}_I)|_{λ=1}`` is implemented using the native
Julia ODE solvers in the library
[`OrdinaryDiffEq.jl`](https://github.com/SciML/OrdinaryDiffEq.jl), using
the recommended method `AutoVern7(Rodas5())`

```@docs
squirrel.gsolve
```

With the geodesic endpoint functions ``x^μ_I=x^μ_I(λ,X_I,{\bf
v}_I)|_{λ=0}`` in hand, one can construct the function ``f=f(X_I,{\bf
v}_I)``:

```@docs
squirrel.zF
```

## Jacobian

Next, one computes the Jacobian of `f`. As mentioned earlier, this is
done by way of automatic differentiation, using the library
[`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl). Here,
the Jacobian matrix of ``x^μ_I=x^μ_I(λ,X_I,{\bf v}_I)|_{λ=1}`` (which one may write schematically as ``\frac{∂x_I}{∂v_J}``) is
computed:

```@docs
squirrel.gejac
```

Given ``\frac{∂x_I}{∂v_J}``, the Jacobian matrix of the function `f` may
be computed by way of the chain rule, as indicated in the schematic
formula:

```math
J = \frac{∂f}{∂x_I}\frac{∂x_I}{∂v_J},
```

The following function applies the chain rule to compute the Jacobian
matrix:

```@docs
squirrel.geocJ
```

## Root finding algorithm

Given some function ``f(x)``, the standard approach is the Newton
method, which finds the roots according to the prescription:

```math
    x_{i+1} = x_i + {\bf J}^{-1}_{i} f(x_i)
```

where ``{\bf J}^{-1}_{i}`` is the inverse of the Jacobian matrix ``{\bf
J}`` for ``f(x)`` evaluated at ``x_i``. However, in situations which the
computation of the Jacobian matrix ``{\bf J}`` becomes expensive, one
may instead employ the Broyden method, which is a quasi-Newton root
finding method. The Broyden algorithm involves first computing the
Jacobian matrix ``{\bf J}`` and its inverse for ``f(x)``. However,
instead of computing the Jacobian matrix at each iteration (as is done
in the Newton method), the Broyden algorithm updates the (inverse)
Jacobian matrix according to the Sherman-Morrison update formula:

```math
    {\bf J}^{-1}_{i+1} 
    = 
        {\bf J}^{-1}_{i}
        +
        \frac{Δv^{T}_i - {\bf J}^{-1}_{i} Δf_i }
        {Δv^{T}_i {\bf J}^{-1}_{i} Δf_i }
        Δv^{T}_i {\bf J}^{-1}_{i},
```

which is implemented in the following function:

```@docs
squirrel.JiSMU
```

The Broyden update formula is implemented in the following function:

```@docs
squirrel.bsolve
```

## Initial data finder

The following function makes use of the preceding functions to construct
the initial data for the four-velocities:

```@docs
squirrel.idf
```

## Locator (Four points)

Finally, the `locator4` function computes the intersection point from a
set of four emission points `X` by first employing the `idf` function to
obtain the initial data for the geodesics, and integrates the geodesics
(up to unit affine parameter ``λ``) to obtain the intersection point:

```@docs
squirrel.locator4
```

## Locator (``n_e>4`` emission points)

If more than four emission points are available, the following function
can generate ``C(n_e,4)`` sets of emission points (where ``C(n,k)`` is
the binomial coefficient):

```@docs
squirrel.combX
```

One feeds the output of `combX` into the `locator4` function to obtain
``C(n_e,4)`` solutions. The following outlier detection formula, which
is based on the comparison with median values, can then be used to
exclude large errors (the assumption here is that the large errors are
infreqent):

```@docs
squirrel.odetc
```

The following function implements the procedure described above, given a matrix of ``4×n_e`` emission points ``X``:

```@docs
squirrel.locator
```
