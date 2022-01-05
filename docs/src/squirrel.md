# Squirrel algorithm

## Description

The squirrel algorithm takes the coordinates of ``n_{\rm e} ≥ 4``
emission points, arranged as column vectors in a ``4×n_{\rm e}`` matrix
``X``,  and computes the intersection of future pointing null geodesics
defined with respect to a (slightly curved) spacetime metric `g`. 

Null geodesics may be described as solutions to the geodesic equation
(implemented here in the form of Hamilton's equations), which may be
formally written as the functions:

```math
    x^μ_I(λ,{X}_I,{\bf v}_I),
```

where the indices ``I∈\{1,2,...,n_{\rm e}\}`` distinguish the emission
points ``X_I`` and their associated null geodesics, and ``{\bf v}_I``
denotes the spatial components of the initial four-velocity vector for
the null geodesic (the time component is determined by the requirement
that the four-velocity is null). 

Given a collection of four geodesic functions ``\{x_1,x_2,x_3,x_4\}``
(each having the form ``x^μ_I=x^μ_I(λ,X_I,{\bf v})``) with endpoints
given by ``{\rm x}^μ_I:=x^μ_I(1,X_I,{\bf v}_I)``, the condition
that they intersect is the vanishing of the following vector-valued
function:

```math
F := \left( {\rm x}_1 - {\rm x}_2 , {\rm x}_1 - {\rm x}_3 
            , {\rm x}_1 - {\rm x}_4 \right),
```

which may be formally written as 
``F=F(X_1,X_2,X_3,X_4,{\bf v}_1,{\bf v}_2,{\bf v}_3,{\bf v}_4)``. 
A variable counting exercise reveals that
``F`` has ``12`` components; since there are a total of 12 undetermined
quantities in the four initial velocities ``{\bf v}_I`` (each of which
have three components), the condition ``F=0`` can thought of as a set of
``12`` equations for the ``12`` unknowns ``{\bf v}_I``. Since the relevant 
variables for the root finding algorithm are ``{\bf v}_I``, one may 
suppress the arguments ``X_I`` to write 
``f(v)=F(X_1,X_2,X_3,X_4,{\bf v}_1,{\bf v}_2,{\bf v}_3,{\bf v}_4)``, where 
``v:=({\bf v}_1,{\bf v}_2,{\bf v}_3,{\bf v}_4)`` represents the 
concatenation of the vectors ``{\bf v}_I``. The squirrel algorithm
seeks to find the roots of the function ``f(v)``.

The squirrel algorithm is then summarized:

1.  Apply the flat spacetime algorithm to the emission points ``X`` to
    obtain an initial guess for the zeros of the function ``f(v)``.

2.  Apply a root finding algorithm to the function ``f(v)`` to obtain the
    initial velocities ``{\bf v}_I``.

3.  Integrate the geodesics with the resulting initial velocities
    ``{\bf v}_I`` and emission points ``X`` to find the intersection point.

A quasi-Newton Broyden algorithm (which will be described in detail
below) is employed to do the root-finding; in such a method, the
Jacobian for ``f(v)`` is computed once in the first iteration of the root
finding algorithm, and is updated in the subsequent iterations. The
function ``f(v)`` is computed by way of numerical integration of geodesics;
if the numerical integration is performed using native Julia libraries,
one can compute the Jacobian by way of automatic differentiation.

## Geodesic endpoint function

To compute the function ``f``, the geodesic endpoint function 
``{\rm x}^μ_I=x^μ_I(λ,X_I,{\bf v}_I)|_{λ=1}`` is implemented using the 
native Julia ODE solvers in the library
[`OrdinaryDiffEq.jl`](https://github.com/SciML/OrdinaryDiffEq.jl), using
the recommended method `AutoVern7(Rodas5())`

```@docs
squirrel.gsolve
```

With the geodesic endpoint functions ``{\rm x}^μ_I=x^μ_I(λ,X_I,{\bf
v}_I)|_{λ=0}`` in hand, one can construct the function 
``f(v)=F(X_1,X_2,X_3,X_4,{\bf v}_1,{\bf v}_2,{\bf v}_3,{\bf v}_4)``:

```@docs
squirrel.zF
```

## Jacobian

Next, one computes the Jacobian of ``f(v)``. As mentioned earlier, this 
is done by way of automatic differentiation, using the library
[`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl). Here,
the Jacobian matrix of ``{\rm x}^μ_I=x^μ_I(λ,X_I,{\bf v}_I)|_{λ=1}`` 
(which one may write schematically as ``\frac{∂{\rm x}_I}{∂v_A}``) is
computed:

```@docs
squirrel.gejac
```

Given ``\frac{∂{\rm x}_I}{∂v_A}`` (note ``A ∈\{1,2,...,12\}``, since
``v`` is a concatenation of the vectors ``{\bf v}_I``), the Jacobian 
matrix of the function ``f(v)`` may be computed by way of the chain rule, 
as indicated in the schematic formula (the Jacobian matrix ``{\bf J}`` is 
not to be confused with the symplectic matrix ``J^{\alpha \beta}``):

```math
{\bf J} = \frac{∂f}{∂{\rm x}_I}\frac{∂{\rm x}_I}{∂v_A},
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
    x_{{\rm i}+1} = x_{\rm i} + {\bf J}^{-1}_{\rm i} f(x_{\rm i})
```

where ``{\bf J}^{-1}_{\rm i}`` is the inverse of the Jacobian matrix
``{\bf J}`` for ``f(x)`` evaluated at ``x_{\rm i}``. However, in
situations which the computation of the Jacobian matrix ``{\bf J}``
becomes expensive, one may instead employ the Broyden method, which is a
quasi-Newton root finding method. The Broyden algorithm involves first
computing the Jacobian matrix ``{\bf J}`` and its inverse for ``f(x)``.
However, instead of computing the Jacobian matrix at each iteration (as
is done in the Newton method), the Broyden algorithm updates the
(inverse) Jacobian matrix according to the Sherman-Morrison update
formula:

```math
    {\bf J}^{-1}_{{\rm i}+1} 
    = 
        {\bf J}^{-1}_{\rm i}
        +
        \frac{Δx^{T}_{\rm i} - {\bf J}^{-1}_{\rm i} Δf_{\rm i} }
        {Δx^{T}_{\rm i} {\bf J}^{-1}_{\rm i} Δf_{\rm i} }
        Δx^{T}_{\rm i} {\bf J}^{-1}_{\rm i},
```

given in the following function:

```@docs
squirrel.JiSMU
```

The Sherman-Morrison update formula is called within the Broyden root
finding function:

```@docs
squirrel.bsolve
```

In the `squirrel.jl` code, the argument `F` in the root finding function
`squirrel.bsolve` corresponds to the function `squirrel.zF` (the
function itself is called as an argument), the argument `J` corresponds
to the initial value of the Jacobian, `f0` and `x0` correspond to the
respective initial guesses for ``f(v)`` and ``v``.

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
can generate ``C(n_{\rm e},4)`` sets of emission points (where ``C(n,k)`` is
the binomial coefficient):

```@docs
squirrel.combX
```

One feeds the output of `combX` into the `locator4` function to obtain
``C(n_{\rm e},4)`` solutions. The following outlier detection formula, which
is based on the comparison with median values, can then be used to
exclude large errors (the assumption here is that the large errors are
infreqent):

```@docs
squirrel.odetc
```

The following function implements the procedure described above, given a matrix of ``4×n_{\rm e}`` emission points ``X``:

```@docs
squirrel.locator
```
