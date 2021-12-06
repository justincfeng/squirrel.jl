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

where the indices ``I∈\\{1,2,...,n_e\\}`` distinguish the emission
points ``X_I`` and their associated null geodesics, and ``{\bf v}_I``
denotes the spatial components of the initial four-velocity vector for
the null geodesic (the time component is determined by the requirement
that the four-velocity is null). 

Given a collection of four geodesic functions ``\\{x_1,x_2,x_3,x_4\\}``,
the condition that they intersect is the vanishing of the following
vector-valued function:

```math
f := \left( x_1 - x_2 , x_1 - x_3 , x_1 - x_4 \right)
```

which may be formally written as ``f=f(X_1,X_2,X_3,X_4,{\bf v}_1,{\bf
v}_2,{\bf v}_3,{\bf v}_4)``. The squirrel algorithm seeks to find the
roots of the function ``f``.

The squirrel algorithm is then summarized:

1.  Apply the flat spacetime algorithm to the emission points `X` to
    obtain an initial guess for the zeros of the function `f`.

2.  Apply a root finding algorithm to the function `f`.

## Geodesic endpoints and Jacobian

To find the roots of the function `f`, the squirrel algorithm first
computes the Jacobian of `f`.

```@docs
squirrel.gsolve
```

```@docs
squirrel.zF
```

```@docs
squirrel.gejac
```

```@docs
squirrel.geocJ
```

```@docs
squirrel.idf
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



## Locator

```@docs
squirrel.locator4
```

```@docs
squirrel.locator
```