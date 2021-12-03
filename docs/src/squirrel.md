# Squirrel algorithm

## Root finding algorithm

Given some function ``f(x)``, the standard approach is the Newton method, which finds the roots according to the prescription:

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

The Broyden update formula is implemented in the following function.

```@docs
squirrel.bsolve
```

## Geodesic integration


## 