# Geodesic solver

## Null Geodesics

Null geodesics in a spacetime geometry described by a metric ``g_{μν}``
are described by the Hamiltonian:

```math
H = \frac{1}{2} g^{μν} p_μ p_ν
```

and the associated Hamilton equations:

```math
\frac{dx^μ}{dλ} = \frac{∂H}{∂p_μ}
```
```math
\frac{dp_μ}{dλ} = - \frac{∂H}{∂x^μ}.
```

The conjugate momenta ``p_μ`` are defined as:

```math
p_μ := g_{μν} \frac{dx^ν}{dλ},
```

and for null geodesics, the initial data satisfies:

```math
\left.δ{^i}{_μ} \frac{dx^μ}{dλ}\right|_{λ=0} = {\bf v}^i
```
```math
\left. g_{μν}({x}_0) 
\frac{dx^μ}{dλ}\frac{dx^ν}{dλ} \right|_{λ=0} = 0.
```

The solution to Hamilton's equations is formally given by
``x^μ=x^μ(λ,{x}_0,{\bf v})``, and since ``λ`` (being an affine
parameter) can be redefined linearly, it is appropriate to set up the
problem so that ``λ∈[0,1]``, with ``λ=0`` being the initial point and
``λ=1`` is the final point.

## Implementation

### Hamiltonian

The Hamiltonian function takes the form:
```@docs
squirrel.HamGeo
```

### Hamilton's equations
Hamilton's equations may be written in terms of the phase space coordinate ``z^α``, where ``z=(x,p)``.

```math
\frac{dz^α}{dλ} = J^{αβ} \frac{∂H}{∂z^β}.
```

where ``J^{αβ}`` is the symplectic matrix, which has the form:

```math
J =
\left[
  \begin{array}{cccc}
     O  &  I  \\
     -I  &  O  
  \end{array}
\right].
```


## Initial data

```@docs
squirrel.nullenforcerp
```

```@docs
squirrel.nullenforcerp
```

## Geodesic functions

```@docs
squirrel.HamGeo
```

