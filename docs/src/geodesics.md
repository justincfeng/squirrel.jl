# Geodesic solver

## Null Geodesics: Hamilton's equations

For a spacetime geometry described by a metric ``g_{μν}`` and its
inverse ``g^{μν}``, null geodesics may be described by the Hamiltonian:

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

The conjugate momenta ``p_μ`` are given by:

```math
p_μ = g_{μν} \frac{dx^ν}{dλ},
```

and for null geodesics, the initial data satisfies:

```math
\left.δ{^i}{_μ} \frac{dx^μ}{dλ}\right|_{λ=0} = {\bf v}^i
```
```math
\left. g_{μν}({x}_0) 
\frac{dx^μ}{dλ}\frac{dx^ν}{dλ} \right|_{λ=0} = 0,
```

where ``{\bf v}`` denotes the spatial components of the initial 
four-velocity vector.

The solution to Hamilton's equations is formally given by
``x^μ=x^μ(λ,{x}_0,{\bf v})`` (with ``x_0`` denoting the initial position
of the geodesic), and since ``λ`` (being an affine parameter) can be
redefined linearly, it is appropriate to set up the problem so that
``λ∈[0,1]``, with ``λ=0`` being the initial point and ``λ=1`` is the
final point.

## Implementation

### Hamiltonian

The Hamiltonian function takes the form:
```@docs
squirrel.HamGeo
```

### Hamilton's equations
Hamilton's equations may be written in terms of the phase space
coordinate ``z^α``, where ``z=(x,p)``,

```math
\frac{dz^α}{dλ} = J^{αβ} \frac{∂H}{∂z^β}.
```

Here ``J^{αβ}`` is the symplectic matrix, which has the block matrix
form:

```math
J^{\cdot\cdot} =
\left[
  \begin{array}{cccc}
     O  &  I  \\
     -I  &  O  
  \end{array}
\right],
```

where ``O`` is a ``4×4`` matrix of zeros and ``I`` is the identity
matrix. The symplectic matrix is implemented as an operator acting on a
vector ``{∂H}/{∂z^β}``:

```@docs
squirrel.Jsympl
```

The quantity ``J^{αβ} \frac{∂H}{∂z^β}`` is evaluated in the following
function: 

```@docs
squirrel.ZdotGeo
```

## Geodesic solver function

Geodesics are solved with the following function, which outputs the
endpoint (``λ=1``) of the solution to Hamilton's equations:

```@docs
squirrel.solveZ
```

Following the recommendations in the [ODE
Solver](https://diffeq.sciml.ai/stable/solvers/ode_solve/#ode_solve)
documentation for the
[`OrdinaryDiffEq.jl`](https://github.com/SciML/OrdinaryDiffEq.jl)
library, the integrators `AutoVern7(Rodas5())` and `AutoVern9(Rodas5())`
are used in `squirrel.jl`.
