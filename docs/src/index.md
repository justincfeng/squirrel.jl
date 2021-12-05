# Home

## Introduction

Relativistic positioning refers to the concept of establishing spacetime
positions from proper time broadcasts emitted by a system of satellites.
Central to relativistic positioning is the relativistic location
location problem, which is the problem of finding the intersection of
future pointing light cones from a collection of at least four emission
points. `squirrel.jl` contains a collection of functions for the
relativistic location problem in slightly curved spacetime geometries.

## Short tutorial

### Setup

The `squirrel.jl` code was written for and tested in Julia 1.6; we
recommend Julia 1.6 or newer. To add the code, run the following command
in the pagkage manager for the Julia `REPL`:

    pkg> add https://github.com/justincfeng/squirrel.jl/

    julia> using squirrel

### Positioning examples

#### Vacuum case

One begins by defining a metric. The Kerr-Schild metric for a rotating object with the mass and angular momentum of the Earth is given by the following function:

    julia> gks = squirrel.metric.ge

The resulting function `gks` takes a four-component vector
`xc=[t,x,y,z]` (representing spacetime coordinates) as an argument, and
returns a ``4×4`` matrix of metric components. The units are normalized
to Earth mass (1 length unit = 0.4435 cm), and chosen so that the speed
of light is 1. One can randomly generate a target point `Xtar` and a set
of ``5`` emission points `X` on the past light cone of `Xtar` with the
following:

    julia> (X,Xtar) = squirrel.seval.pgen(6e9,gks,1e-14,5) ;

The first argument in the function on the left hand side specifies
(roughly) the spatial radius of the emission points from the origin, and
the third argument specified the tolerance for the geodesic integrator.
The target point `Xtar` is placed on the WGS84 reference ellipsoid
(defined with respect to the Cartesian Kerr-Schild coordinates).

The intersection of future light cones from the emission points is computed with the `squirrel.locator` function:

    julia> Xs = squirrel.locator(X,gks,1e-10)

The third argument is the tolerance for the geodesic solvers; the tolerance is looser here to minimize computation time. The accuracy of the result may be estimated by comparing `Xs` and `Xtar`:

    julia> ΔX = Xs-Xtar

Upon multiplying by the conversion factor 0.4435 to obtain the result in units of centimeters (`0.4435*ΔX`), one typically obtains a result well in the submillimeter range. The relative error may be obtained by running:

    julia> squirrel.norm(ΔX)/squirrel.norm(Xtar)

and one typically obtains errors on the order of ``\sim 10^{-12} - 10^{-13}``.

#### Including atmospheric and ionospheric effects

One can incorporate the effects of the atmosphere and ionosphere with the
following metric:

    julia> g = squirrel.metric.g

This metric is the Gordon metric for light propagatation through media; here, a simple model for the atmosphere and ionosphere is implemented. One may repeat the steps of the vacuum case to obtain the errors:

    julia> (X,Xtar) = squirrel.seval.pgen(6e9,g,1e-14,5) ;

    julia> Xs = squirrel.locator(X,g,1e-10)

    julia> ΔX = Xs-Xtar

    julia> 0.4435*ΔX

    julia> squirrel.norm(ΔX)/squirrel.norm(Xtar)

In this case, one typically finds that the errors are larger than the vacuum case by an order of magnitude, though still in the submillimeter range.

### Evaluation

Scripts are provided for a more comprehensive evaluation of the accuracy
and performance of the functions in `squirrel.jl`.

## References

