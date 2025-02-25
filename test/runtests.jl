using Test

using LinearAlgebra, Combinatorics, Statistics
using OrdinaryDiffEq, ForwardDiff, DiffResults
using DoubleFloats

include("../src/type.jl")
include("testfunctions.jl")

include("../src/metrics/Minkowski.jl")
include("../src/metrics/KerrSchild.jl")
include("../src/metric.jl")
include("../src/geosol.jl")
include("../src/broyden.jl")
include("../src/outlier.jl")
include("../src/squirrel.jl")

η = (X,p=Float64[]) -> ημν(typeof(X[1]),p)

@testset "All tests:" begin

    @time @testset "Minkowski tests:" begin 
            include("minkowski_test.jl")
        end

    @time @testset "geosol tests:" begin 
            include("geosol_test.jl")
        end

    @time @testset "broyden tests:" begin 
            include("broyden_test.jl")
        end

    @time @testset "outlier detection tests:" begin 
            include("outlier_test.jl")
        end
    
    @time @testset "squirrel tests:" begin 
            include("squirrel_test.jl") 
        end

    @time @testset "seval tests:" begin 
            include("seval_test.jl") 
    end

end

nothing
