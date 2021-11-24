using Test

using LinearAlgebra, Combinatorics, Statistics
using OrdinaryDiffEq, ForwardDiff, DiffResults

include("../src/type.jl")
include("testfunctions.jl")

include("../src/cereal.jl")
include("../src/metrics/Minkowski.jl")
include("../src/geosol.jl")
include("../src/broyden.jl")
include("../src/outlier.jl")
include("../src/squirrel.jl")
include("../src/evalsq.jl")

@testset "All tests:" begin

    @time @testset "Minkowski tests:" begin 
            include("minkowski_test.jl") 
        end

    @time @testset "cereal tests:" begin 
            include("cereal_test.jl") 
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

    @time @testset "eval tests:" begin 
            include("evalsq_test.jl") 
        end
        
end

nothing
