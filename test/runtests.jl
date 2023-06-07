using HdivBiotElasticityPaper
using Test

@testset "ConvergenceTests.jl" begin
    include("ConvergenceTests.jl")
end

@testset "TransientHdivBiotElasticityTests.jl" begin
    include("TransientHdivBiotElasticityTests.jl")
end

