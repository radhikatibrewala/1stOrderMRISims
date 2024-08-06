include("Initialize.jl")
using Test
include("tests/GradientEchoTest.jl")
include("tests/SpinEchoTest.jl")

#! Do not change filepaths

filepaths = Dict()
filepaths["config_file"] = "tests/GradientEchoTest.yaml"
filepaths["field_file"] = "examples/quadratic_field.mat"
filepaths["image_file"]= "examples/fastmri_052.mat"

@testset "Gradient Echo" begin
    GradientEchoTest(filepaths)
end

filepaths["config_file"] = "tests/SpinEchoTest.yaml"
filepaths["field_file"] = "examples/quadratic_field.mat"
filepaths["image_file"]= "examples/fastmri_052.mat"

@testset "Spin Echo" begin
    SpinEchoTest(filepaths)
end



