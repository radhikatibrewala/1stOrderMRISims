using Pkg
Pkg.activate(".")
using ApproxFun
using Interpolations
using MAT
using Test
using Plots
using FFTW
using JLD2
using ArgParse
using SpecialFunctions
using BenchmarkTools
using Test
using Dates
using StaticArrays
using Printf 
using YAML

include("src/FieldInterpolations.jl")
include("src/FirstOrderSolutions.jl")
include("src/InitializeSequences.jl")
include("src/InputFunctions.jl")
include("src/KspaceFunctions.jl")
include("src/RunSequence.jl")
include("src/SequenceLogic.jl")
include("src/SimulationVariables.jl")
include("src/Structs.jl")
include("src/Plotting.jl")

sequence_functions = Dict(
    "gradient_echo" => gradient_echo,
    "spin_echo" => spin_echo,
    "single_shot_spin_echo_epi"=> single_shot_spin_echo_epi,
    "turbo_spin_echo" => turbo_spin_echo,
    "single_shot_gradient_echo_epi"=> single_shot_gradient_echo_epi,
    "gradient_echo_3d" => gradient_echo_3d
)

function input_variables(config)
    if config["imaging_type"] == "2D"
        params = config["simulation_parameters_2D"]
        params["imaging_type"] = "2D"
    elseif config["imaging_type"] == "3D"
        params = config["simulation_parameters_3D"]
        params["imaging_type"] = "3D"
    else
        print("Incorrect imaging type, should be 2D or 3D [string]")
    end
    return params
end

function parse_config_sim(config)
    input_params = input_variables(config)
    input_params["nominal_field"] = false

    if config["run_sim_nominal"]
        input_params["nominal_field"] = true
    end
    input_params["input_files"] = config["input_files"]
    for keyname in keys(input_params["input_files"])
        input_params["input_files"][keyname] = abspath(input_params["input_files"][keyname])
    end

    return input_params
end

gather_kspace_functions = Dict(
    "2D" => allocate_kspace,
    "3D" => allocate_kspace_3D)

plot_sequence_functions = Dict(
        "2D" => plot_sequence,
        "3D" => plot_sequence_3d)

plot_result_functions = Dict(
    "2D" => plot_image,
    "3D" => plot_image)

rss_functions = Dict(
    "2D" => rss_coil_combine_2D,
    "3D" => rss_coil_combine_3D)