function SpinEchoTest(filepaths)

    config = YAML.load_file(filepaths["config_file"])

    fields = matread(filepaths["field_file"])    
    input_params = config["simulation_parameters_2D"]
    input_params["nominal_field"] = false
    input_params["imaging_type"] = "2D"
    input_params["sequence_name"]  = "spin_echo"

    config["input_files"]["image_file"] = filepaths["image_file"]
    input_grid, simulation_grid = get_required_grids(input_params, fields)   
    B₀, BGx, BGy, BGz, B1₊, B1₋ = initialise_fields(input_params, input_params["nominal_field"], fields, input_params["lower_lims"], input_params["upper_lims"]) 

    ΔB₀ =  calculate_ΔB₀(B₀, input_grid, input_params)      
    @test ΔB₀ ≈ 0.099 rtol = 1e-2
    B₀_I, BGx_I, BGy_I, BGz_I, B1₋_I, B1₊_I = get_interpolation_objects(input_grid, B₀, BGx, BGy, BGz, B1₋, B1₊)
    normalized_BGz =  Interpolations.gradient(BGz_I[3],0,0,0)[3]
    @test normalized_BGz ≈ 1.0 atol = 0

    B₀_, BGx_, BGy_, BGz_, B1₋_, B1₊_ = interpolate_fields(B₀_I, BGx_I, BGy_I, BGz_I, B1₋_I, B1₊_I, simulation_grid)
    ∂B₀, ∂BGx, ∂BGy, ∂BGz = interpolate_field_gradients(B₀_I, BGx_I, BGy_I, BGz_I, simulation_grid)
    normalized_∂BGz = ∂BGz[3][3][70,70,70]
    @test normalized_∂BGz ≈ 1.0 atol = 0

    center_∂B₀ = ∂B₀[3][3][70,70,70]
    @test center_∂B₀ ≈ 0 atol = 0

    m₀ = get_interpolated_m₀(config, simulation_grid)  
    @test m₀[70,70,70] ≈ 2.386e-05 rtol = 1e-03

    B1₊_ = B1₊_/B1₊_I(0,0,0) 
    ω_, ∂ω = calculate_ω_∂ω(input_params, B₀, BGz, input_grid, ΔB₀, simulation_grid)

    @test ω_[70,70,70] ≈ -1.2e-5 rtol = 1e-03
    @test ∂ω[3][70,70,70] ≈ input_params["BGz_"] rtol = 1e-06

    σ =  input_params["BGz_"] * input_params["Δz"] # tesla
    @test σ ≈ 2.0e-5 rtol = 1e-06

    sim_vars_ = store_sim_vars(B₀_, BGx_, BGy_, BGz_, ∂B₀, ∂BGx, ∂BGy, ∂BGz, ω_, ∂ω, m₀, B1₋_, B1₊_,ΔB₀, σ, input_params["order"])

    sequence_params, echo_time = get(sequence_functions, input_params["sequence_name"], ()->error("Invalid sequence name."))(input_params)
    @test sequence_params[1].rf_ == 90.0 
    @test sequence_params[end].indexKx == 140
    @test sequence_params[end].indexKy == 140
    @test sequence_params[16000].gradx ≈ -0.002609 rtol = 1e-03

    new_TR_line_indexes, TR_lines = indexes_per_TR(sequence_params)
    split_TR_line_indexes,  split_TR_lines = makechunks(new_TR_line_indexes, TR_lines, 1)
    current_TR_lines_indexes, current_TR_lines = split_TR_line_indexes[1], split_TR_lines[1] 

    num_timepoints = length(sequence_params)


    mr_signal_measurements = [ zeros(Complex{Float64}, length(sim_vars_.B1₋), num_timepoints)][1]
    @test size(mr_signal_measurements, 2) == 49280
    count = 0
    garbage_handling(count, input_params["grid_cube_size"])

    tr = 1
    flush(stdout)

    ϕ = zeros(Float64,size(sim_vars_.B₀[1]))
    ∂ϕ = [zeros(Float64,size(sim_vars_.B₀[1])),zeros(Float64,size(sim_vars_.B₀[1])),zeros(Float64,size(sim_vars_.B₀[1]))]

    start_idx = 1
    end_idx = 10

    for index in range(start_idx, 10)
        mr_signal_measurements[:,index], ϕ, ∂ϕ =  calculate_kpsace_value(index, sequence_params, sim_vars_, ϕ, ∂ϕ, input_params)
        count += 1
        garbage_handling(count, input_params["grid_cube_size"])   
    end

    @test abs(mr_signal_measurements[1,2]) ≈ 0.0020 rtol = 1e-02
    @test abs(mr_signal_measurements[1,8]) ≈ 8.1454e-5 rtol = 1e-04
    @test abs(mr_signal_measurements[1,10]) ≈ 1.1884e-5 rtol = 1e-4
end    