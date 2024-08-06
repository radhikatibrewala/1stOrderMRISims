"""
get_interpolation_objects(input_grid, B₀, BGx, BGy, BGz, B1₋, B1₊)

Prepares interpolation objects for all the fields required for the simulation 

# Arguments
- `input_grid::Vector{StepRangeLen{Float64, {Float64}, {Float64}, Int64}}`: A vector containing three 1D arrays representing the grid positions in the x, y, and z directions.
- `B₀::Vector{Array{Float64, 3}}`: A vector of the provided main magnetic field x,y,z components
- `BGx::Vector{Array{Float64, 3}}`: A vector of the provided x-gradient field x,y,z components
- `BGy::Vector{Array{Float64, 3}}`: A vector of the provided y-gradient field x,y,z components
- `BGz::Vector{Array{Float64, 3}}`: A vector of the provided z-gradient field x,y,z components
- `B1₋::Array{Float64, 3}`: all ones for the B1 transmit field B1₊
- `B1₊::Vector{Array{ComplexF64, 2}}`: all ones for two recieve coils for the B1 recieve field B1₋.

# Returns
- `B₀_interpolation_object::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the main magnetic field x,y,z components
- `BGx_interpolation_object::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the x-gradient field x,y,z components
- `BGy_interpolation_object::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the y-gradient field x,y,z components
- `BGz_interpolation_object::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the z-gradient field x,y,z components
- `B1₋_interpolation_object::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the values of the B1 receive field B1₋ for each coil
- `B1₊_interpolation_object::Extrapolation{Float64, 3, ...}`:  A vector of the interpolation object for the values of the B1 transmit field B1₊
"""
function get_interpolation_objects(input_grid, B₀, BGx, BGy, BGz, B1₋, B1₊)

    B₀x = Interpolations.extrapolate(scale(interpolate(B₀[1], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())
    B₀y = Interpolations.extrapolate(scale(interpolate(B₀[2], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())
    B₀z = Interpolations.extrapolate(scale(interpolate(B₀[3], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())

    BGx_x = Interpolations.extrapolate(scale(interpolate(BGx[1], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())
    BGx_y = Interpolations.extrapolate(scale(interpolate(BGx[2], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())
    BGx_z = Interpolations.extrapolate(scale(interpolate(BGx[3], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())

    BGy_x = Interpolations.extrapolate(scale(interpolate(BGy[1], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())
    BGy_y = Interpolations.extrapolate(scale(interpolate(BGy[2], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())
    BGy_z = Interpolations.extrapolate(scale(interpolate(BGy[3], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())

    BGz_x = Interpolations.extrapolate(scale(interpolate(BGz[1], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())
    BGz_y = Interpolations.extrapolate(scale(interpolate(BGz[2], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())
    BGz_z = Interpolations.extrapolate(scale(interpolate(BGz[3], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())

    B1₋_r1 = Interpolations.extrapolate(scale(interpolate(B1₋[1], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())
    B1₋_r2 = Interpolations.extrapolate(scale(interpolate(B1₋[2], BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())

    B1₊ = Interpolations.extrapolate(scale(interpolate(B1₊, BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())
    B₀_interpolation_object = [B₀x, B₀y, B₀z]
    BGx_interpolation_object = [BGx_x, BGx_y, BGx_z]
    BGy_interpolation_object = [BGy_x, BGy_y, BGy_z]
    BGz_interpolation_object = [BGz_x, BGz_y, BGz_z]
    B1₋_interpolation_object =  [B1₋_r1, B1₋_r2]
    B1₊_interpolation_object =   B1₊
    return B₀_interpolation_object, BGx_interpolation_object, BGy_interpolation_object, BGz_interpolation_object, B1₋_interpolation_object, B1₊_interpolation_object
end
"""
interpolate_fields(input_grid, B₀, BGx, BGy, BGz, B1₋, B1₊)

Interpolates the fields using the interpolation objects and stores all the values 
# Arguments
- `simulation_grid::Vector{StepRangeLen{Float64, {Float64}, {Float64}, Int64}}`: A vector containing three 1D arrays representing the grid positions in the x, y, and z directions.
- `B₀_I::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the main magnetic field x,y,z components
- `BGx_I::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the x-gradient field x,y,z components
- `BGy_I::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the y-gradient field x,y,z components
- `BGz_I::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the z-gradient field x,y,z components
- `B1₋_I::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the values of the B1 receive field B1₋ for each coil
- `B1₊_I::Extrapolation{Float64, 3, ...}`:  A vector of the interpolation object for the values of the B1 transmit field B1₊

# Returns
- `B₀_interpolated::Vector{Array{Float64, 3}}`: A vector of the interpolated main magnetic field x,y,z components
- `BGx_interpolated::Vector{Array{Float64, 3}}`: A vector of the interpolated x-gradient field x,y,z components
- `BGy_interpolated::Vector{Array{Float64, 3}}`: A vector of the interpolated y-gradient field x,y,z components
- `BGz_interpolated::Vector{Array{Float64, 3}}`: A vector of the interpolated z-gradient field x,y,z components
- `B1₋_interpolated::Array{Float64, 3}`: A vector of the interpolated values of the B1 receive field B1₋ for each coil
- `B1₊_interpolated::Vector{Array{Float64, 3}}`:  A vector of the interpolated values of the B1 transmit field B1₊
"""
function interpolate_fields(B₀_I, BGx_I, BGy_I, BGz_I, B1₋_I, B1₊_I, simulation_grid)
    
    B₀x, B₀y, B₀z = B₀_I[1], B₀_I[2], B₀_I[3]
    BGx_x, BGx_y, BGx_z = BGx_I[1], BGx_I[2], BGx_I[3]
    BGy_x, BGy_y, BGy_z = BGy_I[1], BGy_I[2], BGy_I[3]
    BGz_x, BGz_y, BGz_z = BGz_I[1], BGz_I[2], BGz_I[3]
    
    grid_idcs = CartesianIndices(((length(simulation_grid[1])), length((simulation_grid[2])), length((simulation_grid[3]))))

    B₀x_ = [B₀x(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]
    B₀y_ = [B₀y(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]
    B₀z_ = [B₀z(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]

    BGx_x_ = [BGx_x(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]
    BGx_y_ = [BGx_y(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]
    BGx_z_ = [BGx_z(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]

    BGy_x_ = [BGy_x(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]
    BGy_y_ = [BGy_y(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]
    BGy_z_ = [BGy_z(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]

    BGz_x_ = [BGz_x(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]
    BGz_y_ = [BGz_y(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]
    BGz_z_ = [BGz_z(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]
    
    B1₋_r1 = [B1₋_I[1](simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]
    B1₋_r2 = [B1₋_I[2](simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]

    B1₊_ = [B1₊_I(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]
    
    B₀_interpolated =  [B₀x_, B₀y_, B₀z_]
    BGx_interpolated = [BGx_x_, BGx_y_, BGx_z_]
    BGy_interpolated = [BGy_x_, BGy_y_, BGy_z_]
    BGz_interpolated = [BGz_x_, BGz_y_, BGz_z_]
    B1₋_interpolated = [B1₋_r1, B1₋_r2]
    B1₊_interpolated =   B1₊_
    return B₀_interpolated, BGx_interpolated, BGy_interpolated, BGz_interpolated, B1₋_interpolated, B1₊_interpolated
end
"""
get_interpolated_m₀(input_params, simulation_grid)

Interpolates the fields using the interpolation objects and stores all the values 
# Arguments
- `input_params::Dict{Any, Any}`: A dictionary containing various input parameters including paths to configuration files, nominal field values, and limits for the field initialization.
- `simulation_grid::Vector{StepRangeLen{Float64, {Float64}, {Float64}, Int64}}`: A vector containing three 1D arrays representing the grid positions in the x, y, and z directions.

# Returns
- `m₀::Array{Float64, 3}`: Input image interpolated to create the initial magnetization: m0

"""
function get_interpolated_m₀(input_params, simulation_grid)

    im_3d = load_image_file(input_params["input_files"]["image_file"]);
    img_I = get_m₀(im_3d);

    grid_idcs = CartesianIndices(((length(simulation_grid[1])), length((simulation_grid[2])), length((simulation_grid[3]))))
    m₀  = [img_I(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]

    return m₀
end

"""
interpolate_field_gradients(B₀_I, BGx_I, BGy_I, BGz_I, simulation_grid)

Interpolates the fields using the interpolation objects and stores all the values 
# Arguments
- `simulation_grid::Vector{StepRangeLen{Float64, {Float64}, {Float64}, Int64}}`: A vector containing three 1D arrays representing the grid positions in the x, y, and z directions.
- `B₀_I::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the main magnetic field x,y,z components
- `BGx_I::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the x-gradient field x,y,z components
- `BGy_I::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the y-gradient field x,y,z components
- `BGz_I::Vector{Extrapolation{Float64, 3, ...}`: A vector of the interpolation object for the z-gradient field x,y,z components

# Returns
- `∂B₀_interpolated::Vector{Array{Float64, 3}}`: A vector of the interpolated gradients of the main magnetic field for all x,y,z components w.r.t x-y-z
- `∂BGx_interpolated::Vector{Array{Float64, 3}}`:A vector of the interpolated gradients of the x-gradient field for all x,y,z components w.r.t x-y-z
- `∂BGy_interpolated::Vector{Array{Float64, 3}}`:A vector of the interpolated gradients of the y-gradient field for all x,y,z components w.r.t x-y-z
- `∂BGz_interpolated::Vector{Array{Float64, 3}}`:A vector of the interpolated gradients of the z-gradient field for all x,y,z components w.r.t x-y-z

"""
function interpolate_field_gradients(B₀_I, BGx_I, BGy_I, BGz_I, simulation_grid)
    B₀x_I, B₀y_I, B₀z_I = B₀_I[1], B₀_I[2], B₀_I[3]
    BGx_x_I, BGx_y_I, BGx_z_I = BGx_I[1], BGx_I[2], BGx_I[3]
    BGy_x_I, BGy_y_I, BGy_z_I = BGy_I[1], BGy_I[2], BGy_I[3]
    BGz_x_I, BGz_y_I, BGz_z_I = BGz_I[1], BGz_I[2], BGz_I[3]

    grid_idcs = CartesianIndices(((length(simulation_grid[1])), length((simulation_grid[2])), length((simulation_grid[3]))))

    ∂B₀x = [[Interpolations.gradient(B₀x_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[1] for img_idx ∈ grid_idcs],
    [Interpolations.gradient(B₀x_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[2] for img_idx ∈ grid_idcs],
    [Interpolations.gradient(B₀x_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[3] for img_idx ∈ grid_idcs]]

    ∂B₀y = [[Interpolations.gradient(B₀y_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[1] for img_idx ∈ grid_idcs],
    [Interpolations.gradient(B₀y_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[2] for img_idx ∈ grid_idcs],
    [Interpolations.gradient(B₀y_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[3] for img_idx ∈ grid_idcs]]

    ∂B₀z = [[Interpolations.gradient(B₀z_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[1] for img_idx ∈ grid_idcs],
    [Interpolations.gradient(B₀z_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[2] for img_idx ∈ grid_idcs],
    [Interpolations.gradient(B₀z_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[3] for img_idx ∈ grid_idcs]]

    ∂BGx_x = [[Interpolations.gradient(BGx_x_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[1]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGx_x_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[2]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGx_x_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[3]  for img_idx ∈ grid_idcs]]
    
    ∂BGx_y = [[Interpolations.gradient(BGx_y_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[1]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGx_y_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[2]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGx_y_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[3]  for img_idx ∈ grid_idcs]]

    ∂BGx_z = [[Interpolations.gradient(BGx_z_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[1]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGx_z_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[2]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGx_z_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[3]  for img_idx ∈ grid_idcs]]
    
    ∂BGy_x = [[Interpolations.gradient(BGy_x_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[1]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGy_x_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[2]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGy_x_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[3]  for img_idx ∈ grid_idcs]]
    
    ∂BGy_y = [[Interpolations.gradient(BGy_y_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[1]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGy_y_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[2]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGy_y_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[3]  for img_idx ∈ grid_idcs]]

    ∂BGy_z = [[Interpolations.gradient(BGy_z_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[1]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGy_z_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[2]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGy_z_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[3]  for img_idx ∈ grid_idcs]]

    ∂BGz_x = [[Interpolations.gradient(BGz_x_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[1]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGz_x_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[2]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGz_x_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[3]  for img_idx ∈ grid_idcs]]
    
    ∂BGz_y = [[Interpolations.gradient(BGz_y_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[1]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGz_y_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[2]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGz_y_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[3]  for img_idx ∈ grid_idcs]]

    ∂BGz_z = [[Interpolations.gradient(BGz_z_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[1]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGz_z_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[2]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(BGz_z_I, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[3]  for img_idx ∈ grid_idcs]]

    ∂B₀_interpolated = [∂B₀x, ∂B₀y, ∂B₀z]
    ∂BGx_interpolated = [∂BGx_x, ∂BGx_y, ∂BGx_z]
    ∂BGy_interpolated = [∂BGy_x, ∂BGy_y, ∂BGy_z]
    ∂BGz_interpolated = [∂BGz_x, ∂BGz_y, ∂BGz_z]

    return ∂B₀_interpolated, ∂BGx_interpolated, ∂BGy_interpolated, ∂BGz_interpolated
end
"""
calculate_ω_∂ω(input_params, B₀_, BGz_, input_grid, ΔB₀, simulation_grid)

Interpolates the spatially varying Larmor Frequency at the time of the RF pulse at all the simulation locations  
# Arguments

- `input_grid::Vector{StepRangeLen{Float64, {Float64}, {Float64}, Int64}}`: A vector containing three 1D arrays representing the grid positions in the x, y, and z directions.
-` simulation_grid::Vector{StepRangeLen{Float64, {Float64}, {Float64}, Int64}}`: A vector containing three 1D arrays representing the grid positions in the x, y, and z directions.
- `B₀_::Vector{Array{Float64, 3}}`: A vector of the provided main magnetic field x,y,z components
- `ΔB₀::Float64`: The center frequency value ΔB₀. For 2D imaging, it is the field strength at the specified slice_location in the z-direction. For 3D imaging, it is the field strength deviation at the origin (0, 0, 0).
- `BGz_::Vector{Array{Float64, 3}}`: A vector of the provided z-gradient field x,y,z components
- `input_params::Dict{Any, Any}`: A dictionary containing various input parameters including paths to configuration files, nominal field values, and limits for the field initialization.


# Returns
- `ω_interpolated::Vector{Array{Float64, 3}}`: A vector of the interpolated of the spatially varying Larmor Frequency at the time of the RF pulse
- `∂ω_interpolated::Vector{Array{Float64, 3}}`:A vector of the interpolated gradients of ω for all x,y,z components w.r.t x-y-z


"""
function calculate_ω_∂ω(input_params, B₀_, BGz_, input_grid, ΔB₀, simulation_grid)
    grid_idcs = CartesianIndices(((length(input_grid[1])), length((input_grid[2])), length((input_grid[3]))))
    scale_z = input_params["BGz_"]
    ω_ = (((B₀_[1] + scale_z*BGz_[1]).^2 + (B₀_[2] + scale_z*BGz_[2]).^2 + (B₀_[3] + scale_z*BGz_[3]).^2).^0.5) .- ΔB₀    
    ω = Interpolations.extrapolate(scale(interpolate(ω_, BSpline(Linear())), input_grid[1], input_grid[2],input_grid[3]),Interpolations.Line())
    
    grid_idcs = CartesianIndices(((length(simulation_grid[1])), length((simulation_grid[2])), length((simulation_grid[3]))))
    ω_interpolated  = [ω(simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]],simulation_grid[3][img_idx[3]]) for img_idx ∈ grid_idcs]
    ∂ω_interpolated = [[Interpolations.gradient(ω, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[1]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(ω, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[2]  for img_idx ∈ grid_idcs],
    [Interpolations.gradient(ω, simulation_grid[1][img_idx[1]], simulation_grid[2][img_idx[2]], simulation_grid[3][img_idx[3]])[3]  for img_idx ∈ grid_idcs]]

    return ω_interpolated, ∂ω_interpolated
end