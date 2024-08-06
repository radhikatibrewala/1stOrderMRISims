using Parameters
"""
    read_input_files(input_params)

Reads field configuration data from the specified input file and initializes the field values.

# Arguments
- `input_params::Dict{Any, Any}`: A dictionary containing various input parameters including paths to configuration files, nominal field values, and limits for the field initialization.

# Returns
- `fields::Dict{String, Any}`: The field configuration data loaded from the input file.
- `B₀::Vector{Array{Float64, 3}}`: The initialized main magnetic field.
- `BGx::Vector{Array{Float64, 3}}`: The initialized x-gradient magnetic field.
- `BGy::Vector{Array{Float64, 3}}`: The initialized y-gradient magnetic field.
- `BGz::Vector{Array{Float64, 3}}`: The initialized z-gradient magnetic field.
- `B1₊::Array{Float64, 3}`: The initialized B1+ field component.
- `B1₋::Vector{Array{ComplexF64, number_receive_coils}}`: The initialized B1- field component. 

# Notes:
1. See ../files/field_values.mat for example configuration format and headers
2. lower_lims should be a tuple containing the lower lims of the x,y,z dims of the scanner in meters; example = (-0.21,-0.21,-0.21)
3. upper_lims should be a tuple containing the upper lims of the x,y,z dims of the scanner in meters; example = (0.21,0.21,0.21)

"""

function read_input_files(input_params)
    fields = matread(input_params["input_files"]["field_config_file"])                                      
    B₀, BGx, BGy, BGz, B1₊, B1₋ = initialise_fields(input_params, input_params["nominal_field"], fields, input_params["lower_lims"], input_params["upper_lims"]) 
    
    return fields, B₀, BGx, BGy, BGz, B1₊, B1₋
end

"""
    get_index_lims(fields, lower_lims, upper_lims)

Finds the indexes where the field of interest is - i.e to reduce computation time, we dont want to run the simulation where there is no field/input

# Arguments
- `fields::Dict{String, Any}`: The field configuration data loaded from the input file.
- `lower_lims::Tuple{Float64, Float64, Float64}`: Lower lims of the x,y,z dims of the scanner in meters; example = (-0.21,-0.21,-0.21)
- `upper_lims::Tuple{Float64, Float64, Float64}`: Upper lims of the x,y,z dims of the scanner in meters; example = (-0.21,-0.21,-0.21)

# Returns
- `indexes_lims::Vector{Int64}`: Indexes of limits of desired field coverage
- `coords::Array{Float64, 4}`: Meshgrid of all inputted coordinates
- `num_points::Vector{Int64}`: Number of unique x,y,z locations
"""
function get_index_lims(fields, lower_lims, upper_lims)
    
    indexes_lims = intersect(findall(x -> x >= lower_lims[1], fields["x"]), findall(x -> x <= upper_lims[1], fields["x"]),
    findall(x -> x >= lower_lims[2], fields["y"]), findall(x -> x <= upper_lims[2], fields["y"]),
    findall(x -> x >= lower_lims[3], fields["z"]), findall(x -> x <= upper_lims[3], fields["z"]))
    
    x, y, z =  [fields["x"][indexes_lims], fields["y"][indexes_lims], fields["z"][indexes_lims]]

    num_points = [length(unique(x)), length(unique(y)),length(unique(z))]
    coords = cat(reshape_field(x, num_points), reshape_field(y, num_points), reshape_field(z, num_points), dims = 4)

    return indexes_lims, coords, num_points
end
"""
    get_required_grids(input_params, fields)

Generates the required grids based on the input parameters and field configurations.

# Arguments
- `input_params::Dict{Any, Any}`: A dictionary containing various input parameters including paths to configuration files, nominal field values, and limits for the field initialization.
- `fields::Dict{String, Any}:` The field configuration data loaded from the input file.

# Returns
- `G_pos::3-element Vector{StepRangeLen{Float64, {Float64}, {Float64}, Int64}}`: A vector containing three 1D arrays representing the grid positions in the x, y, and z directions.
- `grid_lims::Vector{StepRangeLen{Float64, {Float64}, {Float64}, Int64}}`: A vector containing three 1D arrays representing the grid limits in the x, y, and z directions, defined by the lower and upper limits with respect to the grid cube size.
"""
function get_required_grids(input_params, fields)
    
    indexes_lims, _, _ = get_index_lims(fields, input_params["lower_lims"], input_params["upper_lims"])
    x, y, z =  [fields["x"][indexes_lims], fields["y"][indexes_lims], fields["z"][indexes_lims]]
    num_points = [length(unique(x)), length(unique(y)),length(unique(z))]
    
    G_pos_x = range(minimum(x), maximum(x), num_points[1])
    G_pos_y = range(minimum(y), maximum(y), num_points[2])
    G_pos_z = range(minimum(z), maximum(z), num_points[3])

    x1d =  input_params["lower_lims"][1]:input_params["grid_cube_size"]:input_params["upper_lims"][1]                     
    y1d =  input_params["lower_lims"][2]:input_params["grid_cube_size"]:input_params["upper_lims"][2]                           
    z1d =  input_params["lower_lims"][3]:input_params["grid_cube_size"]:input_params["upper_lims"][3]                                 

    G_pos = [G_pos_x, G_pos_y, G_pos_z]
    grid_lims = [x1d, y1d, z1d]
    return G_pos, grid_lims
end
"""
    generate_simulation_variables(input_params)

Generates all the necessary variables for the simulation based on the provided input parameters.

# Arguments
- `input_params::Dict`: A dictionary containing all the required input parameters such as field configuration files, grid limits, imaging type, and sequence name.

# Returns
- `sim_vars_::Struct`: A dictionary containing all the generated simulation variables required for the MRI simulation, including field values, gradients, and interpolated values.
# note: see Structs.jl for details on the structure of sim_vars_

# Function Details
1. Reads the input field configurations and initializes the fields.
2. Generates the input and simulation grids based on input file and the simulation grid-cube size.
3. Computes the center frequency value(ΔB₀) and creates interpolation objects for the fields.
4. Interpolates the fields and their gradients to match the simulation grid.
5. Adjusts the B1+ field based on the sequence type, assuming B1+ equals 1 everywhere for spin echo sequences.
6. Calculates additional variables specific to 2D or 3D imaging types, Larmor freq at excitation (ω_) and its gradient (∂ω), or the spectral response (B1₊_Spectral) for 3D.
7. Stores and returns the simulation variables in a structured format depending on the imaging type.

"""
function generate_simulation_variables(input_params)
    fields, B₀, BGx, BGy, BGz, B1₊, B1₋ = read_input_files(input_params)
    input_grid, simulation_grid = get_required_grids(input_params, fields)          
    ΔB₀ =  calculate_ΔB₀(B₀, input_grid, input_params)      
    B₀_I, BGx_I, BGy_I, BGz_I, B1₋_I, B1₊_I = get_interpolation_objects(input_grid, B₀, BGx, BGy, BGz, B1₋, B1₊)
    B₀_, BGx_, BGy_, BGz_, B1₋_, B1₊_ = interpolate_fields(B₀_I, BGx_I, BGy_I, BGz_I, B1₋_I, B1₊_I, simulation_grid)
    ∂B₀, ∂BGx, ∂BGy, ∂BGz = interpolate_field_gradients(B₀_I, BGx_I, BGy_I, BGz_I, simulation_grid)
    m₀ = get_interpolated_m₀(input_params, simulation_grid)    #! work this out
    B1₊_ = B1₊_/B1₊_I(0,0,0) 
    
    """ small tip angle assumption does not hold for spin echo sequences, so assume B1+ = 1 everywhere"""
    if occursin("spin",input_params["sequence_name"])
        B1₊_ = ones(size(B1₊_))
    end

    if input_params["imaging_type"] == "2D"
        ω_, ∂ω = calculate_ω_∂ω(input_params, B₀, BGz, input_grid, ΔB₀, simulation_grid)
        σ =  input_params["BGz_"] * input_params["Δz"] # tesla
        sim_vars_ = store_sim_vars(B₀_, BGx_, BGy_, BGz_, ∂B₀, ∂BGx, ∂BGy, ∂BGz, ω_, ∂ω, m₀, B1₋_, B1₊_,ΔB₀, σ, input_params["order"])
    else
        B1₊_Spectral = get_spectral_response_3D_selection(B₀_, ΔB₀)       
        sim_vars_ = store_sim_vars_3d(B₀_, BGx_, BGy_, BGz_, ∂B₀, ∂BGx, ∂BGy, ∂BGz, m₀, B1₋_, B1₊_, B1₊_Spectral, ΔB₀, input_params["order"])
    end
    
    return sim_vars_
end

"""
    initialise_fields(input_params, nominal_field, fields, lower_lims = [-0.21, -0.21, -0.21], upper_lims = [0.21, 0.21, 0.21])

Initializes and prepares various magnetic field components and their interpolations based on provided parameters and field data. 

# Arguments
- `input_params::Dict{Any, Any}`: A dictionary containing various input parameters including paths to configuration files, nominal field values, and limits for the field initialization.
- `nominal_field::Bool`: A boolean flag indicating whether to use nominal field values (true) or actual field values (false).
-` fields::Dict{String, Any}`: The field configuration data loaded from the input file.
- `lower_lims::Tuple{Float64, Float64, Float64}`: Lower lims of the x,y,z dims of the scanner in meters; example = (-0.21,-0.21,-0.21)
- `upper_lims::Tuple{Float64, Float64, Float64}`: Upper lims of the x,y,z dims of the scanner in meters; example = (-0.21,-0.21,-0.21)


# Returns
- If nominal_field is true, returns a tuple containing:
  - `B₀_nominal::Vector{Array{Float64, 3}}`: A vector of nominal main magnetic field components with the z-component of B₀ set to 0.099 and xy-components set to zero.
  - `BGx_nominal::Vector{Array{Float64, 3}}`: A vector of nominal x-gradient field components with the z-component of calculated based on the grid coordinates and xy-components set to zero.
  - `BGy_nominal::Vector{Array{Float64, 3}}`: A vector of nominal y-gradient field components with the z-component of calculated based on the grid coordinates and xy-components set to zero.
  - `BGz_nominal::Vector{Array{Float64, 3}}`: A vector of nominal z-gradient field components with the z-component of calculated based on the grid coordinates and xy-components set to zero.
  - `B1₊_nominal::Array{Float64, 3}`: all ones for the B1 transmit field B1₊
  - `B1₋_nominal::Vector{Array{ComplexF64, 2}}`: all ones for two recieve coils for the B1 recieve field B1₋.

- If nominal_field is false, returns a tuple containing:
- `B₀::Vector{Array{Float64, 3}}`: A vector of the provided main magnetic field x,y,z components
- `BGx::Vector{Array{Float64, 3}}`: A vector of the provided x-gradient field x,y,z components
- `BGy::Vector{Array{Float64, 3}}`: A vector of the provided y-gradient field x,y,z components
- `BGz::Vector{Array{Float64, 3}}`: A vector of the provided z-gradient field x,y,z components
- `B1₊::Array{Float64, 3}`: the provided values of the B1 transmit field B1₊
- `B1₋::Vector{Array{ComplexF64, 2}}`: the provided values for two recieve coils for the B1 recieve field B1₋.


# Notes
-  See ../files/field_values.mat for example configuration format and headers to get fields
- The function uses interpolation and scaling to adjust the gradient fields based on input grid parameters.
- The nominal field values are used for testing or default scenarios where exact field values are not available.
"""
function initialise_fields(input_params, nominal_field, fields, lower_lims = [-0.21, -0.21, -0.21], upper_lims = [0.21, 0.21, 0.21])
    indexes_lims, coords_3d, num_points = get_index_lims(fields, lower_lims, upper_lims)
    input_grid, _ = get_required_grids(input_params, fields)

    B₀_x, B₀_y, B₀_z       = [fields["b0_x"][indexes_lims], fields["b0_y"][indexes_lims], fields["b0_z"][indexes_lims]]
    BGx_x, BGx_y, BGx_z    = [fields["BGx_x"][indexes_lims], fields["BGx_y"][indexes_lims], fields["BGx_z"][indexes_lims]]
    BGy_x, BGy_y, BGy_z    = [fields["BGy_x"][indexes_lims], fields["BGy_y"][indexes_lims], fields["BGy_z"][indexes_lims]]
    BGz_x, BGz_y, BGz_z    = [fields["BGz_x"][indexes_lims], fields["BGz_y"][indexes_lims], fields["BGz_z"][indexes_lims]]
    B1₊                    =  fields["B1p"][indexes_lims]
    B1₋_r1, B1₋_r2         = [fields["B1m_r1"][indexes_lims], fields["B1m_r2"][indexes_lims]]

    B₀_x, B₀_y, B₀_z       = [reshape_field(B₀_x, num_points), reshape_field(B₀_y, num_points), reshape_field(B₀_z, num_points)]
    BGx_x, BGx_y, BGx_z    = [reshape_field(BGx_x, num_points), reshape_field(BGx_y, num_points), reshape_field(BGx_z, num_points)]
    BGy_x, BGy_y, BGy_z    = [reshape_field(BGy_x, num_points), reshape_field(BGy_y, num_points), reshape_field(BGy_z, num_points)]
    BGz_x, BGz_y, BGz_z    = [reshape_field(BGz_x, num_points), reshape_field(BGz_y, num_points), reshape_field(BGz_z, num_points)]
    B1₊                    =  reshape_field(B1₊, num_points)
    B1₋_r1, B1₋_r2         = [reshape_field(B1₋_r1, num_points), reshape_field(B1₋_r2, num_points)]

    BGx_z_I = Interpolations.extrapolate(scale(interpolate(BGx_z, BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())  
    BGy_z_I = Interpolations.extrapolate(scale(interpolate(BGy_z, BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())  
    BGz_z_I = Interpolations.extrapolate(scale(interpolate(BGz_z, BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())  

    BGx_max, BGy_max, BGz_max  =  [Interpolations.gradient(BGx_z_I, 0, 0, 0)[1], Interpolations.gradient(BGy_z_I, 0, 0, 0)[2], Interpolations.gradient(BGz_z_I, 0, 0, 0)[3]]
    x_scale, y_scale, z_scale  =   1/BGx_max, 1/BGy_max, 1/BGz_max

    BGx_x_s, BGx_y_s, BGx_z_s    = [x_scale .* BGx_x, x_scale .* BGx_y, x_scale .* BGx_z]
    BGy_x_s, BGy_y_s, BGy_z_s    = [y_scale .* BGy_x, y_scale .* BGy_y, y_scale .* BGy_z]
    BGz_x_s, BGz_y_s, BGz_z_s    = [z_scale .* BGz_x, z_scale .* BGz_y, z_scale .* BGz_z]

    gz_z_nominal = coords_3d[:,:,:,3] 
    gy_z_nominal = coords_3d[:,:,:,2] 
    gx_z_nominal = coords_3d[:,:,:,1] 
    b0_z_nominal = ones(size(gx_z_nominal)) .* 0.099 #! change to whatever nominal field strength is neeeded
    nominal_xy_comps = zeros(size(gx_z_nominal))

    B₀_nominal =  [nominal_xy_comps, nominal_xy_comps, b0_z_nominal]
    BGx_nominal = [nominal_xy_comps, nominal_xy_comps, gx_z_nominal]
    BGy_nominal = [nominal_xy_comps, nominal_xy_comps, gy_z_nominal]
    BGz_nominal = [nominal_xy_comps, nominal_xy_comps, gz_z_nominal]
    B1₊_nominal = ones(size(gx_z_nominal))
    B1₋_nominal = [ones(size(gx_z_nominal)), ones(size(gx_z_nominal))]

    B₀  =  [B₀_x, B₀_y, B₀_z]
    BGx =  [BGx_x_s, BGx_y_s, BGx_z_s]
    BGy =  [BGy_x_s, BGy_y_s, BGy_z_s]
    BGz =  [BGz_x_s, BGz_y_s, BGz_z_s]
    B1₋ =  [B1₋_r1, B1₋_r2]

    if nominal_field
        return B₀_nominal, BGx_nominal, BGy_nominal, BGz_nominal, B1₊_nominal, B1₋_nominal
    else
        return B₀, BGx, BGy, BGz, B1₊, B1₋

    end
end


"""
    reshape_field(field_1d_vec, num_points)

Reshapes a 1D vector of field values into a 3D array and permutes its dimensions.

# Arguments
- `field_1d_vec::Vector{Float64}`: A 1-dimensional array containing field values. The length of this array should match the product of the dimensions specified in num_points.
- `num_points::Vector{Int64}`:I A tuple of three integers specifying the dimensions of the desired 3D array. The dimensions are used to reshape and permute the field_1d_vec.

# Returns
- A 3-dimensional array of field values with dimensions specified by num_points, permuted to the order (z, y, x).

# Example
julia
field_1d = 1:27
num_pts = (3, 3, 3)
reshaped_field = reshape_field(field_1d, num_pts)
"""
function reshape_field(field_1d_vec, num_points)
    field_out = permutedims(reshape(field_1d_vec, num_points[1], num_points[2], num_points[3]), (3,2,1))                       
    return field_out
end

"""
    calculate_ΔB₀(B₀, input_grid, input_params)

Calculates the field strength deviation ΔB₀ based on the provided magnetic field components and interpolation grid.

# Arguments
- `B₀::Vector{Array{Float64, 3}}`: A tuple of three 3D arrays representing the x, y, and z components of the main magnetic field. Each array should be of the same dimensions.
- `input_grid::3-element Vector{StepRangeLen{Float64, {Float64}, {Float64}, Int64}}`: A vector containing three 1D arrays representing the grid positions in the x, y, and z directions.
- `input_params::Dict{Any, Any}`: A dictionary containing various input parameters including paths to configuration files, nominal field values, and limits for the field initialization.

# Returns
- `ΔB₀::Float64`: The center frequency value ΔB₀. For 2D imaging, it returns the field strength at the specified slice_location in the z-direction. For 3D imaging, it returns the field strength at the origin (0, 0, 0).

"""
function calculate_ΔB₀(B₀, input_grid, input_params)
    B₀_net = sqrt.(B₀[1].^2 + B₀[2].^2 + B₀[3].^2)

    B₀_net_interp = Interpolations.extrapolate(scale(interpolate(B₀_net, BSpline(Linear())), input_grid[1], input_grid[2], input_grid[3]),Interpolations.Line())
    if input_params["imaging_type"] == "2D"
        ΔB₀ = B₀_net_interp(0,0,input_params["slice_location"])
    else
        ΔB₀ = B₀_net_interp(0,0,0)
    end

    return ΔB₀
end

"""
    get_spectral_response_3D_selection(B₀, ΔB₀)

Simulate a coil spectral response as a Lorentzian curve with a 42kHz FWHM bandwidth centered at ΔB₀
# Arguments
- `B₀::Vector{Array{Float64, 3}}`: A tuple of three 3D arrays representing the x, y, and z components of the main magnetic field. Each array should be of the same dimensions.
- `ΔB₀::Float64`: The center frequency value ΔB₀. For 2D imaging, it is the field strength at the specified slice_location in the z-direction. For 3D imaging, it is the field strength deviation at the origin (0, 0, 0).

# Returns
- `spectral_values::Array{Float64, 3}`: f(ω): spectral response curve values 

"""

function get_spectral_response_3D_selection(B₀, ΔB₀)
    # RF coil properties constants (HFSS (Ansys, Canonsburg, PA) )
    Cp=6021e-12;
    Cs=154.115e-12;
    R=0.03;
    L=222e-9;
    Z0=50;
    w(f) = 2*pi*f
    Zc(f)=R+1im*w(f)*L;
    Z(f) = 1/(1im*w(f)*Cs)+Zc(f)./(1+1im*w(f)*Cp.*Zc(f));
    Ic(f)= (1/(1+1im*w(f)*Cp.*Zc(f))*(2*sqrt(Z0))./(Z(f)+Z0));
    w0_values = 4e6:0.0001e6:4.5e6
    max_val = maximum([abs.(Ic(f)) for f in w0_values])
    new_func(f) = Ic(f)
    jcf = w0_values[argmax([abs.(new_func(f)) for f in w0_values])] 

    B₀_net = sqrt.(B₀[1].^2 + B₀[2].^2 + B₀[3].^2) # Tesla
    gyro =  42.577478518e6   # Hz/T
    net_b0_Hz  = B₀_net*gyro # Hz

    # center curve at the desired center frequency
    center_freq = ΔB₀*gyro # Hz 
    diff = center_freq-jcf 
    new_func_us(f) = Ic(f-diff)
    
    # now we want to ensure that at the CENTER of the scanner - spectral_response_curve(w(r)) = 1
    max_val = abs.(new_func_us(center_freq)) # this is NOT the maximum value, there is a point outside of this that is the best value

    final_curve(f) = abs.(new_func_us(f))/max_val
    # check that it returns 1 at the center #! but this will not be the maximum value, somewhere there will be a higher than one value
    final_curve(center_freq) == 1
    spectral_values = [final_curve(f) for f in net_b0_Hz]
    return spectral_values
end

"""
    makechunks(X::AbstractVector, n::Integer)

Splits an input vector X into n approximately equal-sized chunks.

# Arguments
- X: An AbstractVector (e.g., Array, Vector) to be split into chunks.
- n: An integer specifying the number of chunks to create. The input vector X will be divided into n parts.

# Returns
- A vector of n sub-vectors, where each sub-vector is a contiguous segment of the original vector X. The chunks are created such that they are approximately equal in size. The last chunk may be shorter if length(X) is not perfectly divisible by n.

# Example
julia
X = 1:10
n = 3
chunks = makechunks(X, n)
# chunks will be [[1, 2, 3], [4, 5, 6], [7, 8, 9, 10]]
"""
@views function makechunks(X::AbstractVector, Y::AbstractVector, n::Integer)
    c = length(X) ÷ n
    indexes_return = [X[1+c*k:(k == n-1 ? end : c*k+c)] for k = 0:n-1]
    TR_lines_return = [Y[1+c*k:(k == n-1 ? end : c*k+c)] for k = 0:n-1]
    return indexes_return, TR_lines_return
end
