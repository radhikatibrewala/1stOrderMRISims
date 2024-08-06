using NaNStatistics

"""
    initialise_vectors(num_iters)

Initialize vectors used for gradient fields, RF pulses, time steps, and index values in a simulation.

# Arguments
- `num_iters::Int`: The number of iterations or steps for which the vectors should be initialized. This corresponds to the number of time points in the sequence.

# Returns
- `gradx_::Vector{Float32}`: A vector of zeros with num_iters elements representing the strength of the x-gradient at any iteration during the sequence
- `grady_::Vector{Float32}`: A vector of zeros with num_iters elements representing the strength of the y-gradient at any iteration during the sequence
- `rf_::Vector{Float32}`: A vector filled with NaN32 values, initialized to represent RF pulse values over num_iters.
- `time_::Vector{Float32}`: A vector filled with NaN32 values, initialized to represent time steps over num_iters.
- `indexKx_::Vector{Float32}`: A vector filled with NaN32 values, initialized to represent k-space index values in the x-direction.
- `indexKy_::Vector{Float32}`: A vector filled with NaN32 values, initialized to represent k-space index values in the y-direction.

# Description
This function initializes six vectors that are typically used in MRI simulations or other computational models involving gradient fields and RF pulses. The vectors gradx_ and grady_ are initialized to zeros and are intended to store the x and y components of the gradient field. The vectors rf_, time_, indexKx_, and indexKy_ are initialized with NaN32 values, serving as placeholders for RF pulse data, time steps, and k-space indices, respectively.

# Example
julia
num_iters = 1000
gradx_, grady_, rf_, time_, indexKx_, indexKy_ = initialise_vectors(num_iters)
"""
function initialise_vectors(num_iters)
    gradx_, grady_, rf_, time_, indexKx_, indexKy_ = zeros(Float32, num_iters), zeros(Float32, num_iters),fill(NaN32, num_iters),  fill(NaN32, num_iters), fill(NaN32, num_iters), fill(NaN32, num_iters)
    return gradx_, grady_, rf_, time_, indexKx_, indexKy_
end
"""
    initialise_vectors_3d(num_iters)

Initialize vectors used for gradient fields, RF pulses, time steps, and index values in a simulation.

# Arguments
- `num_iters::Int`: The number of iterations or steps for which the vectors should be initialized. This corresponds to the number of time points in the sequence.

# Returns
- `gradx_::Vector{Float32}`: A vector of zeros with num_iters elements representing the strength of the x-gradient at any iteration during the sequence
- `grady_::Vector{Float32}`: A vector of zeros with num_iters elements representing the strength of the y-gradient at any iteration during the sequence
- `gradz_::Vector{Float32}`: A vector of zeros with num_iters elements representing the strength of the z-gradient at any iteration during the sequence
- `rf_::Vector{Float32}`: A vector filled with NaN32 values, initialized to represent RF pulse values over num_iters.
- `time_::Vector{Float32}`: A vector filled with NaN32 values, initialized to represent time steps over num_iters.
- `indexKx_::Vector{Float32}`: A vector filled with NaN32 values, initialized to represent k-space index values in the x-direction.
- `indexKy_::Vector{Float32}`: A vector filled with NaN32 values, initialized to represent k-space index values in the y-direction.
- `indexKz_::Vector{Float32}`: A vector filled with NaN32 values, initialized to represent k-space index values in the z-direction.

# Description
This function initializes six vectors that are typically used in MRI simulations or other computational models involving gradient fields and RF pulses. The vectors gradx_ and grady_ are initialized to zeros and are intended to store the x and y components of the gradient field. The vectors rf_, time_, indexKx_, and indexKy_ are initialized with NaN32 values, serving as placeholders for RF pulse data, time steps, and k-space indices, respectively.

# Example
julia
num_iters = 1000
gradx_, grady_, gradz_, rf_, time_, indexKx_, indexKy_, indexKz_ = initialise_vectors_3d(num_iters)
"""

function initialise_vectors_3d(num_iters)
    gradx_, grady_, gradz_, rf_, time_, indexKx_, indexKy_, indexKz_= zeros(Float32, num_iters), zeros(Float32, num_iters),zeros(Float32, num_iters),fill(NaN32, num_iters),  fill(NaN32, num_iters), fill(NaN32, num_iters), fill(NaN32, num_iters),fill(NaN32, num_iters)
    return gradx_, grady_, gradz_, rf_, time_, indexKx_, indexKy_,indexKz_
end

"""
    calculate_field_strengths_2d(input_params)

Calculate the gradient field strengths required for a 2D imaging sequence based on input parameters.

# Arguments
- input_params::Dict{String, Float64}: A dictionary containing the necessary input parameters for the calculation:
  - "FOV_x"::Float64: The field of view in the x-direction, in meters.
  - "FOV_y"::Float64: The field of view in the y-direction, in meters.
  - "Δxy"::Float64: The pixel size (resolution) in meters.
  - "ADC_total_time"::Float64: The total acquisition time for the analog-to-digital converter (ADC), in seconds.

# Returns
- matrix_sizes::Tuple{Int, Int}: The matrix sizes in the x and y directions, respectively. Each value is the number of pixels along that axis.
- dwell_times::Tuple{Float64, Float64}: The dwell times for the x and y directions, respectively, in seconds. This is the time spent acquiring each pixel along each axis.
- gradient_field_strengths::Tuple{Float64, Float64}: The gradient field strengths in the x and y directions, respectively, in Tesla per meter (T/m). These values represent the magnetic field gradient per pixel.

# Description
This function calculates the matrix sizes, dwell times, and gradient field strengths required for a 2D imaging sequence based on the input field of view (FOV), pixel size, and total acquisition time. The gradient field strengths are computed using the Larmor equation with a gyromagnetic ratio (gamma_dash) of 42.577 MHz/T.

# Example
julia
input_params = Dict("FOV_x" => 0.42, "FOV_y" => 0.42, "Δxy" => 3e-3, "ADC_total_time" => 0.003)
matrix_sizes, dwell_times, gradient_field_strengths = calculate_field_strengths_2d(input_params)
"""

function calculate_field_strengths_2d(input_params)
    gamma_dash = 42.577e+6                                                    # MHz/T
    matrix_size_x = round(Int, input_params["FOV_x"]/input_params["Δxy"])     
    matrix_size_y = round(Int, input_params["FOV_y"]/input_params["Δxy"])     
    dwelltime_x = input_params["ADC_total_time"]/matrix_size_x                # seconds
    dwelltime_y = input_params["ADC_total_time"]/matrix_size_y                # seconds
    B_Gx_ =  1 / (gamma_dash * input_params["FOV_x"] * dwelltime_x )          # Tesla/m (per pixel)
    B_Gy_ =  1 / (gamma_dash * input_params["FOV_y"] * dwelltime_y)           # Tesla/m (per pixel)
    
    return [matrix_size_x, matrix_size_y], [dwelltime_x, dwelltime_y], [B_Gx_, B_Gy_]
end
"""
    calculate_field_strengths_3d(input_params)

Calculate the gradient field strengths required for a 3D imaging sequence based on input parameters.

# Arguments
- `input_params::Dict{String, Float64}`: A dictionary containing the necessary input parameters for the calculation:
  - "FOV_x"::Float64: The field of view in the x-direction, in meters.
  - "FOV_y"::Float64: The field of view in the y-direction, in meters.
  - "FOV_z"::Float64: The field of view in the z-direction, in meters.
  - "Δxy"::Float64: The pixel size (resolution) in meters.
  - "ADC_total_time"::Float64: The total acquisition time for the analog-to-digital converter (ADC), in seconds.

# Returns
- `matrix_sizes::Tuple{Int, Int, Int}`: The matrix sizes in the x, y and z directions, respectively. Each value is the number of pixels along that axis.
- `dwell_times::Tuple{Float64, Float64, Float64}`: The dwell times for the x, y and z directions, respectively, in seconds. This is the time spent acquiring each pixel along each axis.
- `gradient_field_strengths::Tuple{Float64, Float64, Float64}`: The gradient field strengths in the x, y and z directions, respectively, in Tesla per meter (T/m). These values represent the magnetic field gradient per pixel.

# Description
This function calculates the matrix sizes, dwell times, and gradient field strengths required for a 3D imaging sequence based on the input field of view (FOV), pixel size, and total acquisition time. The gradient field strengths are computed using the Larmor equation with a gyromagnetic ratio (gamma_dash) of 42.577 MHz/T.

# Example
julia
input_params = Dict("FOV_x" => 0.42, "FOV_y" => 0.42,  "FOV_z" => 0.42, "Δxyz" => 3e-3, "ADC_total_time" => 0.003)
matrix_sizes, dwell_times, gradient_field_strengths = calculate_field_strengths_3d(input_params)
"""
function calculate_field_strengths_3d(input_params)
    gamma_dash = 42.577e+6                                                    # MHz/T
    matrix_size_x = round(Int, input_params["FOV_x"]/input_params["Δxyz"])    
    matrix_size_y = round(Int, input_params["FOV_y"]/input_params["Δxyz"])    
    matrix_size_z = round(Int, input_params["FOV_z"]/input_params["Δxyz"])    
    dwelltime_x = input_params["ADC_total_time"]/matrix_size_x                # seconds
    dwelltime_y = input_params["ADC_total_time"]/matrix_size_y                # seconds
    dwelltime_z = input_params["ADC_total_time"]/matrix_size_z                # seconds

    B_Gx_ =  1 / (gamma_dash * input_params["FOV_x"] * dwelltime_x)           # Tesla/m (per pixel)
    B_Gy_ =  1 / (gamma_dash * input_params["FOV_y"] * dwelltime_y)           # Tesla/m (per pixel)
    B_Gz_ =  1 / (gamma_dash * input_params["FOV_z"] * dwelltime_z)           # Tesla/m (per pixel)

    return [matrix_size_x, matrix_size_y, matrix_size_z], [dwelltime_x, dwelltime_y, dwelltime_z], [B_Gx_, B_Gy_, B_Gz_]
end
