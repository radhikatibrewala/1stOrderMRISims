include("Structs.jl")
include("InitializeSequences.jl")
using NaNStatistics

""" 2D Sequences:
    gradient_echo(input_params)
    spin_echo(input_params)
    single_shot_spin_echo_epi(input_params)
    single_shot_gradient_echo_epi(input_params)
    turbo_spin_echo(input_params)

Generate the selected sequence based on provided input parameters for an MRI simulation.

# Arguments
- `input_params::Dict{String, Float64}`: A dictionary containing the necessary input parameters for the calculation:
  - `"FOV_x"`: The field of view in the x-direction, in meters.
  - `"FOV_y"`: The field of view in the y-direction, in meters.
  - `"Δxy"`: The pixel size (resolution) in meters.
  - `"ADC_total_time"`: The total acquisition time for the analog-to-digital converter (ADC), in seconds.

# Returns
- `sequence_params_vector::Vector{Sequence_TimeBlock}`: A vector of `Sequence_TimeBlock` objects representing the gradient echo sequence. Each `Sequence_TimeBlock` contains the parameters for each time block in the sequence:
  - `gradx_`: The strengths of the x-gradient for each time block (T/m)
  - `grady_`: The strengths of the y-gradient for each time block (T/m)
  - `rf_`: The RF pulse value (90 or 180)
  - `time_`: The time step (time, in seconds)
  - `indexKx_`: The k-space index in the x-direction.
  - `indexKy_`: The k-space index in the y-direction.
- `TE::Float64`: The echo time, in seconds, which is the time from the center of the RF pulse to the peak of the echo signal.

# Description
This function generates the sequence parameters for a gradient echo MRI sequence. The process involves:
1. **Calculating Field Strengths**: The function first calculates the matrix sizes, dwell times, and gradient field strengths using `calculate_field_strengths_2d`.
2. **Timing Calculations**: The echo time (TE) and repetition time (TR) are calculated based on the matrix size and dwell time.
3. **Gradient and RF Pulse Initialization**: Vectors for gradient fields, RF pulses, time steps, and k-space indices are initialized using `initialise_vectors`.
4. **Sequence Construction**: The gradient and RF pulse values are filled in based on the timing and sequence requirements. The readout and rewinder gradients are configured, and k-space indices are set accordingly.
5. **Sequence Time Blocks**: The final sequence is returned as a vector of `Sequence_TimeBlock` objects, encapsulating all necessary parameters for each time step in the sequence.

# Examples for each:
```julia
input_params = Dict("FOV_x" => 0.42, "FOV_y" => 0.42, "Δxy" => 3e-3, "ADC_total_time" => 0.003)
sequence_params_vector, TE = gradient_echo(input_params)
sequence_params_vector, TE = spin_echo(input_params)
sequence_params_vector, TE = single_shot_spin_echo_epi(input_params)
sequence_params_vector, TE = single_shot_gradient_echo_epi(input_params)
sequence_params_vector, TE = turbo_spin_echo(input_params) # takes an additional input parameter: tsf = turbo spin factor

"""
# 2D gradient echo sequence 
function gradient_echo(input_params) 
    matrix_sizes, dwelltimes, field_strengths = calculate_field_strengths_2d(input_params)
    m = matrix_sizes[1]
    dt = dwelltimes[1]
    BGx_ = field_strengths[1]
    BGy_ = field_strengths[2]
    TE =  m*dt      
    TR =  (m+(m/2))*dt
    rewinder_time = TE/2
    readout_time = TE
    num_iters = Int(m*(1+(TR/dt))) #! add 1 for the 90
    gradx_, grady_, rf_, time_, indexKx_, indexKy_ = initialise_vectors(num_iters)
    time_[:] .= dt
    rf_[1:Int(TR/dt)+1:end] .= 90    
    
    for (i,index) in enumerate(findall(x->x == 90, rf_))
        rewinder_start = index + 1
        rewinder_end = Int(rewinder_start + rewinder_time/dt -1) 
        readout_start = rewinder_end + 1
        readout_end = Int(readout_start + readout_time/dt -1) 
        
        gradx_[rewinder_start:rewinder_end] .= -BGx_
        gradx_[readout_start:readout_end] .= BGx_
        grady_[rewinder_start:rewinder_end] .= (1 - (i - 1)/((m -1 )/2))*BGy_
        gradx_[readout_start:readout_end] .= BGx_
        indexKx_[readout_start:readout_end] = Int.(1:m)
        indexKy_[readout_start:readout_end] .= Int(i)
    end 
    sequence_params_vector = Array{Sequence_TimeBlock, 1}(undef, num_iters)
    for i in range(1,num_iters)
        sequence_params_vector[i] = Sequence_TimeBlock(gradx_[i], grady_[i], rf_[i], time_[i], indexKx_[i],indexKy_[i])
    end
    return sequence_params_vector, TE

end

# 2D Spin Echo Sequence
function spin_echo(input_params) 
    matrix_sizes, dwelltimes, field_strengths = calculate_field_strengths_2d(input_params)
    m = matrix_sizes[1]
    dt = dwelltimes[1]
    BGx_ = field_strengths[1]
    BGy_ = field_strengths[2]
    TE =  2*m*dt      
    TR =  (2*m*dt) + ((m/2)*dt)
    rewinder_time = TE/4
    readout_time = TE/2
    num_iters = Int(m*(2+(TR/dt))) #! 2 for the 90 + 180
    gradx_, grady_, rf_, time_, indexKx_, indexKy_ = initialise_vectors(num_iters)
    time_[:] .= dt
    rf_[1:Int(round(TR/dt)+2):end] .= 90    
    rf_[Int(round(TE/(2*dt))+2):Int(round(TR/dt)+2):end] .= 180


    for (i,index) in enumerate(findall(x->x == 180, rf_))
        rewinder_start = index + 1
        rewinder_end = Int(rewinder_start + rewinder_time/dt -1) # -1 because non-inclusive
        readout_start = rewinder_end + 1
        readout_end = Int(readout_start + readout_time/dt -1) # -1 because non-inclusive
        
        gradx_[rewinder_start:rewinder_end] .= -BGx_
        gradx_[readout_start:readout_end] .= BGx_
        grady_[rewinder_start:rewinder_end] .= (1 - (i - 1)/((m -1 )/2))*BGy_
        gradx_[readout_start:readout_end] .= BGx_
        indexKx_[readout_start:readout_end] = Int.(1:m)
        indexKy_[readout_start:readout_end] .= Int(i)
    end 
    sequence_params_vector = Array{Sequence_TimeBlock, 1}(undef, num_iters)
    for i in range(1,num_iters)
        sequence_params_vector[i] = Sequence_TimeBlock(gradx_[i], grady_[i], rf_[i], time_[i], indexKx_[i],indexKy_[i])
    end
    return sequence_params_vector, TE

end

# 2D Single Shot Spin Echo Sequence
function single_shot_spin_echo_epi(input_params)
    matrix_sizes, dwelltimes, field_strengths = calculate_field_strengths_2d(input_params)
    m = matrix_sizes[1]
    dt = dwelltimes[1]
    BGx_ = field_strengths[1]
    BGy_ = field_strengths[2]
    rewinder_iters = Int(m/2)
    TE = 2*dt*(Int(((m^2)/2)-1))
    num_iters = Int(Int(TE/(2*dt)) + 2 + (m/2) + (m^2))
    gradx_, grady_, rf_, time_, indexKx_, indexKy_ = initialise_vectors(num_iters)
    time_[:] .= dt
    rf_[1] = 90
    start_180 = Int(TE/(2*dt))+3
    rf_[start_180] = 180
    gradx_[Int(start_180 + 1): Int(start_180 + rewinder_iters)] .= -BGx_
    grady_[Int(start_180 + 1): Int(start_180 + rewinder_iters)] .= BGy_
    indexKx_[Int(start_180 + rewinder_iters)] = 1
    indexKy_[Int(start_180 + rewinder_iters)] = 1
    for row = 1:2:m
        offset = Int(start_180 + rewinder_iters) + 1 
        inds =  offset + (row-1) * m : offset +  (row-1) * m  + m - 2; # positive Gx
        gradx_[inds] .= BGx_
        grady_[inds[end] + 1] = -BGy_
        gradx_[inds .+ m] .= -BGx_
        indexKx_[inds[1]-1:inds[end]] = 1:m
        indexKx_[inds[1]-1 .+ m:inds[end] .+ m] = reverse(1:m)
        indexKy_[inds[1]-1:inds[end]] .= row
        indexKy_[inds[1]-1 .+ m:inds[end] .+ m] .= row+1
        if inds[end] + m + 1 <= num_iters
            grady_[inds[end] + m + 1] = -BGy_
        end
    end
        
    sequence_params_vector = Array{Sequence_TimeBlock, 1}(undef, num_iters) # or `Vector{Coords}(undef, x)`

    for i in range(1,num_iters)
        sequence_params_vector[i] = Sequence_TimeBlock(gradx_[i], grady_[i], rf_[i], time_[i], indexKx_[i],indexKy_[i])
    end
    return sequence_params_vector, TE
end


# 2D Single Shot Gradient Echo Sequence
function single_shot_gradient_echo_epi(input_params)
    matrix_sizes, dwelltimes, field_strengths = calculate_field_strengths_2d(input_params)
    m = matrix_sizes[1]
    dt = dwelltimes[1]
    BGx_ = field_strengths[1]
    BGy_ = field_strengths[2]
    rewinder_iters = Int(m/2)
    TE = dt*(Int(((m^2)/2)-1))
    num_iters = Int((m^2) + (m/2))
    gradx_, grady_, rf_, time_, indexKx_, indexKy_ = initialise_vectors(num_iters)
    time_[:] .= dt
    rf_[1] = 90
    gradx_[Int(1 + 1): Int(1 + rewinder_iters)] .= -BGx_
    grady_[Int(1 + 1): Int(1 + rewinder_iters)] .= BGy_
    indexKx_[Int(1 + rewinder_iters)] = 1
    indexKy_[Int(1 + rewinder_iters)] = 1
    for row = 1:2:m
        offset = Int(1 + rewinder_iters) + 1 
        inds =  offset + (row-1) * m : offset +  (row-1) * m  + m - 2; # positive Gx
        gradx_[inds] .= BGx_
        grady_[inds[end] + 1] = -BGy_
        gradx_[inds .+ m] .= -BGx_
        indexKx_[inds[1]-1:inds[end]] = 1:m
        indexKx_[inds[1]-1 .+ m:inds[end] .+ m] = reverse(1:m)
        indexKy_[inds[1]-1:inds[end]] .= row
        indexKy_[inds[1]-1 .+ m:inds[end] .+ m] .= row+1
        if inds[end] + m + 1 <= num_iters
            grady_[inds[end] + m + 1] = -BGy_
        end
    end
        
    sequence_params_vector = Array{Sequence_TimeBlock, 1}(undef, num_iters) # or `Vector{Coords}(undef, x)`

    for i in range(1,num_iters)
        sequence_params_vector[i] = Sequence_TimeBlock(gradx_[i], grady_[i], rf_[i], time_[i], indexKx_[i],indexKy_[i])
    end
    return sequence_params_vector, TE
end

# 2D Turbo Spin Echo Sequence
function turbo_spin_echo(input_params)
    matrix_sizes, dwelltimes, field_strengths = calculate_field_strengths_2d(input_params)
    m = matrix_sizes[1]
    dt = dwelltimes[1]
    BGx_ = field_strengths[1]
    BGy_ = field_strengths[2]
    tsf = input_params.tsf 
    TE =  2*m*dt
    TR = (m+(tsf*2*m)+tsf)*dt
    rewinder_time = TE/4
    readout_time = TE/2
    num_iters = Int((m/tsf)*(m+(tsf*2*m)+1+tsf))
    gradx_, grady_, rf_, time_, indexKx_, indexKy_ = initialise_vectors(num_iters)
    time_[:] .= dt
    rf_[1:Int(round(TR/dt)+1):end] .= 90   
    rf_[2+m:Int(round(TR/dt)+1):end] .= 180 # first refocusing pulse only


    for (i,index) in enumerate(findall(x->x==180, rf_))
        for tsf_line in range(1,tsf-1)
            rf_[index + tsf_line*((2*m)+1)] = 180 # subsequent refocusing pulses
        end
    end
    iter = 1
    ky_order = zeros(Int32,m) # seems ok
    for a in range(1,Int(m/tsf))
        for k in StepRange(a, Int(m/tsf), Int64(m))
            ky_order[iter] = k
            iter += 1
        end
    end
    for (i,index) in enumerate(findall(x->x==180, rf_))
        rewinder_start = index + 1
        rewinder_end = Int(rewinder_start + rewinder_time/dt -1) # -1 because non-inclusive
        readout_start = rewinder_end + 1
        readout_end = Int(readout_start + readout_time/dt -1) # -1 because non-inclusive
        rewinder_2_start = readout_end + 1
        rewinder_2_end =  Int(rewinder_2_start + rewinder_time/dt -1) # -1 because non-inclusive
        gradx_[rewinder_start:rewinder_end] .= -BGx_
        gradx_[readout_start:readout_end] .= BGx_
        gradx_[rewinder_2_start:rewinder_2_end] .= -BGx_
        grady_[rewinder_start:rewinder_end] .= (1 - (ky_order[i] - 1)/((m -1 )/2))*BGy_
        grady_[rewinder_2_start:rewinder_2_end] .= -(1 - (ky_order[i] - 1)/((m -1 )/2))*BGy_
        indexKx_[readout_start:readout_end] = 1:m
        indexKy_[readout_start:readout_end] .= ky_order[i]
    end      
    sequence_params_vector = Array{Sequence_TimeBlock, 1}(undef, num_iters) # or `Vector{Coords}(undef, x)`

    for i in range(1,num_iters)
        sequence_params_vector[i] = Sequence_TimeBlock(gradx_[i], grady_[i], rf_[i], time_[i], indexKx_[i],indexKy_[i])
    end
    return sequence_params_vector,TE
end
"""
    gradient_echo_3d(input_params)

Generate a 3D gradient echo sequence based on input parameters for an MRI simulation.

# Arguments
- `input_params::Dict{String, Float64}`: A dictionary containing the necessary input parameters for the calculation:
  - `"FOV_x"`: The field of view in the x-direction, in meters.
  - `"FOV_y"`: The field of view in the y-direction, in meters.
  - `"FOV_z"`: The field of view in the z-direction, in meters.
  - `"Δxy"`: The pixel size (resolution) in the x and y directions, in meters.
  - `"Δz"`: The pixel size (resolution) in the z-direction, in meters.
  - `"ADC_total_time"`: The total acquisition time for the analog-to-digital converter (ADC), in seconds.

# Returns
- `sequence_params_vector::Vector{Sequence_TimeBlock_3d}`: A vector of `Sequence_TimeBlock_3d` objects representing the 3D gradient echo sequence. Each `Sequence_TimeBlock_3d` contains the parameters for each time block in the sequence:
  - `gradx_`: The strengths of the x-gradient for each time block (T/m)
  - `grady_`: The strengths of the y-gradient for each time block (T/m)
  - `gradz_`: The strengths of the z-gradient for each time block (T/m)
  - `rf_`: The RF pulse value (90 or 180)
  - `time_`: The time step (time, in seconds)
  - `indexKx_`: The k-space index in the x-direction.
  - `indexKy_`: The k-space index in the y-direction.
  - `indexKz_`: The k-space index in the z-direction.
- `TE::Float64`: The echo time, in seconds, which is the time from the center of the RF pulse to the peak of the echo signal.

# Description
This function generates the sequence parameters for a 3D gradient echo MRI sequence. The process involves:
1. **Calculating Field Strengths**: The function first calculates the matrix sizes, dwell times, and gradient field strengths using `calculate_field_strengths_3d`.
2. **Timing Calculations**: The echo time (TE) and repetition time (TR) are calculated based on the matrix size and dwell time.
3. **Gradient and RF Pulse Initialization**: Vectors for gradient fields in x, y, and z directions, RF pulses, time steps, and k-space indices are initialized using `initialise_vectors_3d`.
4. **Sequence Construction**: The gradient and RF pulse values are filled in based on the timing and sequence requirements. The readout and rewinder gradients are configured, and k-space indices are set accordingly for each direction.
5. **Sequence Time Blocks**: The final sequence is returned as a vector of `Sequence_TimeBlock_3d` objects, encapsulating all necessary parameters for each time step in the sequence.

# Example
```julia
input_params = Dict("FOV_x" => 0.42, "FOV_y" => 0.42, "FOV_z" => 0.42, "Δxyz" => 3e-3, "ADC_total_time" => 0.003)
sequence_params_vector, TE = gradient_echo_3d(input_params)
"""

#3D Gradient Echo Sequence
function gradient_echo_3d(input_params) 
    matrix_sizes, dwelltimes, field_strengths = calculate_field_strengths_3d(input_params)
    m = matrix_sizes[1]
    dt = dwelltimes[1]
    BGx_ = field_strengths[1]
    BGy_ = field_strengths[2]
    BGz_ = field_strengths[3]
    TE =  m*dt      
    TR =  (m+(m/2))*dt
    rewinder_time = TE/2
    readout_time = TE
    num_iters = Int(m + (m/2) + 1)*m*m # m readout, m/2 rewinder, +1 for 90
    gradx_, grady_, gradz_, rf_, time_, indexKx_, indexKy_, indexKz_= initialise_vectors_3d(num_iters)
    time_[:] .= dt
    rf_[1:Int(TR/dt)+1:end] .= 90  
    
    indexes_90 = findall(x->x == 90, rf_)
    weights_pe1 = zeros(length(indexes_90))
    weights_pe2 = zeros(length(indexes_90))
    for i = 1:m
        positions = i*m - (m-1) : i*m 
        weights_pe1[positions] .= i
        weights_pe2[positions] = 1:m
    end
    for (i,index) in enumerate(indexes_90)
        rewinder_start = index + 1
        rewinder_end = Int(rewinder_start + rewinder_time/dt -1) # -1 because non-inclusive
        readout_start = rewinder_end + 1
        readout_end = Int(readout_start + readout_time/dt -1) # -1 because non-inclusive
        
        gradx_[rewinder_start:rewinder_end] .= -BGx_
        grady_[rewinder_start:rewinder_end] .= -(1 - (weights_pe2[i] - 1)/((m -1 )/2))*BGy_
        gradz_[rewinder_start:rewinder_end] .= -(1 - (weights_pe1[i] - 1)/((m -1 )/2))*BGz_

        gradx_[readout_start:readout_end] .= BGx_

        indexKx_[readout_start:readout_end] = Int.(1:m)
        indexKy_[readout_start:readout_end] .= Int(weights_pe2[i])
        indexKz_[readout_start:readout_end] .= Int(weights_pe1[i])

    end 
    sequence_params_vector = Array{Sequence_TimeBlock_3d, 1}(undef, num_iters)
    for i in range(1,num_iters)
        sequence_params_vector[i] = Sequence_TimeBlock_3d(gradx_[i], grady_[i], gradz_[i], rf_[i], time_[i], indexKx_[i],indexKy_[i],indexKz_[i])
    end
    return sequence_params_vector, TE

end

