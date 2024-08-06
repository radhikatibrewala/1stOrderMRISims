"""
    allocate_kspace(measurements_all, input_params, sequence_params)
    allocate_kspace_3D(measurements_all, input_params, sequence_params)

Assembles the k-space data from the given measurements from simulation based on the provided sequence parameters.

# Arguments
- `measurements_all::Array{Complex{Float64}}`: A 2D array where each column has the calculated MRI signal for a specific time point. The dimensions should be `(num_coils, num_timepoints)`.
- `input_params::Dict{String, Any}`: A dictionary containing input parameters. It must (at least) include:
  - `"FOV_x"`: The field of view in the x-direction.
  - `"FOV_y"`: The field of view in the x-direction.
  - `"Δxy"`: Resolution.
- `sequence_params::Array{Sequence_TimeBlock, 1}`: An array of `Sequence_TimeBlock` structs where each element contains:
  - `indexKx::Int`: The k-space index in the x-direction.
  - `indexKy::Int`: The k-space index in the y-direction.

# Returns
- `kspace_return::Array{Complex{Float64}}`: A 3D array representing the assembled k-space data with dimensions `(num_measurements, m, m)`, where `m` is calculated as `Int(input_params["FOV_x"]/input_params["Δxy"])`.

# Description
The function iterates over each sequence parameter and maps the corresponding measurement values to the appropriate k-space locations based on the indices provided in `sequence_params`. It constructs a 3D array (`kspace_return`) where the dimensions correspond to the number of measurements and the spatial resolution defined by the field of view and pixel size.
Use this function once you the entire simulation has run to assemble the signal values into a k-space array (2D/3D)
"""

function allocate_kspace(measurements_all, input_params, sequence_params)
    m = Int(input_params["FOV_x"]/input_params["Δxy"])
    n = Int(input_params["FOV_y"]/input_params["Δxy"])

    kspace_return = zeros(Complex{Float64},(size(measurements_all,1),m,n))
    num_iters = length(sequence_params)
    for i in range(1,num_iters)
        kx = sequence_params[i].indexKx
        ky = sequence_params[i].indexKy
        for j in range(1,size(measurements_all,1))
            if isfinite(kx) & isfinite(ky)
                kspace_return[j, Int(ky),Int(kx)] = measurements_all[j,i]
            end
        end
    end
    return kspace_return 
end

function allocate_kspace_3D(measurements_all, input_params, sequence_params)
    m = Int(input_params["FOV_x"]/input_params["Δxyz"])
    n = Int(input_params["FOV_y"]/input_params["Δxyz"])
    o = Int(input_params["FOV_z"]/input_params["Δxyz"])

    kspace_return = zeros(Complex{Float64},(size(measurements_all,1),m,n,o))
    num_iters = length(sequence_params)
    for i in range(1,num_iters)
        kx = sequence_params[i].indexKx
        ky = sequence_params[i].indexKy
        kz = sequence_params[i].indexKz
        for j in range(1,size(measurements_all,1))
            if isfinite(kx) & isfinite(ky) & isfinite(kz)
                kspace_return[j,Int(kx),Int(ky),Int(kz)] = measurements_all[j,i]
            end
        end
    end
    return kspace_return
end
"""
    indexes_per_TR(sequence_params)

Determines the start and end indices of each repetition time (TR) interval based on the RF pulse angle of 90 degrees in the given sequence parameters.

# Arguments
- `sequence_params::Array{Sequence_TimeBlock, 1}`: An array of `Sequence_TimeBlock` structs where each element contains:
  - `rf_::Float64`: The RF pulse angle in degrees. This function specifically looks for values of 90 degrees.

# Returns
- `vec_list::Array{Vector{Int}, 1}`: A list of vectors where each vector represents the start and end indices of each TR interval. The indices correspond to the positions of the 90-degree RF pulses in `sequence_params`.

# Description
The function identifies the indices in `sequence_params` where the RF pulse angle is 90 degrees. It then groups these indices into intervals, with each interval spanning from one 90-degree RF pulse to the next. The resulting list (`vec_list`) contains vectors specifying the start and end indices of each TR interval, including the final interval extending to the end of the sequence.

"""
function indexes_per_TR(sequence_params)
    ind_90 = zeros(Int, length(sequence_params))
    for i in range(1,length(sequence_params))
        if sequence_params[i].rf_ == 90
            ind_90[i] = 1
        end
    end

    ind_90 = findall(x->x==1, ind_90)
    vec_list = [[] for _ = 1:length(ind_90)]
    tr_lines_list = zeros(Int,length(ind_90))

    for i in range(1,length(ind_90)-1)
        vec_list[i] = [ind_90[i],ind_90[i+1]-1]
        tr_lines_list[i] = i
    end
    vec_list[end] = [ind_90[end], length(sequence_params)]
    tr_lines_list[end] = length(ind_90)
    return vec_list, tr_lines_list
end