using ProgressBars
using LoopVectorization
include("Structs.jl")
"""
    calculate_field_derivatives(sim_vars_, Gx_scale, Gy_scale, Gz_scale)

Calculate the total field and its gradients based on the gradient strengths at any given time point

# Arguments
- `sim_vars_::Struct`: An instance of `store_sim_vars` containing:
    - `B₀::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the main magnetic field - x,y,z components (B0x, B0y, B0z)
    - `BGx::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the x-gradient field  - x,y,z components (BGxx,BGxy,BGxz)
    - `BGy::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the y-gradient field - x,y,z components  (BGyx,BGyy,BGyz)
    - `BGz::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the z-gradient field  - x,y,z components (BGzx,BGzy,BGzz)
    - `B0_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the main magnetic field  - x,y,z components & w.r.t x-y-z directions (∂B0x/∂r), (∂B0y/∂r), (∂B0z/∂r)
    - `B_Gx_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the x-gradient - x,y,z components & w.r.t x-y-z directions (∂BGxx/∂r), (∂BGxy/∂r), (∂BGxz/∂r)
    - `B_Gy_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the y-gradient - x,y,z components & w.r.t x-y-z directions (∂BGyx/∂r), (∂BGyy/∂r), (∂BGyz/∂r)
    - `B_Gz_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the z-gradient - x,y,z components & w.r.t x-y-z directions (∂BGzx/∂r), (∂BGzy/∂r), (∂BGzz/∂r)
    - `ΔB₀::Float64`: The center frequency value (tesla) 
- `Gx_scale::Float64`: The scaling factor for the x-gradient field
- `Gy_scale::Float64`: The scaling factor for the y-gradient field 
- `Gz_scale::Float64`: The scaling factor for the z-gradient field 

# Returns
- `total_field::Array{Float64}`: The accumulated phase (Tesla)
- `∂total_field_x::Array{Float64}`: The gradient of the magnetic field potential with respect to the x direction.
- `∂total_field_y::Array{Float64}`: The gradient of the magnetic field potential with respect to the y direction.
- `∂total_field_z::Array{Float64}`: The gradient of the magnetic field potential with respect to the z direction.

# Description
The function computes the total field and its gradients based on the provided magnetic field components and the gradients.

#Notes:
- See Structs.jl for structure of sim_vars_
"""

function calculate_field_derivatives(sim_vars_, Gx_scale, Gy_scale, Gz_scale)
    total_field = zeros(size(sim_vars_.B₀[1]))
    ∂total_field_x = zeros(size(sim_vars_.B₀[1]))
    ∂total_field_y = zeros(size(sim_vars_.B₀[1]))
    ∂total_field_z = zeros(size(sim_vars_.B₀[1]))
    @inbounds @tturbo for i ∈ eachindex(total_field)
        total_field[i] = (sqrt(
              (sim_vars_.B₀[1][i] + Gx_scale * sim_vars_.BGx[1][i] + Gy_scale * sim_vars_.BGy[1][i] + Gz_scale * sim_vars_.BGz[1][i])^2
            + (sim_vars_.B₀[2][i] + Gx_scale * sim_vars_.BGx[2][i] + Gy_scale * sim_vars_.BGy[2][i] + Gz_scale * sim_vars_.BGz[2][i])^2
            + (sim_vars_.B₀[3][i] + Gx_scale * sim_vars_.BGx[3][i] + Gy_scale * sim_vars_.BGy[3][i] + Gz_scale * sim_vars_.BGz[3][i])^2
        ) - sim_vars_.ΔB₀) 
    end
    @inbounds @tturbo for i ∈ eachindex(∂total_field_x)
        ∂total_field_x[i] =(0.5*(Gx_scale*sim_vars_.BGx[1][i] + Gy_scale*sim_vars_.BGy[1][i] + Gz_scale*sim_vars_.BGz[1][i] + sim_vars_.B₀[1][i])*(2*Gx_scale*sim_vars_.∂BGx[1][1][i] + 2*Gy_scale*sim_vars_.∂BGy[1][1][i] + 2*Gz_scale*sim_vars_.∂BGz[1][1][i] + 2*sim_vars_.∂B₀[1][1][i]) + 0.5*(Gx_scale*sim_vars_.BGx[2][i] + Gy_scale*sim_vars_.BGy[2][i] + Gz_scale*sim_vars_.BGz[2][i] + sim_vars_.B₀[2][i])*(2*Gx_scale*sim_vars_.∂BGx[2][1][i] + 2*Gy_scale*sim_vars_.∂BGy[2][1][i] + 2*Gz_scale*sim_vars_.∂BGz[2][1][i] + 2*sim_vars_.∂B₀[2][1][i]) + 0.5*(Gx_scale*sim_vars_.BGx[3][i] + Gy_scale*sim_vars_.BGy[3][i] + Gz_scale*sim_vars_.BGz[3][i] + sim_vars_.B₀[3][i])*(2*Gx_scale*sim_vars_.∂BGx[3][1][i] + 2*Gy_scale*sim_vars_.∂BGy[3][1][i] + 2*Gz_scale*sim_vars_.∂BGz[3][1][i] + 2*sim_vars_.∂B₀[3][1][i]))*((Gx_scale*sim_vars_.BGx[1][i] + Gy_scale*sim_vars_.BGy[1][i] + Gz_scale*sim_vars_.BGz[1][i] + sim_vars_.B₀[1][i])^2 + (Gx_scale*sim_vars_.BGx[2][i] + Gy_scale*sim_vars_.BGy[2][i] + Gz_scale*sim_vars_.BGz[2][i] + sim_vars_.B₀[2][i])^2 + (Gx_scale*sim_vars_.BGx[3][i] + Gy_scale*sim_vars_.BGy[3][i] + Gz_scale*sim_vars_.BGz[3][i] + sim_vars_.B₀[3][i])^2)^(-0.5)
    end
    @inbounds @tturbo for i ∈ eachindex(∂total_field_y)
        ∂total_field_y[i] =(0.5*(Gx_scale*sim_vars_.BGx[1][i] + Gy_scale*sim_vars_.BGy[1][i] + Gz_scale*sim_vars_.BGz[1][i] + sim_vars_.B₀[1][i])*(2*Gx_scale*sim_vars_.∂BGx[1][2][i] + 2*Gy_scale*sim_vars_.∂BGy[1][2][i] + 2*Gz_scale*sim_vars_.∂BGz[1][2][i] + 2*sim_vars_.∂B₀[1][2][i]) + 0.5*(Gx_scale*sim_vars_.BGx[2][i] + Gy_scale*sim_vars_.BGy[2][i] + Gz_scale*sim_vars_.BGz[2][i] + sim_vars_.B₀[2][i])*(2*Gx_scale*sim_vars_.∂BGx[2][2][i] + 2*Gy_scale*sim_vars_.∂BGy[2][2][i] + 2*Gz_scale*sim_vars_.∂BGz[2][2][i] + 2*sim_vars_.∂B₀[2][2][i]) + 0.5*(Gx_scale*sim_vars_.BGx[3][i] + Gy_scale*sim_vars_.BGy[3][i] + Gz_scale*sim_vars_.BGz[3][i] + sim_vars_.B₀[3][i])*(2*Gx_scale*sim_vars_.∂BGx[3][2][i] + 2*Gy_scale*sim_vars_.∂BGy[3][2][i] + 2*Gz_scale*sim_vars_.∂BGz[3][2][i] + 2*sim_vars_.∂B₀[3][2][i]))*((Gx_scale*sim_vars_.BGx[1][i] + Gy_scale*sim_vars_.BGy[1][i] + Gz_scale*sim_vars_.BGz[1][i] + sim_vars_.B₀[1][i])^2 + (Gx_scale*sim_vars_.BGx[2][i] + Gy_scale*sim_vars_.BGy[2][i] + Gz_scale*sim_vars_.BGz[2][i] + sim_vars_.B₀[2][i])^2 + (Gx_scale*sim_vars_.BGx[3][i] + Gy_scale*sim_vars_.BGy[3][i] + Gz_scale*sim_vars_.BGz[3][i] + sim_vars_.B₀[3][i])^2)^(-0.5)
    end 
    @inbounds @tturbo for i ∈ eachindex(∂total_field_z)   
        ∂total_field_z[i] =(0.5*(Gx_scale*sim_vars_.BGx[1][i] + Gy_scale*sim_vars_.BGy[1][i] + Gz_scale*sim_vars_.BGz[1][i] + sim_vars_.B₀[1][i])*(2*Gx_scale*sim_vars_.∂BGx[1][3][i] + 2*Gy_scale*sim_vars_.∂BGy[1][3][i] + 2*Gz_scale*sim_vars_.∂BGz[1][3][i] + 2*sim_vars_.∂B₀[1][3][i]) + 0.5*(Gx_scale*sim_vars_.BGx[2][i] + Gy_scale*sim_vars_.BGy[2][i] + Gz_scale*sim_vars_.BGz[2][i] + sim_vars_.B₀[2][i])*(2*Gx_scale*sim_vars_.∂BGx[2][3][i] + 2*Gy_scale*sim_vars_.∂BGy[2][3][i] + 2*Gz_scale*sim_vars_.∂BGz[2][3][i] + 2*sim_vars_.∂B₀[2][3][i]) + 0.5*(Gx_scale*sim_vars_.BGx[3][i] + Gy_scale*sim_vars_.BGy[3][i] + Gz_scale*sim_vars_.BGz[3][i] + sim_vars_.B₀[3][i])*(2*Gx_scale*sim_vars_.∂BGx[3][3][i] + 2*Gy_scale*sim_vars_.∂BGy[3][3][i] + 2*Gz_scale*sim_vars_.∂BGz[3][3][i] + 2*sim_vars_.∂B₀[3][3][i]))*((Gx_scale*sim_vars_.BGx[1][i] + Gy_scale*sim_vars_.BGy[1][i] + Gz_scale*sim_vars_.BGz[1][i] + sim_vars_.B₀[1][i])^2 + (Gx_scale*sim_vars_.BGx[2][i] + Gy_scale*sim_vars_.BGy[2][i] + Gz_scale*sim_vars_.BGz[2][i] + sim_vars_.B₀[2][i])^2 + (Gx_scale*sim_vars_.BGx[3][i] + Gy_scale*sim_vars_.BGy[3][i] + Gz_scale*sim_vars_.BGz[3][i] + sim_vars_.B₀[3][i])^2)^(-0.5)
    end
    return total_field, ∂total_field_x, ∂total_field_y, ∂total_field_z
end
"""
    calculate_ϕ_∂ϕ(index, sequence_params, sim_vars_, ϕ, ∂ϕ, gx, gy, gz)

Compute the accumulated phase (`ϕ`) and its gradients (`∂ϕ`) based on the current sequence parameters and simulation variables.

# Arguments
- `index::Int`: The index of the current sequence parameter to be used for computation.
- `sequence_params::Vector{Sequence_TimeBlock}`: A vector of `Sequence_TimeBlock` structs containing:
  - `rf_::Float64`: The RF pulse angle (either 90 or 180 degrees).
  - `time::Float64`: The time parameter associated with the sequence.
  - `indexKx::T`, `indexKy::T`, `indexKz::T`: The indices for the k-space in x, y, and z directions (if applicable).
- `sim_vars_::store_sim_vars`: An instance of `store_sim_vars` containing:
  - Magnetic field components and their gradients and derivatives.
- `ϕ::Array{Float64}`: The array to hold the current phase.
- `∂ϕ::Array{Float64}`: An array of three arrays to hold the computed gradients of the phase.
- `gx::Float64`: The scaling factor for the x-gradient field
- `gy::Float64`: The scaling factor for the y-gradient field 
- `gz::Float64`: The scaling factor for the z-gradient field 


# Returns
- `ϕ::Array{Float64}`: The accumulated phase.
- `∂ϕ::Array{Array{Float64}}`: An array of three arrays containing the gradients of the phase with respect to x, y, and z directions.

# Description
The function accumulates the phase (`ϕ`) and its gradients (`∂ϕ`) based on the sequence parameters and simulation variables.
If a 90 degree RF pulse was applied, the accumulated phase is set to zero
If a 180 degree RF pulse was applied, the accumulated phase is negated
For all other instances, the complete applied field and its derivatives are calculated 
# Example
```julia
ϕ, ∂ϕ = calculate_ϕ_∂ϕ(index, sequence_params, sim_vars, ϕ, ∂ϕ, gx, gy, gz)
"""
function calculate_ϕ_∂ϕ!(index, sequence_params, sim_vars_, ϕ, ∂ϕ, gx, gy, gz)
    γ = 267.52218744e6 # radians/(Tesla * second)
    if sequence_params[index].rf_ == 90
        ϕ = zeros(Float64,size(ϕ))
        ∂ϕ[1] = zeros(Float64,size(ϕ)) #∂ϕ/∂x
        ∂ϕ[2] = zeros(Float64,size(ϕ)) #∂ϕ/∂y
        ∂ϕ[3] = zeros(Float64,size(ϕ)) #∂ϕ/∂z
    elseif sequence_params[index].rf_ == 180
        ϕ      = -ϕ
        ∂ϕ[1]  = -∂ϕ[1]
        ∂ϕ[2]  = -∂ϕ[2]
        ∂ϕ[3]  = -∂ϕ[3]
    else
        total_field, ∂total_field_x, ∂total_field_y, ∂total_field_z =  calculate_field_derivatives(sim_vars_, gx, gy, gz)
        ϕ     += γ * total_field * sequence_params[index].time      # radians 
        ∂ϕ[1] += γ * ∂total_field_x * sequence_params[index].time   # radians/m
        ∂ϕ[2] += γ * ∂total_field_y * sequence_params[index].time   # radians/m
        ∂ϕ[3] += γ * ∂total_field_z * sequence_params[index].time   # radians/m
        
    end
    return ϕ, ∂ϕ 
end
"""
Compute the sum integral of a given function with respect to specified parameters.
    run_zero_or_first_approximation
# Arguments
- `ϕ::Array{Float64}`: The accumulated phase.
- `∂ϕ::Array{Array{Float64}}`: An array of three arrays containing the gradients of the phase with respect to x, y, and z directions.
- `ω::Array{Float64}`: The ω (Larmor frequency at the time of the excitation pulse)
- `∂ω::Array{Array{Float64}}`: Tuple of arrays representing the partial derivatives of `ω` with respect to each spatial dimension.
- `σ::Float64`: Slice thickess for Gaussian slice profile (tesla)
- `grid_cube_size::Float64`: Size of the simulation grid cube
- `order::String`: Determines the type of integral computation. If "first", it uses `analytical_integral` else uses the 0th order equation

# Returns
- `Array{Complex{Float64}}`: An array containing the computed integral values.

# Description
This function computes the sum integral over a grid of values based on the given parameters. The integral computation depends on the `order` parameter:
- If `order` is "first", it uses the 1st order solution to the MR signal to compute the integral for each grid point.
- Otherwise, it computes the MR signal involving `ϕ`, `ω`, and `σ`, multiplied by the cube of `grid_cube_size`.

# Notes
- The function is optimized for performance using multi-threading.
- Ensure that the dimensions of the input arrays are compatible for element-wise operations.
"""

function run_zero_or_first_approximation(ϕ, ∂ϕ, ω, ∂ω, σ, grid_cube_size, order)
    values_all = zeros(Complex{Float64}, size(ϕ))
    if order == "first"
        @inbounds Threads.@threads  for i ∈ eachindex(values_all)
                values_all[i] = analytical_integral(ϕ[i], ∂ϕ[1][i], ∂ϕ[2][i], ∂ϕ[3][i], ω[i], ∂ω[1][i], ∂ω[2][i], ∂ω[3][i], σ, grid_cube_size)
        end
    else
        @inbounds Threads.@threads  for i ∈ eachindex(values_all)
                values_all[i] = exp(-1im*ϕ[i])*exp(-ω[i]^2/(σ)^2)*(grid_cube_size^3) 
        end
    end
    return values_all
end 
"""
    calculate_MR_signal(all_values, m₀, B1₋, B1₊)

Calculates the MR signal given the input values and magnetic field data.

# Arguments
- `all_values::Array{Complex{Float64}, 3}`: A 3D array representing the 1st or 0 order solutions for each location
- `m₀::Array{Float64, 3}`: A 3D array of m0 (input image)
- `B1₋::Array{Complex{Float64}, 4}`: A 4D array representing the receive coil values 
- `B1₊::Array{Float64, 3}`:  A 3D array representing the transmit coil values 

# Returns
- `Array{Float64, 3}`: A 3D array containing the MR signal for that timepoint

# Description
The function computes the MR signal for a 2D sequence for each element using the calculated 1st or 0th order values and the receive/transmit B1 coil profiles as well as m0.
Note: a flip angle of 5 degrees is assumed.
"""
function calculate_MR_signal(all_values, m₀, B1₋, B1₊)
    num_rec_coils = size(B1₋, 4)
    m₀_rep = repeat(m₀, 1, 1, 1, num_rec_coils)  
    B1₊_rep = repeat(B1₊, 1, 1, 1, num_rec_coils)  
    all_values_rep = repeat(all_values, 1, 1, 1, num_rec_coils)  

    MR_signal_all = zeros(Complex{Float64}, size(B1₊_rep))
    @inbounds Threads.@threads for i ∈ eachindex(MR_signal_all)
        MR_signal_all[i] = all_values_rep[i] * m₀_rep[i] * B1₋[i] * 5*pi/180 * B1₊_rep[i]
    end
    return vec(sum(MR_signal_all, dims = (1,2,3)))
end
"""
    calculate_MR_signal_3D(all_values, m₀, B1₋, B1₊)

Calculates the MR signal given the input values and magnetic field data.

# Arguments
- `all_values::Array{Complex{Float64}, 3}`: A 3D array representing the 1st or 0 order solutions for each location
- `m₀::Array{Float64, 3}`: A 3D array of m0 (input image)
- `B1₋::Array{Complex{Float64}, 4}`: A 4D array representing the receive coil values 
- `B1₊::Array{Float64, 3}`:  A 3D array representing the transmit coil values 
- `S::Array{Float64, 3}`:  A 3D array representing the spectral response curve values

# Returns
- `Array{Float64, 3}`: A 3D array containing the MR signal for that timepoint

# Description
The function computes the MR signal for a 3D sequence for each element using the calculated 1st or 0th order values and the receive/transmit B1 coil profiles, m0 and the spectral response curve.
    Note: a flip angle of 90 degrees is assumed.
"""
function calculate_MR_signal_3D(all_values, m₀, B1₋, B1₊, S)
    num_rec_coils = size(B1₋, 4)
    m₀_rep = repeat(m₀, 1, 1, 1, num_rec_coils)  
    B1₊_rep = repeat(B1₊, 1, 1, 1, num_rec_coils)  
    all_values_rep = repeat(all_values, 1, 1, 1, num_rec_coils)  
    S_rep = repeat(S, 1, 1, 1, num_rec_coils)  
    MR_signal_all = zeros(Complex{Float64}, size(B1₊_rep))
    @inbounds Threads.@threads  for i ∈ eachindex(MR_signal_all)
        MR_signal_all[i] = all_values_rep[i] * m₀_rep[i] * B1₋[i] * sin(B1₊_rep[i]*(pi/2)*S_rep[i]) 

    end
    return vec(sum(MR_signal_all, dims = (1,2,3)))
end
"""
    run_zero_or_first_approximation_3D(ϕ, ∂ϕ, grid_cube_size, order)

Computes either a zero-order or first-order approximation to calculate the MR signal 

# Arguments
- `ϕ::Array{Float64, 3}`: A 3D array representing the accumulated phase values
- `∂ϕ::Array{Array{Float64, 3}, 3}`: A tuple of three 3D arrays, each representing partial derivatives of the accumulated phase with respect to x-y-z
- `grid_cube_size::Float64`: The simulation grid cube size 
- `order::String`: Specifies the approximation order. Must be either `"first"` for first-order or "zero" for zero order.

# Returns
- `Array{Complex{Float64}, 3}`: A 3D array of the computed first/zero order values for location.

# Description
The function computes values based on the specified approximation order:
- If `order` is `"first"`, it performs a first-order approximation by evaluating an analytical integral over each voxel. 
- Otherwise, it computes the MR signal involving `ϕ` multiplied by the cube of `grid_cube_size`.

The result is a 3D array of complex values, where each element corresponds to the computed approximation value for the corresponding grid cube in the input field.

# Notes
- The function is optimized for performance using multi-threading.
- Ensure that the dimensions of the input arrays are compatible for element-wise operations.
"""
function run_zero_or_first_approximation_3D(ϕ, ∂ϕ, grid_cube_size, order)
    values_all = zeros(Complex{Float64}, size(ϕ))
    if order == "first"
        @inbounds Threads.@threads  for i ∈ eachindex(values_all) 
                values_all[i] = analytical_integral_3d(ϕ[i], ∂ϕ[1][i], ∂ϕ[2][i], ∂ϕ[3][i], grid_cube_size)
        end
    else
        @inbounds Threads.@threads  for i ∈ eachindex(values_all)
            values_all[i] = exp(-1im*ϕ[i])*(grid_cube_size^3)
        end
    end
    return values_all
end
"""
    calculate_kpsace_value(sequence_index, sequence_params, sim_vars_, ϕ, ∂ϕ, input_params)

Calculates the k-space value for a given sequence index, based on the provided simulation variables and input parameters.

# Arguments
- `sequence_index::Int`: The index of the current sequence in `sequence_params` for which the k-space value is calculated.
- `sequence_params::`: An array containing parameters for each sequence, such as RF pulses and the spatial encoding gradient strengths.
- `sim_vars_::store_sim_vars`: An instance of `store_sim_vars` containing:
    - `B₀::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the main magnetic field - x,y,z components (B0x, B0y, B0z)
    - `BGx::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the x-gradient field  - x,y,z components (BGxx,BGxy,BGxz)
    - `BGy::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the y-gradient field - x,y,z components  (BGyx,BGyy,BGyz)
    - `BGz::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the z-gradient field  - x,y,z components (BGzx,BGzy,BGzz)
    - `B0_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the main magnetic field  - x,y,z components & w.r.t x-y-z directions (∂B0x/∂r), (∂B0y/∂r), (∂B0z/∂r)
    - `B_Gx_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the x-gradient - x,y,z components & w.r.t x-y-z directions (∂BGxx/∂r), (∂BGxy/∂r), (∂BGxz/∂r)
    - `B_Gy_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the y-gradient - x,y,z components & w.r.t x-y-z directions (∂BGyx/∂r), (∂BGyy/∂r), (∂BGyz/∂r)
    - `B_Gz_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the z-gradient - x,y,z components & w.r.t x-y-z directions (∂BGzx/∂r), (∂BGzy/∂r), (∂BGzz/∂r)
    - `ΔB₀::Float64`: The center frequency value (tesla) - `ϕ::Array{Float64, 3}`: A 3D array representing the accrued phase.
- `∂ϕ::Array{Array{Float64, 3}, 3}`: A tuple of three 3D arrays, each representing partial derivatives of the accumulated phase with respect to x-y-z
- `input_params::Dict{String, T}`: A dictionary containing input parameters such as grid cube size (`"grid_cube_size"`), approximation order (`"order"`), and imaging type (2D/3D) (`"imaging_type"`).

# Returns
- `(Array{Complex{Float64}, 1}, Array{Float64, 3}, Array{Array{Float64, 3}, 3})`: 
  - `kval`: The calculated k-space value (complex number).
  - `ϕ::Array{Float64}`: The accumulated phase.
  - `∂ϕ::Array{Array{Float64}}`: An array of three arrays containing the gradients of the phase with respect to x, y, and z directions.

# Description
This function computes the k-space value based on the provided sequence and simulation parameters. Depending on the imaging type (`"3D"` or otherwise), it follows different paths:

- Overall steps:
  1. Extracts gradients (`gradx`, `grady`, `gradz`) for the specified sequence based on pre-calculated sequence parameters.
  2. Accumulates phase and phase derivates (ϕ and ∂ϕ) based on the gradients and RF pulses 
  3. Runs the zero- or first-order approximation based 
  4. Calculates the final MR signal in each grid-cube and sums to find the k-space value 


# Notes
- Ensure that `input_params` contains all required keys (`"grid_cube_size"`, `"order"`, `"imaging_type"`).
- The function assumes that `sim_vars_` contains the necessary fields for B1 fields, magnetization, and other relevant parameters.
-The B1₋ field is initialized and filled before computation, accommodating multi-channel reception. 
"""
function calculate_kpsace_value(sequence_index, sequence_params, sim_vars_, ϕ, ∂ϕ, input_params)
    grid_cube_size = input_params["grid_cube_size"]
    order = input_params["order"]
    imaging_type = input_params["imaging_type"]
    
    B1₋ = zeros(Complex{Float64}, size(sim_vars_.B1₋[1])[1],size(sim_vars_.B1₋[1])[2], size(sim_vars_.B1₋[1])[3],length(sim_vars_.B1₋))
    for i in range(1,size(B1₋,4))
        B1₋[:,:,:,i] = sim_vars_.B1₋[i][:,:,:]
    end

    if imaging_type == "3D"
        gx, gy, gz =  sequence_params[sequence_index].gradx, sequence_params[sequence_index].grady, sequence_params[sequence_index].gradz
        ϕ, ∂ϕ  = calculate_ϕ_∂ϕ!(sequence_index, sequence_params, sim_vars_, ϕ, ∂ϕ, gx, gy, gz)
        value =  run_zero_or_first_approximation_3D(ϕ, ∂ϕ, grid_cube_size, order) 
        kval =   calculate_MR_signal_3D(value, sim_vars_.m₀ , B1₋, sim_vars_.B1₊, sim_vars_.B1₊_Spectral)

    else
        gx, gy, gz =  sequence_params[sequence_index].gradx, sequence_params[sequence_index].grady, 0
        ϕ, ∂ϕ  = calculate_ϕ_∂ϕ!(sequence_index, sequence_params, sim_vars_, ϕ, ∂ϕ, gx, gy, gz)   
        value =  run_zero_or_first_approximation(ϕ, ∂ϕ, sim_vars_.ω, sim_vars_.∂ω, sim_vars_.σ, grid_cube_size, order)
        kval  =  calculate_MR_signal(value, sim_vars_.m₀,  B1₋, sim_vars_.B1₊)

    end
    
    return kval,  ϕ, ∂ϕ

end
"""
    run_sequence(sequence_params, sim_vars_, input_params, excite_pulse_timings)

Executes the entire MRI sequence simulation and computes MR signal measurements over time.

# Arguments
- `sequence_params::`: An array containing parameters for each sequence, such as RF pulses and the spatial encoding gradient strengths.
- `sim_vars_::store_sim_vars`: An instance of `store_sim_vars` containing:
    - `B₀::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the main magnetic field - x,y,z components (B0x, B0y, B0z)
    - `BGx::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the x-gradient field  - x,y,z components (BGxx,BGxy,BGxz)
    - `BGy::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the y-gradient field - x,y,z components  (BGyx,BGyy,BGyz)
    - `BGz::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the z-gradient field  - x,y,z components (BGzx,BGzy,BGzz)
    - `B0_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the main magnetic field  - x,y,z components & w.r.t x-y-z directions (∂B0x/∂r), (∂B0y/∂r), (∂B0z/∂r)
    - `B_Gx_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the x-gradient - x,y,z components & w.r.t x-y-z directions (∂BGxx/∂r), (∂BGxy/∂r), (∂BGxz/∂r)
    - `B_Gy_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the y-gradient - x,y,z components & w.r.t x-y-z directions (∂BGyx/∂r), (∂BGyy/∂r), (∂BGyz/∂r)
    - `B_Gz_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the z-gradient - x,y,z components & w.r.t x-y-z directions (∂BGzx/∂r), (∂BGzy/∂r), (∂BGzz/∂r)
    - `ΔB₀::Float64`: The center frequency value (tesla) - `input_params::Dict{String, T}`: A dictionary containing input parameters such as grid cube size (`"grid_cube_size"`) and other control variables.
- `excite_pulse_timings::Array{Tuple{Int, Int}, 1}`: An array of tuples, where each tuple specifies the start and end indices of the excitation pulse timings for each repetition time (TR).

# Returns
- `Array{Complex{Float64}, 2}`: A 2D array containing the computed MR signal measurements for each time point and reception coil. The dimensions are `(num_receive_coils, num_timepoints)`.

# Description
This function runs the entire MRI sequence by iterating over all excitation pulses defined by `excite_pulse_timings`. For each TR index:
1. It initializes the  `ϕ` phase and its partial derivatives `∂ϕ`.
2. For each time point within the current TR interval, it calculates the k-space value using `calculate_kpsace_value`.
3. The resulting MR signal measurement is stored in `mr_signal_measurements`.
4. A `garbage_handling` function is periodically called to manage memory and resources.

The function provides progress updates during the simulation and measures the time taken for each TR index. It returns the MR signal measurements collected over all time points.

# Note
- This function needs to run sequentially since the phase and its gradients are accumulated over time. By dividing into excite pulse timings, the function is already parallelized over different TRs (i.e when ϕ = 0)
"""

function run_sequence(sequence_params, sim_vars_, input_params, excite_pulse_timings)
    num_timepoints = length(sequence_params)
    
    mr_signal_measurements = [ zeros(Complex{Float64}, length(sim_vars_.B1₋), num_timepoints)][1]
    count = 0
    garbage_handling(count, input_params["grid_cube_size"])
    
    for tr ∈ ProgressBar(eachindex(excite_pulse_timings))
        println("On TR index ", tr, " of ", length(excite_pulse_timings))
        flush(stdout)

        ϕ = zeros(Float64,size(sim_vars_.B₀[1]))
        ∂ϕ = [zeros(Float64,size(sim_vars_.B₀[1])),zeros(Float64,size(sim_vars_.B₀[1])),zeros(Float64,size(sim_vars_.B₀[1]))]
    
        start_idx = excite_pulse_timings[tr][1]
        end_idx = excite_pulse_timings[tr][2]
    
        for index in range(start_idx, end_idx)
            mr_signal_measurements[:,index], ϕ, ∂ϕ =  calculate_kpsace_value(index, sequence_params, sim_vars_, ϕ, ∂ϕ, input_params)
            count += 1
            garbage_handling(count, input_params["grid_cube_size"])   
        end
    end
    flush(stdout)
    return mr_signal_measurements
end
"""
    garbage_handling(count, grid_cube_size)

Manages garbage collection during the simulation to optimize memory usage based on the current iteration count and grid cube size.

# Arguments
- `count::Int`: The current iteration count, used to determine when to trigger garbage collection.
- `grid_cube_size::Float64`: The size of the grid cube, used to determine the frequency of garbage collection.

# Description
This function is designed to manage memory during intensive simulations by triggering garbage collection (GC) at appropriate intervals. The frequency of garbage collection is determined based on the `grid_cube_size`:

- If `grid_cube_size` is less than `3/1000`, garbage collection is triggered every 10 iterations.
- Otherwise, garbage collection is triggered every 200 iterations.

The function checks if the current iteration count (`count`) is a multiple of the determined number of iterations (`num_iters_collect`). If it is, garbage collection is enabled, executed, and then disabled again to ensure that the simulation runs smoothly without excessive memory usage.

# Note
- The function uses Julia's `GC` module to enable, trigger, and disable garbage collection as needed.
- This function is typically called periodically during long-running simulations to manage memory effectively.
"""

function garbage_handling(count, grid_cube_size)
    if grid_cube_size < 3/1000
        flush(stdout)
        num_iters_collect = 10
    else
        num_iters_collect = 200
    end
    if mod(count, num_iters_collect) == 0
        GC.enable(true)
        GC.gc(true)
        GC.enable(false)
    end
end
