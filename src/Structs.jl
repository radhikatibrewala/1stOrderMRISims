using Parameters
"""
    struct Sequence_TimeBlock{T}

A structure representing a time block in a gradient echo sequence.

# Fields
- `gradx::T`: The strength of the x-gradient
- `grady::T`: The strength of the y-gradient
- `rf_::T`: The RF pulse value for this time block.
- `time::T`: The time step value.
- `indexKx::T`: The k-space index in the x-direction.
- `indexKy::T`: The k-space index in the y-direction.

# Description
The `Sequence_TimeBlock` struct is used to encapsulate the parameters for a single time block in a 2D sequence. Each field represents a different component of the sequence, including gradient fields, RF pulses, timing, and k-space indices.
"""
struct Sequence_TimeBlock{T} 
    gradx::T
    grady::T
    rf_::T
    time::T
    indexKx::T
    indexKy::T
end
"""
    struct Sequence_TimeBlock_3d{T}

A structure representing a time block in a gradient echo sequence.

# Fields
- `gradx::T`: The strength of the x-gradient
- `grady::T`: The strength of the y-gradient
- `gradz::T`: The strength of the z-gradient
- `rf_::T`: The RF pulse value for this time block.
- `time::T`: The time step value.
- `indexKx::T`: The k-space index in the x-direction.
- `indexKy::T`: The k-space index in the y-direction.
- `indexKz::T`: The k-space index in the z-direction.

# Description
The `Sequence_TimeBlock` struct is used to encapsulate the parameters for a single time block in a 3D sequence. Each field represents a different component of the sequence, including gradient fields, RF pulses, timing, and k-space indices.
"""
struct Sequence_TimeBlock_3d{T} 
    gradx::T
    grady::T
    gradz::T
    rf_::T
    time::T
    indexKx::T
    indexKy::T
    indexKz::T
end

"""
    struct store_sim_vars

A structure for storing simulation variables in a 2D imaging context.

# Fields
- `B0::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the main magnetic field - x,y,z components (B0x, B0y, B0z)
- `B_Gx::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the x-gradient field  - x,y,z components (BGxx,BGxy,BGxz)
- `B_Gy::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the y-gradient field - x,y,z components  (BGyx,BGyy,BGyz)
- `B_Gz::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the z-gradient field  - x,y,z components (BGzx,BGzy,BGzz)
- `B0_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the main magnetic field  - x,y,z components & w.r.t x-y-z directions (∂B0x/∂r), (∂B0y/∂r), (∂B0z/∂r)
- `B_Gx_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the x-gradient - x,y,z components & w.r.t x-y-z directions (∂BGxx/∂r), (∂BGxy/∂r), (∂BGxz/∂r)
- `B_Gy_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the y-gradient - x,y,z components & w.r.t x-y-z directions (∂BGyx/∂r), (∂BGyy/∂r), (∂BGyz/∂r)
- `B_Gz_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the z-gradient - x,y,z components & w.r.t x-y-z directions (∂BGzx/∂r), (∂BGzy/∂r), (∂BGzz/∂r)
- `ω::Array{Float64,3}`: A 3D array representing ω - sum of main magnetic field and the z-gradient field (ωx,ωy,ωz)
- `ω_d::Vector{Array{Float64, 3}}`: A vector of 3D arrays representing the gradient of ω w.r.t to x-y-z directions (∂ωx/∂r,∂ωy/∂r,∂ωz/∂r)
- `IMG::Array{Float64}`: the input image (m0)
- `B1M::Vector{Array{ComplexF64, 3}}`: A vector of 3D arrays representing the B1- field for each coil 
- `B1P::Array{Float64}`: An array representing the B1+ field 
- `ΔB0::Float64`: The center frequency value (tesla) 
- `σ::Float64`:The slice thickess, in tesla (slice_thickness [m] * BGz [T/m])
- `order::String`: "first" or "zero" order simulation

# Description
The `store_sim_vars` struct is used to store various variables related to a 2D imaging simulation, including magnetic field values, gradient fields, and image data. 
"""

struct store_sim_vars
    B₀::Vector{Array{Float64,3}}
    BGx::Vector{Array{Float64,3}}
    BGy::Vector{Array{Float64,3}}
    BGz::Vector{Array{Float64,3}}
    ∂B₀::Vector{Vector{Array{Float64, 3}}}
    ∂BGx::Vector{Vector{Array{Float64, 3}}}
    ∂BGy::Vector{Vector{Array{Float64, 3}}}
    ∂BGz::Vector{Vector{Array{Float64, 3}}}
    ω::Array{Float64,3}
    ∂ω::Vector{Array{Float64, 3}}
    m₀::Array{Float64}
    B1₋::Vector{Array{ComplexF64, 3}} #! Adjust 3 to number of receive coils 
    B1₊::Array{Float64}
    ΔB₀::Float64
    σ::Float64
    order::String
end


"""
    struct store_sim_vars

A structure for storing simulation variables in a 2D imaging context.

# Fields
- `B0::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the main magnetic field - x,y,z components (B0x, B0y, B0z)
- `B_Gx::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the x-gradient field  - x,y,z components (BGxx,BGxy,BGxz)
- `B_Gy::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the y-gradient field - x,y,z components  (BGyx,BGyy,BGyz)
- `B_Gz::Vector{Array{Float64,3}}`: A vector of 3D arrays representing the z-gradient field  - x,y,z components (BGzx,BGzy,BGzz)
- `B0_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the main magnetic field  - x,y,z components & w.r.t x-y-z directions (∂B0x/∂r), (∂B0y/∂r), (∂B0z/∂r)
- `B_Gx_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the x-gradient - x,y,z components & w.r.t x-y-z directions (∂BGxx/∂r), (∂BGxy/∂r), (∂BGxz/∂r)
- `B_Gy_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the y-gradient - x,y,z components & w.r.t x-y-z directions (∂BGyx/∂r), (∂BGyy/∂r), (∂BGyz/∂r)
- `B_Gz_d::Vector{Vector{Array{Float64, 3}}}`: A vector of vectors of 3D arrays representing the gradient of the z-gradient - x,y,z components & w.r.t x-y-z directions (∂BGzx/∂r), (∂BGzy/∂r), (∂BGzz/∂r)
- `IMG::Array{Float64}`: the input image (m0)
- `B1M::Vector{Array{ComplexF64, 3}}`: A vector of 3D arrays representing the B1- field for each coil 
- `B1P::Array{Float64}`: An array representing the B1+ field 
- `B1P_Spectral::Array{Float64}`: The weighting from the spectral response curve 

- `ΔB0::Float64`: The center frequency value (tesla) 
- `Δz_tesla::Float64`:The slice thickess, in tesla (slice_thickness [m] * BGz [T/m])
- `order::String`: "first" or "zero" order simulation

# Description
The `store_sim_vars_3d` struct is used to store various variables related to a 2D imaging simulation, including magnetic field values, gradient fields, and image data. 
    """
struct store_sim_vars_3d
    B₀::Vector{Array{Float64,3}}
    BGx::Vector{Array{Float64,3}}
    BGy::Vector{Array{Float64,3}}
    BGz::Vector{Array{Float64,3}}
    ∂B₀::Vector{Vector{Array{Float64, 3}}}
    ∂BGx::Vector{Vector{Array{Float64, 3}}}
    ∂BGy::Vector{Vector{Array{Float64, 3}}}
    ∂BGz::Vector{Vector{Array{Float64, 3}}}
    m₀::Array{Float64}
    B1₋::Vector{Array{ComplexF64, 3}} #! Adjust 3 to number of receive coils 
    B1₊::Array{Float64}
    B1₊_Spectral::Array{Float64}
    ΔB₀::Float64
    order::String
end
