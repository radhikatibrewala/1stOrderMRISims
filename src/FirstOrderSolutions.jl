"""
    analytical_integral(ϕ, ∂ϕx, ∂ϕy, ∂ϕz, ω, ∂ωx, ∂ωy, ∂ωz, σ, grid_cube_size)

This function computes the first-order analytical solution to the 2D MRI signal equation, assuming a Gaussian Slice Profile

# Arguments
- `ϕ::Float64`: Accumulated phase in the grid-cube [radians]
- `∂ϕx::Float64`: Derivative of accumulated phase w.r.t x in the grid-cube [radians/m]
- `∂ϕy::Float64`: Derivative of accumulated phase w.r.t y in the grid-cube [radians/m]
-` ω::Float64`: Larmor frequency at the time of the excitation pulse in the grid-cube [tesla]
- `∂ωx::Float64`: Derivative of ω w.r.t to x in the grid-cube [tesla/m]
- `∂ωy::Float64`: Derivative of ω w.r.t to y in the grid-cube [tesla/m]
- `∂ωz::Float64`: Derivative of ω w.r.t to z in the grid-cube [tesla/m]
- `σ::Float64`: Spectral representation of the slice thickess [tesla]
- `grid_cube_size::Float64`:Simulation grid cube size [meters]

# Returns
- `integrated_value::Float64`: Integrated value of accumulated phase over grid-cube for 2D sequence

    # Notes
- Due to numerical instability, care needs to be taken while using ∂ωx, ∂ωy, ∂ωz values that tend to zero.
"""
function analytical_integral(ϕ, ∂ϕx, ∂ϕy, ∂ϕz, ω, ∂ωx, ∂ωy, ∂ωz, σ, grid_cube_size)
    d = grid_cube_size/2
    s = σ
    
    a1 = ∂ϕx*∂ωy*∂ωz* exp(1im * ( ∂ϕx*d + ∂ϕy*d + (∂ϕz/(∂ωz))*(ω - ∂ωx*d - ∂ωy*d)))
    a2 = ∂ϕx*∂ωy*∂ωz* exp(1im * ( ∂ϕx*d + ∂ϕz*d + (∂ϕy/(∂ωy))*(ω - ∂ωx*d - ∂ωz*d)))
    a3 = ∂ϕx*∂ωy*∂ωz* exp(1im * ( ∂ϕx*d - ∂ϕz*d + (∂ϕy/(∂ωy))*(ω - ∂ωx*d + ∂ωz*d)))
    a4 = ∂ϕx*∂ωy*∂ωz* exp(1im * ( ∂ϕx*d - ∂ϕy*d + (∂ϕz/(∂ωz))*(ω - ∂ωx*d + ∂ωy*d)))
    a5 = ∂ϕx*∂ωy*∂ωz* exp(1im * ( ∂ϕy*d - ∂ϕz*d + (∂ϕx/(∂ωx))*(ω - ∂ωy*d + ∂ωz*d)))
    a6 = ∂ϕx*∂ωy*∂ωz* exp(1im * (-∂ϕx*d + ∂ϕy*d + (∂ϕz/(∂ωz))*(ω + ∂ωx*d - ∂ωy*d)))
    a7 = ∂ϕx*∂ωy*∂ωz* exp(1im * (-∂ϕx*d + ∂ϕz*d + (∂ϕy/(∂ωy))*(ω + ∂ωx*d - ∂ωz*d)))
    a8 = ∂ϕx*∂ωy*∂ωz* exp(1im * (-∂ϕx*d - ∂ϕz*d + (∂ϕy/(∂ωy))*(ω + ∂ωx*d + ∂ωz*d)))
    a9 = ∂ϕx*∂ωy*∂ωz* exp(1im * (-∂ϕx*d - ∂ϕy*d + (∂ϕz/(∂ωz))*(ω + ∂ωx*d + ∂ωy*d)))
    a10 = ∂ϕy*∂ωx*∂ωz*exp(1im * ( ∂ϕx*d + ∂ϕy*d + (∂ϕz/(∂ωz))*(ω - ∂ωx*d - ∂ωy*d)))
    a11 = ∂ϕy*∂ωx*∂ωz*exp(1im * ( ∂ϕx*d - ∂ϕy*d + (∂ϕz/(∂ωz))*(ω - ∂ωx*d + ∂ωy*d)))
    a12 = ∂ϕy*∂ωx*∂ωz*exp(1im * ( ∂ϕy*d + ∂ϕz*d + (∂ϕx/(∂ωx))*(ω - ∂ωy*d - ∂ωz*d)))
    a13 = ∂ϕy*∂ωx*∂ωz*exp(1im * ( ∂ϕy*d - ∂ϕz*d + (∂ϕx/(∂ωx))*(ω - ∂ωy*d + ∂ωz*d)))
    a14 = ∂ϕy*∂ωx*∂ωz*exp(1im * (-∂ϕy*d + ∂ϕz*d + (∂ϕx/(∂ωx))*(ω + ∂ωy*d - ∂ωz*d)))
    a15 = ∂ϕy*∂ωx*∂ωz*exp(1im * (-∂ϕy*d - ∂ϕz*d + (∂ϕx/(∂ωx))*(ω + ∂ωy*d + ∂ωz*d)))
    a16 = ∂ϕy*∂ωx*∂ωz*exp(1im * (-∂ϕx*d + ∂ϕy*d + (∂ϕz/(∂ωz))*(ω + ∂ωx*d - ∂ωy*d)))
    a17 = ∂ϕy*∂ωx*∂ωz*exp(1im * (-∂ϕx*d - ∂ϕy*d + (∂ϕz/(∂ωz))*(ω + ∂ωx*d + ∂ωy*d)))
    a18 = ∂ϕz*∂ωx*∂ωy*exp(1im * ( ∂ϕx*d + ∂ϕz*d + (∂ϕy/(∂ωy))*(ω - ∂ωx*d - ∂ωz*d)))
    a19 = ∂ϕz*∂ωx*∂ωy*exp(1im * ( ∂ϕx*d - ∂ϕz*d + (∂ϕy/(∂ωy))*(ω - ∂ωx*d + ∂ωz*d)))
    a20 = ∂ϕz*∂ωx*∂ωy*exp(1im * ( ∂ϕy*d + ∂ϕz*d + (∂ϕx/(∂ωx))*(ω - ∂ωy*d - ∂ωz*d)))
    a21 = ∂ϕz*∂ωx*∂ωy*exp(1im * ( ∂ϕy*d - ∂ϕz*d + (∂ϕx/(∂ωx))*(ω - ∂ωy*d + ∂ωz*d)))
    a22 = ∂ϕz*∂ωx*∂ωy*exp(1im * (-∂ϕy*d + ∂ϕz*d + (∂ϕx/(∂ωx))*(ω + ∂ωy*d - ∂ωz*d)))
    a23 = ∂ϕz*∂ωx*∂ωy*exp(1im * (-∂ϕy*d - ∂ϕz*d + (∂ϕx/(∂ωx))*(ω + ∂ωy*d + ∂ωz*d)))
    a24 = ∂ϕz*∂ωx*∂ωy*exp(1im * (-∂ϕx*d + ∂ϕz*d + (∂ϕy/(∂ωy))*(ω + ∂ωx*d - ∂ωz*d)))
    a25 = ∂ϕz*∂ωx*∂ωy*exp(1im * (-∂ϕx*d - ∂ϕz*d + (∂ϕy/(∂ωy))*(ω + ∂ωx*d + ∂ωz*d)))

    p1 = 
    a1 * error_exponential_approximation(-∂ϕz*σ/(2*(∂ωz)),  ∂ωx*d/σ + ∂ωy*d/σ + ∂ωz*d/σ - ω/σ)  + a1 * error_exponential_approximation( ∂ϕz*σ/(2*(∂ωz)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ + ω/σ) -
    a2 * error_exponential_approximation(-∂ϕy*σ/(2*(∂ωy)),  ∂ωx*d/σ + ∂ωy*d/σ + ∂ωz*d/σ - ω/σ)  - a2 * error_exponential_approximation( ∂ϕy*σ/(2*(∂ωy)), -∂ωx*d/σ + ∂ωy*d/σ - ∂ωz*d/σ + ω/σ) -
    a3 * error_exponential_approximation( ∂ϕy*σ/(2*(∂ωy)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ + ω/σ)  - a3 * error_exponential_approximation(-∂ϕy*σ/(2*(∂ωy)),  ∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ - ω/σ) +
    a4 * error_exponential_approximation( ∂ϕz*σ/(2*(∂ωz)), -∂ωx*d/σ + ∂ωy*d/σ - ∂ωz*d/σ + ω/σ)  + a4 * error_exponential_approximation(-∂ϕz*σ/(2*(∂ωz)),  ∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ - ω/σ) -
    a5 * error_exponential_approximation( ∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ + ω/σ)  + a5 * error_exponential_approximation( ∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ + ω/σ) +
    a6 * error_exponential_approximation(-∂ϕz*σ/(2*(∂ωz)), -∂ωx*d/σ + ∂ωy*d/σ - ∂ωz*d/σ - ω/σ)  + a6 * error_exponential_approximation( ∂ϕz*σ/(2*(∂ωz)),  ∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ + ω/σ) -
    a7 * error_exponential_approximation( ∂ϕy*σ/(2*(∂ωy)),  ∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ + ω/σ)  - a7 * error_exponential_approximation(-∂ϕy*σ/(2*(∂ωy)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ - ω/σ) -
    a8 * error_exponential_approximation(-∂ϕy*σ/(2*(∂ωy)), -∂ωx*d/σ + ∂ωy*d/σ - ∂ωz*d/σ - ω/σ)  + a8 * error_exponential_approximation(-∂ϕy*σ/(2*(∂ωy)), -∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ - ω/σ) +
    a9 * error_exponential_approximation(-∂ϕz*σ/(2*(∂ωz)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ - ω/σ)  - a9 * error_exponential_approximation(-∂ϕz*σ/(2*(∂ωz)), -∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ - ω/σ) -
    a10* error_exponential_approximation(-∂ϕz*σ/(2*(∂ωz)),  ∂ωx*d/σ + ∂ωy*d/σ + ∂ωz*d/σ - ω/σ)  - a10* error_exponential_approximation( ∂ϕz*σ/(2*(∂ωz)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ + ω/σ)-
    a11* error_exponential_approximation( ∂ϕz*σ/(2*(∂ωz)), -∂ωx*d/σ + ∂ωy*d/σ - ∂ωz*d/σ + ω/σ)  - a11* error_exponential_approximation(-∂ϕz*σ/(2*(∂ωz)),  ∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ - ω/σ)+
    a12* error_exponential_approximation(-∂ϕx*σ/(2*(∂ωx)),  ∂ωx*d/σ + ∂ωy*d/σ + ∂ωz*d/σ - ω/σ)  + a12* error_exponential_approximation( ∂ϕx*σ/(2*(∂ωx)),  ∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ + ω/σ)+
    a13* error_exponential_approximation( ∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ + ω/σ)  + a13* error_exponential_approximation(-∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ + ∂ωy*d/σ - ∂ωz*d/σ - ω/σ) +
    a14* error_exponential_approximation( ∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ + ∂ωy*d/σ - ∂ωz*d/σ + ω/σ)  + a14* error_exponential_approximation(-∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ - ω/σ) -
    a15* error_exponential_approximation( ∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ + ∂ωy*d/σ + ∂ωz*d/σ + ω/σ)  - a15* error_exponential_approximation(-∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ - ω/σ)-
    a16* error_exponential_approximation(-∂ϕz*σ/(2*(∂ωz)), -∂ωx*d/σ + ∂ωy*d/σ - ∂ωz*d/σ - ω/σ)  - a16* error_exponential_approximation( ∂ϕz*σ/(2*(∂ωz)),  ∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ + ω/σ) -
    a17* error_exponential_approximation(-∂ϕz*σ/(2*(∂ωz)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ - ω/σ)  + a17* error_exponential_approximation(-∂ϕz*σ/(2*(∂ωz)), -∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ - ω/σ) +
    a18* error_exponential_approximation(-∂ϕy*σ/(2*(∂ωy)),  ∂ωx*d/σ + ∂ωy*d/σ + ∂ωz*d/σ - ω/σ)  + a18* error_exponential_approximation( ∂ϕy*σ/(2*(∂ωy)), -∂ωx*d/σ + ∂ωy*d/σ - ∂ωz*d/σ + ω/σ) +
    a19* error_exponential_approximation( ∂ϕy*σ/(2*(∂ωy)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ + ω/σ)  + a19* error_exponential_approximation(-∂ϕy*σ/(2*(∂ωy)),  ∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ - ω/σ) -
    a20* error_exponential_approximation(-∂ϕx*σ/(2*(∂ωx)),  ∂ωx*d/σ + ∂ωy*d/σ + ∂ωz*d/σ - ω/σ)  - a20* error_exponential_approximation( ∂ϕx*σ/(2*(∂ωx)),  ∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ + ω/σ) -
    a21* error_exponential_approximation(-∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ + ∂ωy*d/σ - ∂ωz*d/σ - ω/σ)  - a21* error_exponential_approximation( ∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ + ω/σ) -
    a22* error_exponential_approximation( ∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ + ∂ωy*d/σ - ∂ωz*d/σ + ω/σ)  - a22* error_exponential_approximation(-∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ - ω/σ) +
    a23* error_exponential_approximation( ∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ + ∂ωy*d/σ + ∂ωz*d/σ + ω/σ)  + a23* error_exponential_approximation(-∂ϕx*σ/(2*(∂ωx)), -∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ - ω/σ) +
    a24* error_exponential_approximation( ∂ϕy*σ/(2*(∂ωy)),  ∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ + ω/σ)  + a24* error_exponential_approximation(-∂ϕy*σ/(2*(∂ωy)), -∂ωx*d/σ - ∂ωy*d/σ + ∂ωz*d/σ - ω/σ)+
    a25* error_exponential_approximation(-∂ϕy*σ/(2*(∂ωy)), -∂ωx*d/σ + ∂ωy*d/σ - ∂ωz*d/σ - ω/σ)  - a25* error_exponential_approximation(-∂ϕy*σ/(2*(∂ωy)), -∂ωx*d/σ - ∂ωy*d/σ - ∂ωz*d/σ - ω/σ)
    
    denominator_term = (∂ϕy * ∂ωx - ∂ϕx * ∂ωy) * (∂ϕz * ∂ωx - ∂ϕx * ∂ωz) * (∂ϕz * ∂ωy - ∂ϕy * ∂ωz)
    denominator_term = adjust_inputs_varargs(denominator_term)[1]
    integrated_value = p1*exp(-1im * ϕ) * sqrt(π) * σ /(2 * denominator_term) 

    return integrated_value
end



"""Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables,  7.1.23 on p. 298  """
function error_exponential_approximation(A, B)
    complex_input = A*1im + B
    result = exp(-(A^2))*erf(complex_input)
    if ~isfinite(real(result)) && abs(angle(complex_input)) < 3pi/4
        result = (exp(-(A^2))-(exp((-2*A*B*1im)-(B^2))/((π^0.5)*(A*1im + B))))
    end
    return result
end


"""
analytical_integral_3d(ϕ, ∂ϕx, ∂ϕy, ∂ϕz,  grid_cube_size)

This function computes the first-order analytical solution to the 3D MRI signal equation

# Arguments
- `ϕ::Float64`: Accumulated phase in the grid-cube [radians]
- `∂ϕx::Float64`: Derivative of accumulated phase w.r.t x in the grid-cube [radians/m]
- `∂ϕy::Float64`: Derivative of accumulated phase w.r.t y in the grid-cube [radians/m]
- `σ::Float64`: Spectral representation of the slice thickess [tesla]
- `grid_cube_size::Float64`:Simulation grid cube size [meters]

# Returns
- `integrated_value::Float64`: Integrated value of accumulated phase over grid-cube for 3D sequence

# Notes
# Julia default is to calculate the normalized sinc function defined as sinc(x) = sin(pi*x)/(pi*x) 
# For our purposes, we need the  un-normalized sinc function defined as sinc(x) = sin(x)/(x)
# Therefore, we divide the input by pi to return the desired value
"""

function analytical_integral_3d(ϕ, ∂ϕx, ∂ϕy, ∂ϕz, grid_cube_size)
    d = grid_cube_size/2
    integrated_value = 8 * d^3 * exp(-1im * ϕ) * sinc(d * ∂ϕx / pi) * sinc(d * ∂ϕy / pi) * sinc(d * ∂ϕz / pi)
    return integrated_value
end

"""
    adjust_inputs_varargs(args...)

Adjusts input arguments by adding a small epsilon to any zero values.

# Arguments
- args...: A variable number of input arguments, each expected to be a numeric value or an array of numeric values.

# Returns
- Array: An array containing the adjusted input arguments. For each input, if a value is exactly 0.0, it is incremented by a small epsilon (1e-10). Non-zero values are returned unchanged.

# Description
This function is used to prevent issues that can arise when input arguments contain zero values, which might cause problems in subsequent computations.

# Example
julia
adjusted_values = adjust_inputs_varargs(0.0, 1.5, -2.3, 0.0)
# Returns: [1.0e-10, 1.5, -2.3, 1.0e-10]
"""
function adjust_inputs_varargs(args...)

    epsilon = 1e-10
    adjusted_args = map(x -> x == 0.0 ? x + epsilon : x, args)
    return adjusted_args
end
