"""
    load_image_file(image_path)

Loads and processes an image from the specified file path.

# Arguments
- `image_path::String`: The file path to the image file to be loaded.

# Returns
- `input_img::Array{Float64, 3}`: A 3D array representing the processed image data. The image is read from the file, permuted, and zeroed at the boundaries to handle extrapolation.

# Description
This function loads an image file using `matread`, processes it by reversing dimensions, permuting them, and ensuring that extrapolation values are set to zero along the boundaries. This helps in avoiding artifacts due to extrapolation.
"""

function load_image_file(image_path)
    input_img = permutedims(reverse(matread(image_path)["im"],dims=2), (1,3,2));
    input_img = permutedims(input_img, (2,3,1));

    input_img[  1,:,:] .= 0;                                 # ensure that extrapolation is zero in AP and LR direction
    input_img[end,:,:] .= 0;                                 # ensure that extrapolation is zero in AP and LR direction
    input_img[:,  1,:] .= 0;                                 # ensure that extrapolation is zero in AP and LR direction
    input_img[:,end,:] .= 0;                        

    return input_img
end


"""
get_m0(input_img)

Interpolates and transforms a 3D image to a new coordinate system.

# Arguments
- `input_img::Array{Float64, 3}`: A 3D array representing the input image data.

# Returns
- `img::Function`: A function that performs interpolation on the image data. The function maps 3D coordinates to image values with specific transformations.

# Description
This function creates an interpolated version of the input image using cubic B-spline interpolation. It then defines and returns a function `m0(x, y, z)` that provides the interpolated image values at given coordinates, with additional transformations applied.
"""
function get_m₀(input_img)
    img_posx = 0.56/1000 * (-size(input_img,1)/2+1:size(input_img,1)/2) # m   #!  change based on recon params
    img_posy = 0.56/1000 * (-size(input_img,2)/2+1:size(input_img,2)/2) # m   #!  change based on recon params
    img_posz = 3/1000 * (-size(input_img,3)/2+1:size(input_img,3)/2)    # m   #!  change based on recon params
    m₀ = Interpolations.extrapolate(scale(interpolate(input_img, BSpline(Cubic(Interpolations.Line(OnGrid())))), img_posx, img_posy, img_posz), Interpolations.Flat());
    return m₀
end

"""
    ift2c(kspace)

Performs an inverse 2D Fourier transform on the input k-space data.

# Arguments
- `kspace::Array{Complex{T}, 2}`: A 2D array representing the k-space data to be transformed, where `T` is a numeric type (e.g., `Float64`).

# Returns
- `out::Array{Complex{T}, 2}`: A 2D array containing the result of the inverse Fourier transform.

# Description
This function applies an inverse 2D Fourier transform to the provided k-space data. The transformation includes shifting the k-space data and then performing the inverse FFT to return the spatial domain representation.
"""
function ift2c(kspace)
    out = fftshift(ifft(ifftshift(kspace)))         
    return out
end

"""
    ft2c(im)

Performs a forward 2D Fourier transform on the input image data.

# Arguments
- `im::Array{Float64, 2}`: A 2D array representing the image data to be transformed.

# Returns
- `kout::Array{Complex{Float64}, 2}`: A 2D array containing the result of the forward Fourier transform.

# Description
This function applies a forward 2D Fourier transform to the provided image data. The transformation includes shifting the image data and then performing the FFT to return the k-space representation.
"""
function ft2c(im)
    kout = ifftshift(fft(fftshift(im)))         
    return kout
end
"""
    rss_coil_combine(im, dims)

Combines multiple coil images into a single image using root-sum-of-squares (RSS) method (2D image).

# Arguments
- `im::Array{Float64, N}`: 3-D array where the first dimension represents different coil images, and the other dimensions represent the 2D image.

# Returns
- `output::Array{Float64, 2}`: A 2D array representing the combined image.

# Description
This function combines multiple coil images into a single  image using the root-sum-of-squares method. The first dimension of the input array represents different coils. The output image is computed by summing the squares of the individual coil images and taking the square root of the result.
"""
function rss_coil_combine_2D(im)
    """ coil in first dim"""
    num_coils = size(im,1)
    output = zeros(size(im,2),size(im,3))
    for i in range(1,num_coils)
        output += im[i,:,:].^2
    end
    output = sqrt.(output)
    return output
end
"""
    rss_coil_combine_3D(im)

Combines multiple coil images into a single image using root-sum-of-squares (RSS) method (3D image).

# Arguments
- `im::Array{Float64, N}`: 4-D array where the first dimension represents different coil images, and the other dimensions represent the 3D image.

# Returns
- `output::Array{Float64, 2}`: A 3D array representing the combined image.

# Description
This function combines multiple coil images into a single  image using the root-sum-of-squares method. The first dimension of the input array represents different coils. The output image is computed by summing the squares of the individual coil images and taking the square root of the result.
"""
function rss_coil_combine_3D(im)
    """ coil in first dim"""
    num_coils = size(im,1)
    output = zeros(size(im,2),size(im,3),size(im,4))

    for i in range(1,num_coils)
        output += im[i,:,:].^2
    end
    output = sqrt.(output)
    return output
end
