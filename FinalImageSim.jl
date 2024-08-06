include("Initialize.jl")
include("src/Plotting.jl")
""" Inputs 
`result_dir::String`: Path to directory where the outputs from main.jl/main_slurm.jl have been saved. Ensure that only outputs from one experiment are there.
"""
result_dir = "/gpfs/data/johnsonplab/data/radhika_simulation/mri_inhomogeneous_field_spatial_encoding/FIG_7_2d_TSE_rc"

# Grab the files from the result directory and find experiment parameters
allfiles = readdir(result_dir, join=true)
example = load(allfiles[1])
input_params = example["params"]
sequence_params = example["sequence_params"]

# Plot the first 5 ms of the sequence 
plot_sequence(sequence_params, example["echo_time"], [1,5], true, true, false)

# Sum up all the data if spread across files 
final_kspace = zeros(Complex{Float64},size(example["k_space"]))
for f in allfiles
    out = load(f)
    ks = out["k_space"]
    final_kspace += ks
end

# Perform an IFFT on each coil image 
simulated_image_multicoil = zeros(size(final_kspace))

for i in range(1, size(final_kspace,1))
    simulated_image_multicoil[i,:,:] = abs.(ift2c(final_kspace[i,:,:]))
end

# Get final image by Root sum of squares 
simulated_image = get(rss_functions, input_params["imaging_type"], ()->error("Imaging type should be 2D or 3D"))(simulated_image_multicoil)

# Plot the image (If 3D, it will plot the middle slice)
get(plot_result_functions, input_params["imaging_type"], ()->error("Imaging type should be 2D or 3D"))(simulated_image)
