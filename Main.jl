include("Initialize.jl")
#! It is highly recommended to use the main_slurm.jl since this code uses threads and heavy parallization techniques and is primarily built to run over high performance computing resources 

""" Inputs 
`config_file::String`: Path to .yaml file containing the field configuration files, the image input file and the simulation parameters
`filesave::String`: Path (without extension) for resulting MR values to be stored in
"""

config_file = "examples/example_config.yaml"
filesave = "example_out"

# 1. Make a directory to save the final outputs
if ~isdir(dirname(filesave))
    mkdir(dirname(filesave))
end

# 2. Load the config file and get the selected order (1st or 0th order simulation)
config = YAML.load_file(config_file)
input_params = parse_config_sim(config)
order = input_params["order"]

# 3. Generate the simulation variables
sim_vars_ = generate_simulation_variables(input_params)

# 4. Obtain the sequence parameters
sequence_params, echo_time = get(sequence_functions, input_params["sequence_name"], ()->error("Invalid sequence name."))(input_params)

# 5. Calculate the indexes for the TRs in this simulation (simulate all TRs in this example, see main_slurm.jl for parallization over TRs)
new_TR_line_indexes, TR_lines = indexes_per_TR(sequence_params)
split_TR_line_indexes,  split_TR_lines = makechunks(new_TR_line_indexes, TR_lines, 1) #! 2 should be num jobs
current_TR_lines_indexes, current_TR_lines = split_TR_line_indexes[1], split_TR_lines[1] #! 2 should be array num
println("These phase encode lines are being acquired:", current_TR_lines)

# 6. Run the simulation for the sequence
k_space_all = run_sequence(sequence_params, sim_vars_, input_params, current_TR_lines_indexes)

# 7. Once the MR signal is calculated for each sequence timepoint, gather the 2D/3D k-space
kspace_final = get(gather_kspace_functions, input_params["imaging_type"], ()->error("Invalid sequence name."))(k_space_all, input_params, sequence_params)

# 8. Save the outputted files
save_loc = filesave*"_"*order*".jld2"
println("Saving and displaying final image...")
println(save_loc)
jldsave(save_loc; k_space=kspace_final, params = input_params, sequence_params = sequence_params, echo_time = echo_time)

# 9. Load the saved file and display the image after performing an IFFT on each coil and running a root sum of squares combination
perform_ifft_coil_combo_image(filesave)


