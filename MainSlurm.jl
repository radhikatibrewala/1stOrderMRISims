include("Initialize.jl")

""" Inputs 
`config_file::String`: Path to .yaml file containing the field configuration files, the image input file and the simulation parameters
`filesave::String`: Path (without extension) for resulting MR values to be stored in
`num_jobs::Int`: "Number of bash slurm jobs this simulation is divided into
`array_num::Int`: "Array num of bash job"
"""

function run_sim(config_file, filesave, num_jobs, array_num)

    # 1. Make a directory to save the final outputs
    if ~isdir(dirname(filesave))
        mkdir(dirname(filesave))
    end
    save_loc = filesave*".jld2"
    println(save_loc)
    
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
    split_TR_line_indexes,  split_TR_lines = makechunks(new_TR_line_indexes, TR_lines, num_jobs) 
    current_TR_lines_indexes, current_TR_lines = split_TR_line_indexes[array_num], split_TR_lines[array_num] 
    println("These phase encode lines are being acquired:", current_TR_lines)

    # 6. Run the simulation for the sequence
    k_space_all = run_sequence(sequence_params, sim_vars_, input_params, current_TR_lines_indexes)

    # 7. Once the MR signal is calculated for each sequence timepoint, gather the 2D/3D k-space
    kspace_final = get(gather_kspace_functions, input_params["imaging_type"], ()->error("Invalid sequence name."))(k_space_all, input_params, sequence_params)

    # 8. Save the outputted files
    jldsave(save_loc; k_space=kspace_final, params = input_params, sequence_params = sequence_params, echo_time = echo_time)

end

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin
        "--array_num"
            help = "Slurm ARRAY number"
            arg_type = Int
        "--config_file_path"
            help = "Path to config file for this experiment"
            arg_type = String
        "--filesave"
            help = "Result filename without extension"
            arg_type = String
        "--num_jobs"
            help = "Number of jobs this simulation is divided into"
            arg_type = Int
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println(parsed_args)

    run_sim(parsed_args["config_file_path"], parsed_args["filesave"], parsed_args["num_jobs"], parsed_args["array_num"])
end

main()
