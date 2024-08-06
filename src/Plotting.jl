plotlyjs()
"""
Plots the simulated image `im_out`.

# Arguments
- `im_out::Array{<:Real}`: The image data to be plotted, a 2D array.
"""
function plot_image(im_out)
    heatmap(rot180(im_out),cmap=:grays,colorbar=false, axis=([],false),size=(300,300))
end

"""
Plots middle slice of the simulated image `im_out`.

# Arguments
- `im_out::Array{<:Real}`: The image data to be plotted, a 3D array.
"""

function plot_image_3d(im_out)
    slice_num = Int(ceil(size(im_out,3)/2))
    heatmap(im_out[:,:,slice_num],cmap=:grays,colorbar=false)

end

"""
    plot_sequence(sequence_params, echo_time, timing, plot_timed, plot_echo, plot_whole)

Plots various parameters from a sequence of measurements over time for a 2D sequence

# Arguments
- `sequence_params::`: An array containing parameters for each sequence, such as RF pulses and the spatial encoding gradient strengths.
- `echo_time::Float64`: The echo time to be plotted as a vertical line (can be obtained from the saved output of the simulation)
- `timing::Tuple{Float64, Float64}`: A tuple specifying the start and end times for the plot range (in ms)
- `plot_timed::Bool`: Flag indicating whether to plot a specific time range defined by `timing`.
- `plot_echo::Bool`: Flag indicating whether to indicate the location of the echo in the plot.
- `plot_whole::Bool`: Flag indicating whether to plot the entire range of data (not recommended, slow)

# Description
Basic sequence diagram 

"""
function plot_sequence(sequence_params, echo_time, timing, plot_timed, plot_echo, plot_whole)
    time_plot = zeros(1,length(sequence_params))
    gradx_plot = zeros(length(sequence_params))
    grady_plot = zeros(length(sequence_params))

    adcx_plot = fill(NaN32,length(sequence_params))
    adcy_plot = fill(NaN32,length(sequence_params))

    rf_plot = zeros(length(sequence_params))

    for i in range(1,length(time_plot))
        time_plot[i] = sequence_params[i].time*1000
        gradx_plot[i] = sequence_params[i].gradx
        grady_plot[i] = sequence_params[i].grady

        if  isfinite(sequence_params[i].indexKx) 
            adcx_plot[i] = Int(sequence_params[i].indexKx )           
        end
        if  isfinite(sequence_params[i].indexKy) 
            adcy_plot[i] = Int(sequence_params[i].indexKy)
        end

        if  isfinite(sequence_params[i].rf_) 
            rf_plot[i] = sequence_params[i].rf_*(pi/180)
            time_plot[i] = 0
        end
    end
    unique_y_lines = unique(adcy_plot)
    time_plot_ = cumsum(time_plot, dims = 2)
    if plot_whole
        index_to_plot_till = length(time_plot_)
        index_start = 1
    elseif plot_timed
        time_start = timing[1]
        time_end = timing[2]
        index_start = argmin(abs.(time_plot_ .- time_start))[2]
        index_to_plot_till = argmin(abs.(time_plot_ .- time_end))[2]
    else
        index_to_plot_till = findall(x->x==unique_y_lines[3], adcy_plot)[end]
        index_start = 1
    end

    fs = 12

    p1 = bar(time_plot_[1,index_start:index_to_plot_till],rf_plot[index_start:index_to_plot_till],ylabel="Flip angle[rad]",linewidth = 8, color=:green,xtickfontsize=fs,)
    p2 = plot(time_plot_[1,index_start:index_to_plot_till],gradx_plot[index_start:index_to_plot_till], ylabel = "BGx [T/m]", linewidth = 3,xtickfontsize=fs)
    p3 = plot(time_plot_[1,index_start:index_to_plot_till],grady_plot[index_start:index_to_plot_till] , ylabel = "BGy [T/m]",linewidth = 3,xtickfontsize=fs)
    p5 = plot(time_plot_[1,index_start:index_to_plot_till],adcx_plot[index_start:index_to_plot_till] , ylabel = "ADC_X",linewidth = 3,xtickfontsize=fs)
    p6 = plot(time_plot_[1,index_start:index_to_plot_till],adcy_plot[index_start:index_to_plot_till] , ylabel = "ADC_Y",xlabel = "Time[ms]",linewidth = 3,xtickfontsize=fs)

    plot(p1,p2,p3,p5,p6,layout = (5,1) ,size=(1100,800),legend=false,sharex=true, xlims=(time_plot_[index_start], time_plot_[index_to_plot_till]))

    if plot_echo
        echo_time_plot = zeros(length(time_plot_))
        ind_echo = argmin(abs.(time_plot_ .- echo_time*1000))[2]
        echo_time_plot[ind_echo] = 1;
        plot!(time_plot_[1,index_start:index_to_plot_till], echo_time_plot[index_start:index_to_plot_till],linewidth = 3,xtickfontsize=fs,sharex=true)
    end
end

"""
plot_sequence_3D(sequence_params, echo_time, timing, plot_timed, plot_echo, plot_whole)

Plots various parameters from a sequence of measurements over time for a 3D sequence

# Arguments
- `sequence_params::`: An array containing parameters for each sequence, such as RF pulses and the spatial encoding gradient strengths.
- `echo_time::Float64`: The echo time to be plotted as a vertical line (can be obtained from the saved output of the simulation)
- `timing::Tuple{Float64, Float64}`: A tuple specifying the start and end times for the plot range (in ms)
- `plot_timed::Bool`: Flag indicating whether to plot a specific time range defined by `timing`.
- `plot_echo::Bool`: Flag indicating whether to indicate the location of the echo in the plot.
- `plot_whole::Bool`: Flag indicating whether to plot the entire range of data (not recommended, slow)

# Description
Basic sequence diagram 

"""

function plot_sequence_3D(sequence_params, echo_time, timing, plot_timed, plot_echo, plot_whole)
    time_plot = zeros(1,length(sequence_params))
    gradx_plot = zeros(length(sequence_params))
    grady_plot = zeros(length(sequence_params))
    gradz_plot = zeros(length(sequence_params))

    adcx_plot = fill(NaN32,length(sequence_params))
    adcy_plot = fill(NaN32,length(sequence_params))
    adcz_plot = fill(NaN32,length(sequence_params))

    rf_plot = zeros(length(sequence_params))

    for i in range(1,length(time_plot))
        time_plot[i] = sequence_params[i].time*1000
        gradx_plot[i] = sequence_params[i].gradx
        grady_plot[i] = sequence_params[i].grady
        gradz_plot[i] = sequence_params[i].gradz

        if  isfinite(sequence_params[i].indexKx) 
            adcx_plot[i] = Int(sequence_params[i].indexKx )           
        end
        if  isfinite(sequence_params[i].indexKy) 
            adcy_plot[i] = Int(sequence_params[i].indexKy)
        end
        if  isfinite(sequence_params[i].indexKz) 
            adcz_plot[i] = Int(sequence_params[i].indexKz)
        end
        if  isfinite(sequence_params[i].rf_) 
            rf_plot[i] = sequence_params[i].rf_*(pi/180)
            time_plot[i] = 0
        end
    end
    unique_y_lines = unique(adcy_plot)
    time_plot_ = cumsum(time_plot, dims = 2)
    if plot_whole
        index_to_plot_till = length(time_plot_)
        index_start = 1
    elseif plot_timed
        time_start = timing[1]
        time_end = timing[2]
        index_start = argmin(abs.(time_plot_ .- time_start))[2]
        index_to_plot_till = argmin(abs.(time_plot_ .- time_end))[2]
    else
        index_to_plot_till = findall(x->x==unique_y_lines[3], adcy_plot)[end]
        index_start = 1
    end

    fs = 12

    p1 = bar(time_plot_[1,index_start:index_to_plot_till],rf_plot[index_start:index_to_plot_till],ylabel="Flip angle[rad]",linewidth = 8, color=:green,xtickfontsize=fs,)
    p2 = plot(time_plot_[1,index_start:index_to_plot_till],gradx_plot[index_start:index_to_plot_till], ylabel = "BGx [T/m]", linewidth = 3,xtickfontsize=fs)
    p3 = plot(time_plot_[1,index_start:index_to_plot_till],grady_plot[index_start:index_to_plot_till] , ylabel = "BGy [T/m]",linewidth = 3,xtickfontsize=fs)
    p4 = plot(time_plot_[1,index_start:index_to_plot_till],gradz_plot[index_start:index_to_plot_till] , ylabel = "BGz [T/m]",linewidth = 3,xtickfontsize=fs)
    p5 = plot(time_plot_[1,index_start:index_to_plot_till],adcx_plot[index_start:index_to_plot_till] , ylabel = "ADC_X",linewidth = 3,xtickfontsize=fs)
    p6 = plot(time_plot_[1,index_start:index_to_plot_till],adcy_plot[index_start:index_to_plot_till] , ylabel = "ADC_Y",xlabel = "Time[ms]",linewidth = 3,xtickfontsize=fs)
    p7 = plot(time_plot_[1,index_start:index_to_plot_till],adcz_plot[index_start:index_to_plot_till] , ylabel = "ADC_Z",xlabel = "Time[ms]",linewidth = 3,xtickfontsize=fs)

    plot(p1,p2,p3,p4,p5,p6,p7,layout = (7,1) ,size=(900,800),legend=false,sharex=true, xlims=(time_plot_[index_start], time_plot_[index_to_plot_till]))

    if plot_echo
        echo_time_plot = zeros(length(time_plot_))
        ind_echo = argmin(abs.(time_plot_ .- echo_time*1000))[2]
        echo_time_plot[ind_echo] = 1;
        plot!(time_plot_[1,index_start:index_to_plot_till], echo_time_plot[index_start:index_to_plot_till],linewidth = 3,xtickfontsize=fs,sharex=true)
    end
end

