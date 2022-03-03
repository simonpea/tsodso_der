using Pkg
Pkg.activate("RunDSO")

println("Loading required packages: ProgressMeter, CSV, DataFrames, Dates, JuMP, Gurobi.")
using ProgressMeter
using CSV
using DataFrames
using Dates
using JuMP
using Gurobi
println("Loading of packages complete.")


###################
# Terminal ########
###################

    usedayperiod = false
    startday = 1
    lastday = 2
    starthour = 1
    lasthour = 8760
	x_ag = 0.3
	x_re = 0.3
	x_ind = 0.4
	s_ag = 0.05
	s_re = 0.05
	s_ind = 0.1
    workhours = 9:16

    ### Create Dict for average RD Price
    datapath = "data/"
    inputpath = "Processed_output/BAU/"
    output = "output/"
    outputpath = joinpath(output, "output_LoadShift/")

    isworkhour = Dict()
    for i in 0:24
        if i in workhours isworkhour[i] = 1
        else isworkhour[i] = 0
        end
    end
    if usedayperiod == true
        timeperiod = string((startday-1)*24+1)*":"*string((lastday*24))
    else
        timeperiod = string(starthour)*":"*string(lasthour)
    end

    timeperiod_str = replace(string(timeperiod), ":" => "-")


    println()
    # println("Case >" * case * "< selected.")
    println("Calculating for new load profiles.")
    println()

###################
# Run CM ##########
###################

    include("model/src/Sesam_ptdf_LoadShift.jl")

    using .Sesam

###################
# Data ############
###################
# Read nodes and load from .csv
    println("Reading data from csv...")

	nodes_df = CSV.read(joinpath(datapath, "nodes.csv"), DataFrame)

	load_df = CSV.read(joinpath(datapath, "load.csv"), DataFrame)

	exchange_df = CSV.read(joinpath(datapath, "exchange.csv"), DataFrame)

    nodal_hourlycost_df = CSV.read(joinpath(inputpath, "hourly_nodal_RD_cost.csv"), DataFrame)

    nodal_hourlyvolume_df = CSV.read(joinpath(inputpath, "hourly_nodal_RD_volume.csv"), DataFrame)

    marketprice_df = CSV.read(joinpath(inputpath, "ED_price.csv"), DataFrame)
                                    
    println("Reading of data complete.")
    println()
    println("Creating dictionaries for costs, volumes and marketprice...")


    nod_vol_dict = Dict((nodal_hourlyvolume_df[i, 1] => Dict((j => nodal_hourlyvolume_df[i, j+1]) for j in 1:8760)) for i in 1:nrow(nodal_hourlyvolume_df))
    nod_cost_dict = Dict((nodal_hourlycost_df[i, 1] => Dict((j => nodal_hourlycost_df[i, j+1]) for j in 1:8760)) for i in 1:nrow(nodal_hourlycost_df))
    av_el_price_dict = Dict((i => (sum(marketprice_df[(i-1)*24+1:i*24, 2]))/24) for i in 1:365)
###################
# Model setup #####
###################

    println("Preparing nodes...")

    # Create sets
        # Create set of nodes
        nodes = Nodes(nodes_df, load_df, exchange_df, 1000, 1)

        N = nodes.id

    println("Preparation complete.")
    println()

    session_time = Dates.now()
    session = Dates.format(session_time, "yyyy-mm-dd__HH-MM-SS")*"_newload"


Newload_df = DataFrame(node = N)
LoadShift_df = DataFrame(hour = Int64[], node=Int64[], OriginalLoad=Float64[], NewLoad=Float64[], ShiftUp=Float64[], ShiftDn=Float64[])
CSV.write(joinpath(outputpath, "Loadshift.csv"), LoadShift_df)

prog = Progress(length(dayslicer(timeperiod)))


for i in 1:365
    hours = dayslicer(timeperiod)[i]
    println("Shifting load for hours " * string(hours) * "...")

    av_el_price = av_el_price_dict[i] 
    av_rd_price_dict = Dict((n => Dict()) for n in N)
    LoadShift_df = DataFrame(hour = Int64[], node=Int64[], OriginalLoad=Float64[], NewLoad=Float64[], ShiftUp=Float64[], ShiftDn=Float64[])

    for n in N
        for t in hours
            if nod_vol_dict[n][t] == 0
                av_rd_price_dict[n][t] = 0
            else
                av_rd_price_dict[n][t] = nod_cost_dict[n][t]/nod_vol_dict[n][t]
            end
        end
    end

        global    Load_shifted_opt, ShiftUp_opt, ShiftDn_opt, ShiftAgUp_opt, ShiftAgDn_opt, ShiftReUp_opt, ShiftReDn_opt, ShiftIndUp_opt, ShiftIndDn_opt = 
            DSO_min_LoadShift(hours, nodes, av_el_price, av_rd_price_dict, x_ag, x_re, x_ind, s_ag, s_re, s_ind, isworkhour)

    for t in hours
        Newload_df[!, string(t)] = ((Load_shifted_opt[t, n] for n in N) |> collect)
        for n in N
            newload = Load_shifted_opt[t, n]
            origload = nodes.load[n][t]
            if (newload - origload > 0)
                shiftup = newload - origload
                shiftdn = 0
            else
                shiftup = 0
                shiftdn = origload - newload
            end
            newload = Load_shifted_opt[t, n]
            origload = nodes.load[n][t]
            push!(LoadShift_df, [t n origload newload shiftup shiftdn])
        end
    end
    CSV.write(joinpath(outputpath, "Loadshift.csv"), LoadShift_df, append = true)    
    println()
    next!(prog)
    println()
end 

CSV.write(joinpath(outputpath, "Newload.csv"), Newload_df)

Newload_df
#=
sum(nodes.load[6][1:24])
sum((Load_shifted_opt[i, 6]) for i in 1:24)
nodes.load[6][1]
Load_shifted_opt[1, 6]
=#

