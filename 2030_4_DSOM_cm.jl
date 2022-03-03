""" This file executes the CM with the previously calculated, new load profiles"""

using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end

###################
# Terminal ########
###################

    usedayperiod = true
    startday = 1
    lastday = 365
    starthour = 1
    lasthour = 8760
	trm = 0	# reduces - or at negative value, increases - available capacity of transmission lines
	use221 = true	# whether the fictional generators with mc 221 are used

    datapath = "data/"
	newpath = "NewData/"
	inputpath = "output_2030/Output_LoadShift/Sens_quart/"
	outputpath = "output_2030/Output_DSOM/DSOM_Sens_quart_trm_0/"



   if usedayperiod == true
        timeperiod = string((startday-1)*24+1)*":"*string((lastday*24))
   else
        timeperiod = string(starthour)*":"*string(lasthour)
   end
    timeperiod_str = replace(string(timeperiod), ":" => "-")

    println()
    # println("Case >" * case * "< selected.")
    println("Calculating for time period >" * timeperiod * "<.")
    println()

###################
# Run CM ##########
###################

    println("Loading required packages: ProgressMeter, CSV, DataFrames, Dates, JuMP, Gurobi.")
    using ProgressMeter
    using CSV
    using DataFrames
    using Dates
    using JuMP
    using Gurobi
    println("Loading of packages complete.")

    include("model/src/Sesam_ptdf_Loadshift.jl") #WHEN USING PTDF
    # include("model/src/Sesam.jl") WHEN USING NO PTDF
    using .Sesam

###################
# Data ############
###################

	println("Reading and processing Data from CSV...")

    # nodes_df_30 =   # Dataframe mit ID, Lat, longitude
    nodes_df_30 = CSV.read(joinpath(newpath, "nodeCoordinates.csv"))
    nodes_df_30 = rename!(nodes_df_30, :node => :id)[:, 1:3]
#	nodes_df_30[:id] = nodes_df_30[:id] .+1

    # load_df_30 =    # Dataframe mit systemload und knotenweisen shares über alle stunden
    load_df_30 = CSV.read(joinpath(newpath, "demandHourly.csv"))
    rename!(load_df_30, 1 => :systemload)
    b = CSV.read(joinpath(newpath, "loadShare.csv"))
#   	b[:Column1] = b[:Column1] .+1
    for i in 1:nrow(b)
        load_df_30[Symbol(b[i, 1])] = repeat([b[i, 2]], 8760)
    end
	
    exchange_df_30 = CSV.read(joinpath(datapath, "exchange.csv"))

    # lines_df_30 =    # Dataframe mit id, from, to, voltage, reactance, resistance, circuits, pmax
    lines_df_30 = CSV.read(joinpath(newpath, "lines.csv"))
#	lines_df_30[:id] = lines_df_30[:id] .+1
#	lines_df_30[:from] = lines_df_30[:from] .+1
#	lines_df_30[:to] = lines_df_30[:to] .+1

    # powerplants_df_30 = # Dataframe mit unit, node, technology, fuel, capacity, efficiency, emission, varcost, rup, rdn, latitude, longitude
    temp = CSV.read(joinpath(newpath, "convGenData.csv"))
#	temp[:node] = temp[:node] .+ 1
    powerplants_df_30 = DataFrame(unit = String[], node = Int64[], technology = String[], fuel = String[], capacity = Float64[], efficiency = Float64[],
                                    emission = Float64[], varcost = Float64[], rup = Float64[], rdn = Float64[])
	use221 == true ? techs = names(temp)[2:end] : techs = names(temp)[2:end-1]			
    for i in 1:nrow(temp)
        for class in names(temp)[2:end]
            if !ismissing(temp[class][i]) && temp[class][i] > 0
                pp_id = "NODE"*string(temp[i,1])*"C"*string(class)
                cap = temp[class][i]
                emission = 0
                efficiency = 1
                varcost = 0
                rup = 0
                rdn = 0
                push!(powerplants_df_30, [pp_id temp[i, 1] string(class) string(class) cap efficiency emission varcost rup rdn])
            end
        end
    end

    # avail_powerplants_df_30
    avail_powerplants_df_30 = DataFrame() 
    for i in powerplants_df_30.unit
        avail_powerplants_df_30[Symbol(i)] = ones(8760)
    end    

    # mustrun_powerplants_df_30 = # Dataframe mit überall nullen pro generator
    mustrun_powerplants_df_30 = DataFrame()
    for i in powerplants_df_30.unit
        mustrun_powerplants_df_30[Symbol(i)] = zeros(8760)
    end

    # fuelcost_df_30
    fuelcost_df_30 = DataFrame()
    for i in names(temp)[2:end]
        fuelcost_df_30[i] = [parse(Int, string(i))]
    end
    fuelcost_df_30[:CO2] = [0]

    # renewables_df_30 = # Dataframe mit Node, WindOnshore, WindOffshore, SolarPV, RoR, Geothermal
    solar_df_30 = CSV.read(joinpath(newpath, "solarHourly.csv"))
    solarshare_df_30 = CSV.read(joinpath(newpath, "solarShareNodes.csv"))
#	solarshare_df_30[:node] = solarshare_df_30[:node] .+1
    wind_df_30 = CSV.read(joinpath(newpath, "windHourly.csv"))
    windshare_df_30 = CSV.read(joinpath(newpath, "windShareNodes.csv"))
#	windshare_df_30[:node] = windshare_df_30[:node] .+ 1

    renewables_df_30 = DataFrame(node=nodes_df_30[:id])
    renewables_df_30[:SolarPV] = [maximum(solar_df_30[:solarhourly]*solarshare_df_30[i,2]) for i in 1:nrow(solarshare_df_30)]
    renewables_df_30[:Wind] = [maximum(wind_df_30[:windhourly]*windshare_df_30[i,2]) for i in 1:nrow(windshare_df_30)]
    renewables_df_30

    # avail_solar_pv_30 = # Dataframe mit nodalen, stündlichen relativen Verfügbarkeiten

    # avail_solar_pv_df
    avail_solar_pv_df_30 = DataFrame()
    for n in 1:nrow(solarshare_df_30)
        temparray = []
        for h in 1:8760
            renewables_df_30[:SolarPV][n] > 0 ? x = (solarshare_df_30[n,2]*solar_df_30[h,1])/renewables_df_30[n,2] : x = 0
            push!(temparray, x)
        end
        avail_solar_pv_df_30[Symbol(n)] = temparray	
    end
	avail_solar_pv_df_30

    avail_wind_df_30 = DataFrame()
    for n in 1:nrow(windshare_df_30)
        temparray = []
        for h in 1:8760
            renewables_df_30[:Wind][n] > 0 ? x = (windshare_df_30[n,2]*wind_df_30[h,1])/renewables_df_30[n,3] : x = 0
            push!(temparray, x)
        end
        avail_wind_df_30[Symbol(n)] = temparray	
    end
    avail_wind_df_30

    renewables_avail_dict_30 = Dict(:SolarPV => avail_solar_pv_df_30,
        :Wind => avail_wind_df_30)

		println("Reading of data complete.")
		println()

###################
# Model setup #####
###################

		println("Preparing model technologies, nodes, and grid...")
		prog = Progress(6)

	# Create sets
		# Create set of nodes and lines
		nodes = Nodes(nodes_df_30, load_df_30, exchange_df_30, 1000, 0)
			next!(prog)
		lines = Lines(lines_df_30, 1000, trm)
			next!(prog)
		# Create set of power plants	ED_L_S_opt, ED_p
		pp = PowerPlants(powerplants_df_30, avail_powerplants_df_30, mustrun_powerplants_df_30, fuelcost_df_30)
			next!(prog)
		# Create set of renewables
		res = Renewables(renewables_df_30, renewables_avail_dict_30)
			next!(prog)
		# Create a set of storages
		storages_df = CSV.read(joinpath(datapath, "storages.csv"))
			next!(prog)
		storages_df[:storage] = zeros(32)
		storages_df[:power] = zeros(32)
		storages = Storages(storages_df)
			next!(prog)
		sum(sum((res.infeed[i]) for i in keys(res.infeed)))
	# Create Dict Zone to Gen
		N = nodes.id
		P = pp.unit
		R = res.unit
		S = storages.unit
		Node2PP = Dict((i, []) for i in nodes.id)
		for i in keys(pp.node)
			push!(Node2PP[pp.node[i]], i)
		end
		Node2Res = Dict((i, []) for i in nodes.id)
		for i in keys(res.node)
			push!(Node2Res[res.node[i]], i)
		end
		Node2Stor = Dict((i, []) for i in nodes.id)
		for i in keys(storages.node)
			push!(Node2Stor[storages.node[i]], i)
		end


	println("Preparation complete.")
	println()

	println("Loading new load profiles.")
	Newloads_df = CSV.read(joinpath(inputpath, "Newload.csv"))

	newloads_dict = Dict((Newloads_df[r,1] => Newloads_df[r, 2:end]) for r in 1:nrow(Newloads_df))
	nodes.load = newloads_dict

		session_time = Dates.now()
		session = Dates.format(session_time, "yyyy-mm-dd__HH-MM-SS")*"_trm_" * string(trm)


	# Create session folder for output
		mkdir(joinpath(outputpath, session*"__"*timeperiod_str))


	prog = Progress(length(dayslicer(timeperiod)))


	### Parameter outputs
		println("Exporting parameters.")
		println("Exporting merit order...")
		# Merit order
		output_ED_merit_order_pp = DataFrame(unit=Symbol[],fuel=String[],capacity=Float64[],mc=Float64[])

		try
			for g in pp.unit
					mc = pp.mc[g]
					fuel = pp.fuel[g]
					capacity = pp.capacity[g]
					push!(output_ED_merit_order_pp, [g fuel capacity mc])
			end

			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_merit_order_pp.csv"), output_ED_merit_order_pp)
		catch e
			println("The marginal costs for power plants are not available.")
		end

		output_ED_merit_order_res = DataFrame(unit=Symbol[],fuel=Symbol[],capacity=Float64[],mc=Float64[])

		try
			for r in res.unit
					mc = 0
					fuel = res.fuel[r]
					capacity = res.capacity[r]
					push!(output_ED_merit_order_res, [r fuel capacity mc])
			end

			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_merit_order_res.csv"), output_ED_merit_order_res)
		catch e
			println("The marginal costs for renewables are not available.")
		end

		output_ED_merit_order_storages = DataFrame(unit=Symbol[],fuel=String[],capacity=Float64[],mc=Float64[])

		try
			for s in storages.unit
					mc = 0
					fuel = storages.fuel[s]
					capacity = storages.power[s]
					push!(output_ED_merit_order_storages, [s fuel capacity mc])
			end

			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_merit_order_storages.csv"), output_ED_merit_order_storages)
		catch e
			println("The marginal costs for renewables are not available.")
		end

	println("Export complete.")
	println()



### Initial run of first dayslice
	println("Initial run. Calculating first dayslice: " * string(dayslicer(timeperiod)[1]))
	println()
	println()

	ED_mod, ED_P_opt, ED_P_R_opt, ED_P_R_opt_dist, ED_P_S_opt, ED_D_S_opt,
		ED_L_S_opt, ED_price = EconomicDispatch80(dayslicer(timeperiod)[1], nodes, pp, res, storages)

	CM_mod, CM_ΔP_up_opt, CM_ΔP_dn_opt, CM_ΔP_R_up_opt, CM_ΔP_R_dn_opt,
		CM_Θ_opt, CM_P_flow_opt, CM_price, CM_P_gen_lost_opt, CM_P_load_lost_opt, CM_P_S_up_opt = CongestionManagementPHS80(dayslicer(timeperiod)[1],
		nodes, lines, pp, res, storages, ED_price, ED_P_opt, ED_P_R_opt_dist,
		ED_P_S_opt, ED_D_S_opt, ED_L_S_opt)

	next!(prog)

		C_V_RD_N = DataFrame(time=Int64[], node=Int64[], CostPP=Float64[], CostRes=Float64[], CostStor=Float64[], CostLostLoad=Float64[], CostLostGen=Float64[], CostTotal=Float64[],
			VolumePP=Float64[], VolumeRes=Float64[], VolumeStor=Float64[], VolumeLostLoad=Float64[], VolumeLostGen=Float64[], VolumeTotal=Float64[])
			for t in dayslicer(timeperiod)[1]
				for n in N
					if isempty(Node2PP[n])
						costpp = 0 
					else
						costpp = sum((pp.mc[g][1]*CM_ΔP_up_opt[t, g] + (ED_price[t] - pp.mc[g][]) * CM_ΔP_dn_opt[t, g]) for g in Node2PP[n])
					end
					if isempty(Node2Res[n])
						costres = 0 
					else
						costres = sum(((ED_price[t]) * CM_ΔP_R_dn_opt[t, res]) for res in Node2Res[n])
					end
					if isempty(Node2Stor[n])
						coststor = 0 
					else
						coststor = sum(((ED_price[t]/storages.efficiency[stor]) * CM_P_S_up_opt[t, stor]) for stor in Node2Stor[n])
					end
						
					costlostload = 1000*CM_P_load_lost_opt[t, n]
					costlostgen = 1000*CM_P_gen_lost_opt[t, n] # ÜBERALL NULL!
					totalcost = costpp + costres + coststor + costlostload + costlostgen

					if isempty(Node2PP[n])
						volumepp = 0 
					else
						volumepp = sum((CM_ΔP_up_opt[t, g] - CM_ΔP_dn_opt[t, g]) for g in Node2PP[n])
					end
					if isempty(Node2Res[n])
						volumeres = 0
					else
						volumeres = sum(-(CM_ΔP_R_dn_opt[t, res]) for res in Node2Res[n])
					end
					if isempty(Node2Stor[n])
						volumestor = 0
					else
						volumestor = sum((CM_P_S_up_opt[t, stor]) for stor in Node2Stor[n])
					end
		
					volumelostload = CM_P_load_lost_opt[t, n]
					volumelostgen = CM_P_gen_lost_opt[t, n]
					totalvolume = volumepp + volumeres + volumestor + volumelostload + volumelostgen
		
					push!(C_V_RD_N, [t n costpp costres coststor costlostload costlostgen totalcost volumepp volumeres volumestor volumelostload volumelostgen totalvolume])
				end
			end
			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "C_V_RD_N__"*timeperiod_str*".csv"), C_V_RD_N)

	# Save ED Data 
		# Dispatchable power plants
		output_ED_P_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])

		try
			for t in dayslicer(timeperiod)[1]
				for g in pp.unit
					val = ED_P_opt[t, g]
					fuel =  pp.fuel[g]
					push!(output_ED_P_opt, [t g fuel val pp.node[g]])
				end
			end
			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_P.csv"), output_ED_P_opt)
		catch e
			println("ED_P_opt does not exist.")
		end

		# RES generation
		output_ED_P_R_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=Symbol[],output=Float64[],node=Int64[])
		try
			for t in dayslicer(timeperiod)[1]
				for r in res.unit
					val = ED_P_R_opt[t, r]
					fuel =  res.fuel[r]
					push!(output_ED_P_R_opt, [t r fuel val res.node[r]])
				end
			end
			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_P_R.csv"), output_ED_P_R_opt)
		catch e
			println("ED_P_R_opt does not exist.")
		end

		# Maximum available RES infeed from input data
		output_MaxRES = DataFrame(time=Int64[],unit=Symbol[],fuel=Symbol[],output=Float64[],node=Int64[])
		try
			for t in dayslicer(timeperiod)[1]
				for r in res.unit
					output = res.infeed[r][t]
					fuel =  res.fuel[r]
					push!(output_MaxRES, [t r fuel output res.node[r]])
				end
			end
			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "MaxRES.csv"), output_MaxRES)
		catch e
			println("res.infeed does not exist.")
		end

		# Generation from storages
		output_ED_P_S_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])
		try
			for t in dayslicer(timeperiod)[1]
				for s in storages.unit
					val = ED_P_S_opt[t, s]
					fuel =  storages.fuel[s]
					push!(output_ED_P_S_opt, [t s fuel val storages.node[s]])
				end
			end
			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_P_S.csv"), output_ED_P_S_opt)
		catch e
			println("ED_P_S_opt does not exist.")
		end

		# Electricity demand from storages
		output_ED_D_S_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])
		try
			for t in dayslicer(timeperiod)[1]
				for s in storages.unit
					val = -1*ED_D_S_opt[t, s]
					fuel =  storages.fuel[s]
					push!(output_ED_D_S_opt, [t s fuel val storages.node[s]])
				end
			end
			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_D_S.csv"), output_ED_D_S_opt)
		catch e
			println("ED_D_S_opt does not exist.")
		end

		# Storage level
		output_ED_L_S_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],storage=Float64[],node=Int64[])
		try
			for t in dayslicer(timeperiod)[1]
				for s in storages.unit
					val = ED_L_S_opt[t, s]
					fuel =  storages.fuel[s]
					push!(output_ED_L_S_opt, [t s fuel val storages.node[s]])
				end
			end
			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_L_S.csv"), output_ED_L_S_opt)
		catch e
			println("ED_L_S_opt does not exist.")
		end

		# Load
		output_load = DataFrame(time=Int64[],load=Float64[])
		try
			for t in dayslicer(timeperiod)[1]
				load = nodes.systemload[t]
				push!(output_load, [t load])
			end
			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "load.csv"), output_load)
		catch e
			println("nodes.systemload does not exist.")
		end

		# Prices
		output_ED_price = DataFrame(time=Int64[],price=Float64[])
		try
			for t in dayslicer(timeperiod)[1]
				price =  ED_price[t]
				push!(output_ED_price, [t price])
			end
			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_price.csv"), output_ED_price)
		catch e
			println("ED_price does not exist.")
		end

		# Must-run
		output_pmin = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])
		try
			for t in dayslicer(timeperiod)[1]
				for g in pp.unit
					val = pp.pmin[g][t]
					fuel =  pp.fuel[g]
					push!(output_pmin, [t g fuel val pp.node[g]])
				end
			end
			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "pmin.csv"), output_pmin)
		catch e
			println("pp.pmin does not exist.")
		end

		### Model parameters
		output_ED_info = DataFrame(timeperiod=String[],
										termination_status=String[],
										primal_status=String[],
										dual_status=String[],
										objective_value=Float64[],
										solve_time=Float64[])

		push!(output_ED_info, [string(dayslicer(timeperiod)[1]) string(termination_status(ED_mod)) string(primal_status(ED_mod)) string(dual_status(ED_mod)) objective_value(ED_mod) solve_time(ED_mod)])
		CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_info.csv"), output_ED_info)

		# Nodal exchange
		output_exchange = DataFrame(time=Int64[],node=Int64[],output=Float64[])
		try
			for t in dayslicer(timeperiod)[1]
				for j in nodes.id
					push!(output_exchange, [t j nodes.exchange[j][t]])
				end
			end
			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "exchange.csv"), output_exchange)
		catch e
			println("output_exchange does not exist.")
		end

###################
# CM ##############
###################
# Save CM Data
	# Dispatchable power plants
	output_CM_ΔP_up_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])
	try
		for t in dayslicer(timeperiod)[1]
			for g in pp.unit
				val = CM_ΔP_up_opt[t, g]
				fuel =  pp.fuel[g]
				push!(output_CM_ΔP_up_opt, [t g fuel val pp.node[g]])
			end
		end
		CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_up.csv"), output_CM_ΔP_up_opt)
	catch e
		println("CM_ΔP_up_opt does not exist.")
	end

	output_CM_ΔP_dn_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])
	try
		for t in dayslicer(timeperiod)[1]
			for g in pp.unit
				val = -1*CM_ΔP_dn_opt[t, g]
				fuel =  pp.fuel[g]
				push!(output_CM_ΔP_dn_opt, [t g fuel val pp.node[g]])
			end
		end
		CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_dn.csv"), output_CM_ΔP_dn_opt)
	catch e
		println("CM_ΔP_dn_opt does not exist.")
	end

	# RES generation
	output_CM_ΔP_R_up_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=Symbol[],output=Float64[],node=Int64[])
	try
		for t in dayslicer(timeperiod)[1]
			for r in res.unit
				val = CM_ΔP_R_up_opt[t, r]
				fuel =  res.fuel[r]
				push!(output_CM_ΔP_R_up_opt, [t r fuel val res.node[r]])
			end
		end
		CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_R_up.csv"), output_CM_ΔP_R_up_opt)
	catch e
		println("CM_ΔP_R_up_opt does not exist.")
	end

	output_CM_ΔP_R_dn_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=Symbol[],output=Float64[],node=Int64[])
	try
		for t in dayslicer(timeperiod)[1]
			for r in res.unit
				val = -1*CM_ΔP_R_dn_opt[t, r]
				fuel =  res.fuel[r]
				push!(output_CM_ΔP_R_dn_opt, [t r fuel val res.node[r]])
			end
		end
		CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_R_dn.csv"), output_CM_ΔP_R_dn_opt)
	catch e
		println("CM_ΔP_R_dn_opt does not exist.")
	end

	### Model parameters
	output_CM_info = DataFrame(timeperiod=String[],
									   termination_status=String[],
									   primal_status=String[],
									   dual_status=String[],
									   objective_value=Float64[],
									   solve_time=Float64[])

	push!(output_CM_info, [string(dayslicer(timeperiod)[1]) string(termination_status(CM_mod)) string(primal_status(CM_mod)) string(dual_status(CM_mod)) objective_value(CM_mod) solve_time(CM_mod)])
	CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_info.csv"), output_CM_info)

	### Power flow results
	output_CM_P_flow_opt = DataFrame(time=Int64[],
									 line=Int64[],
									 from=Int64[],
									 to=Int64[],
									 pmax=Float64[],
									 pflow=Float64[])

	 try
 		for t in dayslicer(timeperiod)[1]
 			for l in lines.id
				from = parse(Int, string(lines.from[l]))
				to = parse(Int, string(lines.to[l]))
				pmax = lines.pmax[l]
				pflow = CM_P_flow_opt[t, l]
 				push!(output_CM_P_flow_opt, [t l from to pmax pflow])
 			end
 		end
 		CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_flow.csv"), output_CM_P_flow_opt)
 	catch e
 		println("CM_P_flow_opt does not exist.")
 	end

	# Lost generation
	output_CM_P_gen_lost = DataFrame(time=Int64[],node=Int64[],output=Float64[])
	try
		for t in dayslicer(timeperiod)[1]
			for j in nodes.id
				output = CM_P_gen_lost_opt[t, j]
				push!(output_CM_P_gen_lost, [t, j, output])
			end
		end
		CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_gen_lost.csv"), output_CM_P_gen_lost)
	catch e
		println("CM_P_gen_lost does not exist.")
	end

	# Lost load
	output_CM_P_load_lost = DataFrame(time=Int64[],node=Int64[],output=Float64[])
	try
		for t in dayslicer(timeperiod)[1]
			for j in nodes.id
				output = CM_P_load_lost_opt[t, j]
				push!(output_CM_P_load_lost, [t, j, output])
			end
		end
		CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_load_lost.csv"), output_CM_P_load_lost)
	catch e
		println("CM_P_load_lost does not exist.")
	end

	output_CM_price = DataFrame(time=Int64[],node=Int64[],price=Float64[])
	 try
		 for t in dayslicer(timeperiod)[1]
			 for j in nodes.id
				 price =  CM_price[j, t]
				 push!(output_CM_price, [t j price])
			 end
		 end
		 CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_price.csv"), output_CM_price)
	 catch e
		 println("CM_price does not exist.")
	 end

	# PHS output increase
	output_CM_P_S_up_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])
	try
		for t in dayslicer(timeperiod)[1]
			for s in storages.unit
				val = CM_P_S_up_opt[t, s]
				fuel =  storages.fuel[s]
				push!(output_CM_P_S_up_opt, [t s fuel val storages.node[s]])
			end
		end
		CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_S_up.csv"), output_CM_P_S_up_opt)
	catch e
		println("CM_P_S_up_opt does not exist.")
	end

### Remaining periods: Initialising first hour with last hour of previous day
# 	and append data

	for hours in dayslicer(timeperiod)[2:end]
		hours_extended = [hours[1]-1; hours]

	println()
	println()
	println("Calculating economic dispatch: " * string(hours))

	global ED_mod, ED_P_opt, ED_P_R_opt, ED_P_R_opt_dist, ED_P_S_opt, ED_D_S_opt,
		ED_L_S_opt, ED_price = EconomicDispatch80(hours, nodes, pp, res, storages)

	    ###################
    # ED ##############
	    ###################

	    # Dispatchable power plants
	    output_ED_P_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])

	    try
	    	for t in hours
	    		for g in pp.unit
	    			val = ED_P_opt[t, g]
	    			fuel =  pp.fuel[g]
	    			push!(output_ED_P_opt, [t g fuel val pp.node[g]])
	    		end
	    	end
	    	CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_P.csv"), output_ED_P_opt, append = true)
	    catch e
	    	println("ED_P_opt does not exist.")
	    end

	    # RES generation
	    output_ED_P_R_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=Symbol[],output=Float64[],node=Int64[])
	    try
	    	for t in hours
	    		for r in res.unit
	    			val = ED_P_R_opt[t, r]
	    			fuel =  res.fuel[r]
	    			push!(output_ED_P_R_opt, [t r fuel val res.node[r]])
	    		end
	    	end
	    	CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_P_R.csv"), output_ED_P_R_opt, append = true)
	    catch e
	    	println("ED_P_R_opt does not exist.")
	    end

	    # Maximum available RES infeed from input data
	    output_MaxRES = DataFrame(time=Int64[],unit=Symbol[],fuel=Symbol[],output=Float64[],node=Int64[])
	    try
	    	for t in hours
	    		for r in res.unit
	    			output = res.infeed[r][t]
	    			fuel =  res.fuel[r]
	    			push!(output_MaxRES, [t r fuel output res.node[r]])
	    		end
	    	end
	    	CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "MaxRES.csv"), output_MaxRES, append = true)
	    catch e
	    	println("res.infeed does not exist.")
	    end

	    # Generation from storages
	    output_ED_P_S_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])
	    try
	    	for t in hours
	    		for s in storages.unit
	    			val = ED_P_S_opt[t, s]
	    			fuel =  storages.fuel[s]
	    			push!(output_ED_P_S_opt, [t s fuel val storages.node[s]])
	    		end
	    	end
	    	CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_P_S.csv"), output_ED_P_S_opt, append = true)
	    catch e
	    	println("ED_P_S_opt does not exist.")
	    end

	    # Electricity demand from storages
	    output_ED_D_S_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])
	    try
	    	for t in hours
	    		for s in storages.unit
	    			val = -1*ED_D_S_opt[t, s]
	    			fuel =  storages.fuel[s]
	    			push!(output_ED_D_S_opt, [t s fuel val storages.node[s]])
	    		end
	    	end
	    	CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_D_S.csv"), output_ED_D_S_opt, append = true)
	    catch e
	    	println("ED_D_S_opt does not exist.")
	    end

	    # Storage level
	    output_ED_L_S_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],storage=Float64[],node=Int64[])
	    try
	    	for t in hours
	    		for s in storages.unit
	    			val = ED_L_S_opt[t, s]
	    			fuel =  storages.fuel[s]
	    			push!(output_ED_L_S_opt, [t s fuel val storages.node[s]])
	    		end
	    	end
	    	CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_L_S.csv"), output_ED_L_S_opt, append = true)
	    catch e
	    	println("ED_L_S_opt does not exist.")
	    end

	    # Load
	    output_load = DataFrame(time=Int64[],load=Float64[])
	    try
	    	for t in hours
	    		load = nodes.systemload[t]
	    		push!(output_load, [t load])
	    	end
	    	CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "load.csv"), output_load, append = true)
	    catch e
	    	println("nodes.systemload does not exist.")
	    end

	    # Prices
	    output_ED_price = DataFrame(time=Int64[],price=Float64[])
	    try
	    	for t in hours
	    		price =  ED_price[t]
	    		push!(output_ED_price, [t price])
	    	end
	    	CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_price.csv"), output_ED_price, append = true)
	    catch e
	    	println("ED_price does not exist.")
	    end

	    # Must-run
	    output_pmin = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])
	    try
	    	for t in hours
	    		for g in pp.unit
	    			val = pp.pmin[g][t]
	    			fuel =  pp.fuel[g]
	    			push!(output_pmin, [t g fuel val pp.node[g]])
	    		end
	    	end
	    	CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "pmin.csv"), output_pmin, append = true)
	    catch e
	    	println("pp.pmin does not exist.")
	    end

		### Model parameters
		output_ED_info = DataFrame(timeperiod=String[],
										   termination_status=String[],
										   primal_status=String[],
										   dual_status=String[],
										   objective_value=Float64[],
										   solve_time=Float64[])

		push!(output_ED_info, [string(hours) string(termination_status(ED_mod)) string(primal_status(ED_mod)) string(dual_status(ED_mod)) objective_value(ED_mod) solve_time(ED_mod)])
		CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "ED_info.csv"), output_ED_info, append = true)


		###################
# CM ##############
		###################

		println()
		println()
		println("Calculating congestion management: " * string(hours))

			global CM_mod, CM_ΔP_up_opt, CM_ΔP_dn_opt, CM_ΔP_R_up_opt, CM_ΔP_R_dn_opt,
				CM_Θ_opt, CM_P_flow_opt, CM_price, CM_P_gen_lost_opt, CM_P_load_lost_opt, CM_P_S_up_opt = CongestionManagementPHS80(hours, nodes, lines, pp, res, storages, ED_price,
				ED_P_opt,  ED_P_R_opt_dist, ED_P_S_opt, ED_D_S_opt, ED_L_S_opt)


		C_V_RD_N = DataFrame(time=Int64[], node=Int64[], CostPP=Float64[], CostRes=Float64[], CostStor=Float64[], CostLostLoad=Float64[], CostLostGen=Float64[], CostTotal=Float64[],
			VolumePP=Float64[], VolumeRes=Float64[], VolumeStor=Float64[], VolumeLostLoad=Float64[], VolumeLostGen=Float64[], VolumeTotal=Float64[])
	    for t in hours
			for n in N
				if isempty(Node2PP[n])
					costpp = 0 
				else
					costpp = sum((pp.mc[g][1]*CM_ΔP_up_opt[t, g] + (ED_price[t] - pp.mc[g][]) * CM_ΔP_dn_opt[t, g]) for g in Node2PP[n])
				end
				if isempty(Node2Res[n])
					costres = 0 
				else
					costres = sum(((ED_price[t]) * CM_ΔP_R_dn_opt[t, res]) for res in Node2Res[n])
				end
				if isempty(Node2Stor[n])
					coststor = 0 
				else
					coststor = sum(((ED_price[t]/storages.efficiency[stor]) * CM_P_S_up_opt[t, stor]) for stor in Node2Stor[n])
				end
					
				costlostload = 1000*CM_P_load_lost_opt[t, n]
				costlostgen = 1000*CM_P_gen_lost_opt[t, n] # ÜBERALL NULL!
				totalcost = costpp + costres + coststor + costlostload + costlostgen

				if isempty(Node2PP[n])
					volumepp = 0 
				else
					volumepp = sum((CM_ΔP_up_opt[t, g] - CM_ΔP_dn_opt[t, g]) for g in Node2PP[n])
	    #            VolumePP = sum((pp.mc[g][1]*CM_ΔP_up_opt[t, g] + (ED_price[t] - pp.mc[g][]) * CM_ΔP_dn_opt[t, g]) for g in Node2PP[n])
				end
				if isempty(Node2Res[n])
					volumeres = 0
				else
					volumeres = sum(-(CM_ΔP_R_dn_opt[t, res]) for res in Node2Res[n])
				end
				if isempty(Node2Stor[n])
					volumestor = 0
				else
					volumestor = sum((CM_P_S_up_opt[t, stor]) for stor in Node2Stor[n])
				end
	
				volumelostload = CM_P_load_lost_opt[t, n]
				volumelostgen = CM_P_gen_lost_opt[t, n]
				totalvolume = volumepp + volumeres + volumestor + volumelostload + volumelostgen
	
				push!(C_V_RD_N, [t n costpp costres coststor costlostload costlostgen totalcost volumepp volumeres volumestor volumelostload volumelostgen totalvolume])
			end
		end
		CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "C_V_RD_N__"*timeperiod_str*".csv"), C_V_RD_N, append = true)

			# Dispatchable power plants
			output_CM_ΔP_up_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])
			try
				for t in hours
					for g in pp.unit
						val = CM_ΔP_up_opt[t, g]
						fuel =  pp.fuel[g]
						push!(output_CM_ΔP_up_opt, [t g fuel val pp.node[g]])
					end
				end
				CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_up.csv"), output_CM_ΔP_up_opt, append = true)
			catch e
				println("CM_ΔP_up_opt does not exist.")
			end

			output_CM_ΔP_dn_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])
			try
				for t in hours
					for g in pp.unit
						val = -1*CM_ΔP_dn_opt[t, g]
						fuel =  pp.fuel[g]
						push!(output_CM_ΔP_dn_opt, [t g fuel val pp.node[g]])
					end
				end
				CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_dn.csv"), output_CM_ΔP_dn_opt, append = true)
			catch e
				println("CM_ΔP_dn_opt does not exist.")
			end

			# RES generation
			output_CM_ΔP_R_up_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=Symbol[],output=Float64[],node=Int64[])
			try
				for t in hours
					for r in res.unit
						val = CM_ΔP_R_up_opt[t, r]
						fuel =  res.fuel[r]
						push!(output_CM_ΔP_R_up_opt, [t r fuel val res.node[r]])
					end
				end
				CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_R_up.csv"), output_CM_ΔP_R_up_opt, append = true)
			catch e
				println("CM_ΔP_R_up_opt does not exist.")
			end

			output_CM_ΔP_R_dn_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=Symbol[],output=Float64[],node=Int64[])
			try
				for t in hours
					for r in res.unit
						val = -1*CM_ΔP_R_dn_opt[t, r]
						fuel =  res.fuel[r]
						push!(output_CM_ΔP_R_dn_opt, [t r fuel val res.node[r]])
					end
				end
				CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_R_dn.csv"), output_CM_ΔP_R_dn_opt, append = true)
			catch e
				println("CM_ΔP_R_dn_opt does not exist.")
			end

			### Model parameters
			output_CM_info = DataFrame(timeperiod=String[],
											   termination_status=String[],
											   primal_status=String[],
											   dual_status=String[],
											   objective_value=Float64[],
											   solve_time=Float64[])

			push!(output_CM_info, [string(hours) string(termination_status(CM_mod)) string(primal_status(CM_mod)) string(dual_status(CM_mod)) objective_value(CM_mod) solve_time(CM_mod)])
			CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_info.csv"), output_CM_info, append = true)

			### Power flow results
			output_CM_P_flow_opt = DataFrame(time=Int64[],
											 line=Int64[],
											 from=Int64[],
											 to=Int64[],
											 pmax=Float64[],
											 pflow=Float64[])

			try
			   for t in hours
				   for l in lines.id
					   from = parse(Int, string(lines.from[l]))
					   to = parse(Int, string(lines.to[l]))
					   pmax = lines.pmax[l]
					   pflow = CM_P_flow_opt[t, l]
					   push!(output_CM_P_flow_opt, [t l from to pmax pflow])
				   end
			   end
			   CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_flow.csv"), output_CM_P_flow_opt, append = true)
		   catch e
			   println("CM_P_flow_opt does not exist.")
		   end

		   # Lost generation
		   output_CM_P_gen_lost = DataFrame(time=Int64[],node=Int64[],output=Float64[])
		   try
			   for t in hours
				   for j in nodes.id
					   output = CM_P_gen_lost_opt[t, j]
					   push!(output_CM_P_gen_lost, [t, j, output])
				   end
			   end
			   CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_gen_lost.csv"), output_CM_P_gen_lost, append = true)
		   catch e
			   println("CM_P_gen_lost does not exist.")
		   end

		   # Lost load
		   output_CM_P_load_lost = DataFrame(time=Int64[],node=Int64[],output=Float64[])
		   try
			   for t in hours
				   for j in nodes.id
					   output = CM_P_load_lost_opt[t, j]
					   push!(output_CM_P_load_lost, [t, j, output])
				   end
			   end
			   CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_load_lost.csv"), output_CM_P_load_lost, append = true)
		   catch e
			   println("CM_P_load_lost does not exist.")
		   end

		   # Noral Price
		   output_CM_price = DataFrame(time=Int64[],node=Int64[],price=Float64[])
	   	    try
	   	    	for t in hours
					for j in nodes.id
		   	    		price =  CM_price[j, t]
		   	    		push!(output_CM_price, [t j price])
					end
	   	    	end
	   	    	CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_price.csv"), output_CM_price, append = true)
	   	    catch e
	   	    	println("CM_price does not exist.")
	   	    end

			output_CM_P_S_up_opt = DataFrame(time=Int64[],unit=Symbol[],fuel=String[],output=Float64[],node=Int64[])
			try
				for t in hours
					for s in storages.unit
						val = CM_P_S_up_opt[t, s]
						fuel =  storages.fuel[s]
						push!(output_CM_P_S_up_opt, [t s fuel val storages.node[s]])
					end
				end
				CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "CM_P_S_up.csv"), output_CM_P_S_up_opt, append = true)
			catch e
				println("CM_P_S_up_opt does not exist.")
			end

			# Nodal exchange
			output_exchange = DataFrame(time=Int64[],node=Int64[],output=Float64[])
			try
				for t in hours
					for j in nodes.id
						push!(output_exchange, [t j nodes.exchange[j][t]])
					end
				end
				CSV.write(joinpath(outputpath, session*"__"*timeperiod_str, "exchange.csv"), output_exchange, append = true)
			catch e
				println("output_exchange does not exist.")
			end

	next!(prog)

end

println()
println("Calculation complete.")
