"""File for processing the data as CSV output from executable .jl files"""

"""DO NOT run this file as a whole"""
"""This file is meant for block-by-block execution in REPL"""
"""Source files may need to be manually transferred to folders or source paths changed""" 

using Pkg
Pkg.activate("Process")

using CSV
using DataFrames
using ProgressMeter

EDdatapath = "EDdata_2030" # CAN BE CHANGED TO BAUPATH
processpath = "Processed_output_2030"

BAUpath = "output_2015/Output_BAU/BAU_1-8760/"
DSOMpath = "output_2015/Output_DSOM/DSOM_1-8760/"
TSOMpath = "output_2015/Output_TSOM/TSOM_1-8760/"
LoadShiftpath = "output_2015/Output_LoadShift/"
ED = "ED"
BAU = "BAU"
TSOM = "TSOM"
DSOM = "DSOM"
Global = "Global"

# dailyfunctions - calculate corresponding hours for each day
    function dayhours(day::Int64)
        starthour = (day-1)*24+1    
        endhour = (day)*24
        return (starthour:endhour)
    end
    function starthour(day::Int64)
        starthour = (day-1)*24+1
        return starthour
    end
    function endhour(day::Int64)
        endhour = (day)*24
        return endhour
    end
# dailyfunctions END
# calculates day by day the cost over next 14 days and saves the highest one. returns cost and startday
    function highestcost2weeks(daily_df::DataFrame) # 
        highcost = 0
        index = Int64
        for i in 1:(nrow(daily_df)-14) 
            currentcost = sum(daily_df[(i:i+13), :Global_Cost])
            if currentcost > highcost
                highcost = currentcost
                index = i
            end   
        end
        return [highcost,index]
    end

### Read Nodes and info files
        numberofnodes = 485


### Read ED generations, calculate generation by tech after ED over the whole year and save

        ED_P_df = CSV.read(joinpath(EDdatapath, "ED_P.csv"), DataFrame)
        ED_P_R_df = CSV.read(joinpath(EDdatapath, "ED_P_R.csv"), DataFrame)
        ED_P_S_df = CSV.read(joinpath(EDdatapath, "ED_P_S.csv"), DataFrame)
        ED_all_df = vcat(ED_P_df, ED_P_R_df, ED_P_S_df)
        pp = names(unstack(ED_P_df, :fuel, :output))[4:end]
        res = names(unstack(ED_P_R_df, :fuel, :output))[4:end]
        stor = names(unstack(ED_P_S_df, :fuel, :output))[4:end]
        g = vcat(pp,res,stor)
        ED_yearlong_tech_dict = Dict()
        for p in pp
            ED_yearlong_tech_dict[p] = sum(subset(ED_P_df, :fuel => x -> x .== p).output)/(1*10^6)
        end
        for r in res
            ED_yearlong_tech_dict[r] = sum(subset(ED_P_R_df, :fuel => x -> x .== r).output)/(1*10^6)
        end
        for s in stor
            ED_yearlong_tech_dict[s] = sum(subset(ED_P_S_df, :fuel => x -> x .== s).output)/(1*10^6)
        end
        ED_yearlong_tech_df = DataFrame(tech = g, TotalGen = (ED_yearlong_tech_dict[p] for p in g) |> collect)

    CSV.write(joinpath(processpath, ED, "ED_yearlong_Generation_by_Tech_in_TWh.csv"), ED_yearlong_tech_df)     # Save

### Read CM changes per Case, calculate generation by tech after CM over the whole year
        ED_P_df = CSV.read(joinpath(EDdatapath, "ED_P.csv"), DataFrame)
        ED_P_R_df = CSV.read(joinpath(EDdatapath, "ED_P_R.csv"), DataFrame)
        ED_P_S_df = CSV.read(joinpath(EDdatapath, "ED_P_S.csv"), DataFrame)
        ED_all_df = vcat(ED_P_df, ED_P_R_df, ED_P_S_df)

        CM_P_up_BAU_df = CSV.read(joinpath(BAUpath, "CM_P_up.csv"), DataFrame)
        CM_P_dn_BAU_df = CSV.read(joinpath(BAUpath, "CM_P_dn.csv"), DataFrame)
        CM_P_R_up_BAU_df = CSV.read(joinpath(BAUpath, "CM_P_R_up.csv"), DataFrame)
        CM_P_R_dn_BAU_df = CSV.read(joinpath(BAUpath, "CM_P_R_dn.csv"), DataFrame)
        CM_P_S_up_BAU_df = CSV.read(joinpath(BAUpath, "CM_P_S_up.csv"), DataFrame)

        pp = names(unstack(CM_P_up_BAU_df, :fuel, :output))[4:end]
        res = names(unstack(CM_P_R_up_BAU_df, :fuel, :output))[4:end]
        stor = names(unstack(CM_P_S_up_BAU_df, :fuel, :output))[4:end]
        g = vcat(pp,res,stor)

        Gen_BAU_P_df = copy(ED_P_df)
        Gen_BAU_P_df.output = ED_P_df.output + CM_P_up_BAU_df.output + CM_P_dn_BAU_df.output
        Gen_BAU_P_R_df = copy(ED_P_R_df)
        Gen_BAU_P_R_df.output = ED_P_R_df.output + CM_P_R_up_BAU_df.output + CM_P_R_dn_BAU_df.output
        Gen_BAU_P_S_df = copy(ED_P_S_df)
        Gen_BAU_P_S_df.output = ED_P_S_df.output + CM_P_S_up_BAU_df.output

        sum(Gen_BAU_P_R_df.output)
        sum(ED_P_R_df.output)
        sum(Gen_BAU_P_df.output)
        sum(ED_P_df.output)

        Gen_BAU_Total_Dict = Dict()
        for p in pp
            Gen_BAU_Total_Dict[p] = sum(subset(Gen_BAU_P_df, :fuel => x -> x .== p).output)/(1*10^6)
        end
        for r in res
            Gen_BAU_Total_Dict[r] = sum(subset(Gen_BAU_P_R_df, :fuel => x -> x .== r).output)/(1*10^6)
        end
        for s in stor
            Gen_BAU_Total_Dict[s] = sum(subset(Gen_BAU_P_S_df, :fuel => x -> x .== s).output)/(1*10^6)
        end
        BAU_Gen_yearlong_tech_df = DataFrame(tech = g, TotalGen = (Gen_BAU_Total_Dict[p] for p in g) |> collect)

        CM_P_up_TSOM_df = CSV.read(joinpath(TSOMpath, "CM_P_up.csv"), DataFrame)
        CM_P_dn_TSOM_df = CSV.read(joinpath(TSOMpath, "CM_P_dn.csv"), DataFrame)
        CM_P_R_up_TSOM_df = CSV.read(joinpath(TSOMpath, "CM_P_R_up.csv"), DataFrame)
        CM_P_R_dn_TSOM_df = CSV.read(joinpath(TSOMpath, "CM_P_R_dn.csv"), DataFrame)
        CM_P_S_up_TSOM_df = CSV.read(joinpath(TSOMpath, "CM_P_S_up.csv"), DataFrame)

        Gen_TSOM_P_df = copy(ED_P_df)
        Gen_TSOM_P_df.output = ED_P_df.output + CM_P_up_TSOM_df.output + CM_P_dn_TSOM_df.output
        Gen_TSOM_P_R_df = copy(ED_P_R_df)
        Gen_TSOM_P_R_df.output = ED_P_R_df.output + CM_P_R_up_TSOM_df.output + CM_P_R_dn_TSOM_df.output
        Gen_TSOM_P_S_df = copy(ED_P_S_df)
        Gen_TSOM_P_S_df.output = ED_P_S_df.output + CM_P_S_up_TSOM_df.output

        Gen_TSOM_Total_Dict = Dict()
        for p in pp
            Gen_TSOM_Total_Dict[p] = sum(subset(Gen_TSOM_P_df, :fuel => x -> x .== p).output)/(1*10^6)
        end
        for r in res
            Gen_TSOM_Total_Dict[r] = sum(subset(Gen_TSOM_P_R_df, :fuel => x -> x .== r).output)/(1*10^6)
        end
        for s in stor
            Gen_TSOM_Total_Dict[s] = sum(subset(Gen_TSOM_P_S_df, :fuel => x -> x .== s).output)/(1*10^6)
        end
        TSOM_Gen_yearlong_tech_df = DataFrame(tech = g, TotalGen = (Gen_TSOM_Total_Dict[p] for p in g) |> collect)

        CM_P_up_DSOM_df = CSV.read(joinpath(DSOMpath, "CM_P_up.csv"), DataFrame)
        CM_P_dn_DSOM_df = CSV.read(joinpath(DSOMpath, "CM_P_dn.csv"), DataFrame)
        CM_P_R_up_DSOM_df = CSV.read(joinpath(DSOMpath, "CM_P_R_up.csv"), DataFrame)
        CM_P_R_dn_DSOM_df = CSV.read(joinpath(DSOMpath, "CM_P_R_dn.csv"), DataFrame)
        CM_P_S_up_DSOM_df = CSV.read(joinpath(DSOMpath, "CM_P_S_up.csv"), DataFrame)

        Gen_DSOM_P_df = copy(ED_P_df)
        Gen_DSOM_P_df.output = ED_P_df.output + CM_P_up_DSOM_df.output + CM_P_dn_DSOM_df.output
        Gen_DSOM_P_R_df = copy(ED_P_R_df)
        Gen_DSOM_P_R_df.output = ED_P_R_df.output + CM_P_R_up_DSOM_df.output + CM_P_R_dn_DSOM_df.output
        Gen_DSOM_P_S_df = copy(ED_P_S_df)
        Gen_DSOM_P_S_df.output = ED_P_S_df.output + CM_P_S_up_DSOM_df.output

        Gen_DSOM_Total_Dict = Dict()
        for p in pp
            Gen_DSOM_Total_Dict[p] = sum(subset(Gen_DSOM_P_df, :fuel => x -> x .== p).output)/(1*10^6)
        end
        for r in res
            Gen_DSOM_Total_Dict[r] = sum(subset(Gen_DSOM_P_R_df, :fuel => x -> x .== r).output)/(1*10^6)
        end
        for s in stor
            Gen_DSOM_Total_Dict[s] = sum(subset(Gen_DSOM_P_S_df, :fuel => x -> x .== s).output)/(1*10^6)
        end

        DSOM_Gen_yearlong_tech_df = DataFrame(tech = g, TotalGen = (Gen_DSOM_Total_Dict[p] for p in g) |> collect)
    CSV.write(joinpath(processpath, "BAU_Gen_yearlong_tech.csv"), BAU_Gen_yearlong_tech_df)
    CSV.write(joinpath(processpath, "TSOM_Gen_yearlong_tech.csv"), TSOM_Gen_yearlong_tech_df)
    CSV.write(joinpath(processpath, "DSOM_Gen_yearlong_tech.csv"), DSOM_Gen_yearlong_tech_df)


    sum(ED_yearlong_tech_df.TotalGen)
    sum(BAU_Gen_yearlong_tech_df.TotalGen)
    sum(TSOM_Gen_yearlong_tech_df.TotalGen)
    sum(DSOM_Gen_yearlong_tech_df.TotalGen)

### Calculate Daily and hourly nodal Cost and Volume of RD (BAU)
    
            working_df = CSV.read(joinpath(BAUpath, "C_V_RD_N__1-8760.csv"), DataFrame)
            nodal_hourlycost_df = unstack(working_df[!,[:node, :time, :CostTotal]], :time, :CostTotal)
            nodal_hourlyvolume_df = unstack(working_df[!,[:node, :time, :VolumeTotal]], :time, :VolumeTotal)
            nodal_dailycost_df = DataFrame(node=nodal_hourlycost_df.node)
            for i in 1:365
                hours=(i-1)*24+2:i*24+1 
                nodal_dailycost_df[!, string(i)] = sum(eachcol(nodal_hourlycost_df[!,hours]))
            end
            
            nodal_dailyvolume_df = DataFrame(node=nodal_hourlyvolume_df.node)
            for i in 1:365
                hours=(i-1)*24+2:i*24+1
                nodal_dailyvolume_df[!, string(i)] = sum(eachcol(nodal_hourlyvolume_df[!,hours]))
            end

            nodal_twoweek_cost = DataFrame(node=Int64[], RDcost = Float64[])
            for i in 1:nrow(nodal_dailycost_df)
                n = nodal_dailycost_df[i,1]
                rdcost = sum(nodal_dailycost_df[i,322:335])
                push!(nodal_twoweek_cost, [n rdcost])
            end
            nodal_twoweek_RDcost = CSV.read("NewData/nodeCoordinates.csv", DataFrame)
            nodal_twoweek_RDcost.RDcost = nodal_twoweek_cost.RDcost
            CSV.write(joinpath(processpath, BAU, "nodal_twoweek_RDcost.csv"), nodal_twoweek_RDcost)

            
            nodal_annual_cost = DataFrame(node=Int64[], RDcost = Float64[])
            for i in 1:nrow(nodal_dailycost_df)
                n = nodal_dailycost_df[i,1]
                rdcost = sum(nodal_dailycost_df[i,2:end])
                push!(nodal_annual_cost, [n rdcost])
            end
            nodal_annual_RDcost = CSV.read("NewData/nodeCoordinates.csv", DataFrame)
            nodal_annual_RDcost.RDcost = nodal_annual_cost.RDcost
            CSV.write(joinpath(processpath, BAU, "nodal_annual_RDcost.csv"), nodal_annual_RDcost)


            

    nodal_dailyvolume_df
    CSV.write(joinpath(processpath, BAU, "hourly_nodal_RD_volume.csv"), nodal_hourlyvolume_df)
    CSV.write(joinpath(processpath, BAU, "hourly_nodal_RD_cost.csv"), nodal_hourlycost_df)
    CSV.write(joinpath(processpath, BAU, "daily_nodal_RD_volume.csv"), nodal_dailyvolume_df)
    CSV.write(joinpath(processpath, BAU, "daily_nodal_RD_cost.csv"), nodal_dailycost_df)

### Load total cost of RD for scenarios TOTAL COST OVER YEAR

    BAU_info_df = CSV.read(joinpath(BAUpath, "CM_info.csv"), DataFrame)
    TSOM_info_df = CSV.read(joinpath(TSOMpath, "CM_info.csv"), DataFrame)
    DSOM_info_df = CSV.read(joinpath(DSOMpath, "CM_info.csv"), DataFrame)

    Total_daily_RD_cost_bycase_df = DataFrame(day=1:365, BAU =BAU_info_df.objective_value, TSOM = TSOM_info_df.objective_value, DSOM = DSOM_info_df.objective_value)
    Total_daily_RD_cost_bycase_df = DataFrame(day=1:365, BAU =BAU_info_df.objective_value, DSOM = DSOM_info_df.objective_value)

    Total_annual_RD_cost_bycase_df = DataFrame(case=["BAU", "DSOM", "TSOM"], cost = [sum(Total_daily_RD_cost_bycase_df.BAU), sum(Total_daily_RD_cost_bycase_df.DSOM), sum(Total_daily_RD_cost_bycase_df.TSOM),])
    CSV.write(joinpath(processpath, Global, "Total_annual_RD_cost_bycase.csv"), Total_annual_RD_cost_bycase_df)

    (sum(Total_daily_RD_cost_bycase_df.BAU)-sum(Total_daily_RD_cost_bycase_df.TSOM))/sum(Total_daily_RD_cost_bycase_df.BAU)*100
    (sum(Total_daily_RD_cost_bycase_df.BAU)-sum(Total_daily_RD_cost_bycase_df.DSOM))/sum(Total_daily_RD_cost_bycase_df.BAU)*100
    sum(Total_daily_RD_cost_bycase_df.BAU)/(10^6)
    sum(Total_daily_RD_cost_bycase_df.DSOM)/(10^6)
    sum(Total_daily_RD_cost_bycase_df.TSOM)/(10^6)
    sum(Total_daily_RD_cost_bycase_df[321:334, :BAU])
    (sum(Total_daily_RD_cost_bycase_df[321:334, :BAU])-sum(Total_daily_RD_cost_bycase_df[321:334, :TSOM]))/sum(Total_daily_RD_cost_bycase_df[321:334, :BAU])
    (sum(Total_daily_RD_cost_bycase_df[321:334, :BAU])-sum(Total_daily_RD_cost_bycase_df[321:334, :DSOM]))/sum(Total_daily_RD_cost_bycase_df[321:334, :BAU])

### Output new Load profiles, Load Duration Curves

    LoadShift_TSOM_twoweeks_df = subset(CSV.read(joinpath(TSOMpath, "CM_LoadShift.csv"), DataFrame), :time => x -> 7681 .<= x .<= 8016)
    LoadShift_TSOM_twoweeks_df.ShiftUp = zeros(nrow(LoadShift_TSOM_twoweeks_df))
    LoadShift_TSOM_twoweeks_df.ShiftDn = zeros(nrow(LoadShift_TSOM_twoweeks_df))
    for i in 1:nrow(LoadShift_TSOM_twoweeks_df)
        if LoadShift_TSOM_twoweeks_df[i, 5] > 0
                LoadShift_TSOM_twoweeks_df[i, 6] = LoadShift_TSOM_twoweeks_df[i, 5]
        elseif LoadShift_TSOM_twoweeks_df[i, 5] < 0
                LoadShift_TSOM_twoweeks_df[i, 7] = -LoadShift_TSOM_twoweeks_df[i, 5]
        end
    end
    sum(LoadShift_TSOM_twoweeks_df.ShiftUp)
    sum(LoadShift_TSOM_twoweeks_df.ShiftDn)
    LoadShift_DSOM_twoweeks_df = subset(CSV.read(joinpath(LoadShiftpath, "Loadshift.csv"), DataFrame), :hour => x -> 7681 .<= x .<= 8016)
    sum(LoadShift_DSOM_twoweeks_df.ShiftUp)
    sum(LoadShift_DSOM_twoweeks_df.ShiftDn)
    LoadShift_DSOM_df = CSV.read(joinpath(LoadShiftpath, "LoadShift.csv"), DataFrame)
    LoadShift_TSOM_df = CSV.read(joinpath(TSOMpath, "CM_Loadshift.csv"), DataFrame)

    Loads_yearlong_BYCASE_df = DataFrame(time = 1:8760, 
        OriginalLoad = (sum(LoadShift_DSOM_df[(h-1)*numberofnodes+1:h*numberofnodes, :OriginalLoad]) for h in 1:8760) |> collect,
        LoadDSOM = (sum(LoadShift_DSOM_df[(h-1)*numberofnodes+1:h*numberofnodes, :NewLoad]) for h in 1:8760) |> collect
        )

    
    LoadDuration_yearlong_BYCASE_df = DataFrame(id=1:8760,
        BAU = sort((sum(LoadShift_DSOM_df[(h-1)*numberofnodes+1:h*numberofnodes, :OriginalLoad]) for h in 1:8760) |> collect, rev = true),
        DSOM = sort((sum(LoadShift_DSOM_df[(h-1)*numberofnodes+1:h*numberofnodes, :NewLoad]) for h in 1:8760) |> collect, rev = true),
        TSOM = sort((sum(LoadShift_TSOM_df[(h-1)*numberofnodes+1:h*numberofnodes, :ShiftedLoad]) for h in 1:8760) |> collect, rev = true)
        )


    CSV.write(joinpath(processpath, Global, "LoadDurationCurve_yearlong_BYCASE.csv"), LoadDuration_yearlong_BYCASE_df)

    Loads_twoweeks_BYCASE_df = DataFrame(time = 7681:8016, 
        OriginalLoad = (sum(LoadShift_DSOM_twoweeks_df[(h-1)*numberofnodes+1:h*numberofnodes, :OriginalLoad]) for h in 1:336) |> collect, 
        LoadTSOM = (sum(LoadShift_TSOM_twoweeks_df[(h-1)*numberofnodes+1:h*numberofnodes, :ShiftedLoad]) for h in 1:336) |> collect,
        LoadDSOM = (sum(LoadShift_DSOM_twoweeks_df[(h-1)*numberofnodes+1:h*numberofnodes, :NewLoad]) for h in 1:336) |> collect
        )

    CSV.write(joinpath(processpath, Global, "Loads_twoweeks_BYCASE.csv"), Loads_twoweeks_BYCASE_df)


### Calculate daily RD for days 321-334 by fuel, GESAMTVOLUMEN RD, RD Profiles

    CM_P_up_BAU_df = CSV.read(joinpath(BAUpath, "CM_P_up.csv"), DataFrame)
    CM_P_dn_BAU_df = CSV.read(joinpath(BAUpath, "CM_P_dn.csv"), DataFrame)
    CM_P_R_up_BAU_df = CSV.read(joinpath(BAUpath, "CM_P_R_up.csv"), DataFrame)
    CM_P_R_dn_BAU_df = CSV.read(joinpath(BAUpath, "CM_P_R_dn.csv"), DataFrame)
    # CM_P_S_up_BAU_df = CSV.read(joinpath(BAUpath, "CM_P_S_up.csv"), DataFrame)
    CM_LL_BAU_df = CSV.read(joinpath(BAUpath, "CM_P_load_lost.csv"), DataFrame)
    CM_all_up_BAU_df = vcat(CM_P_up_BAU_df, CM_P_R_up_BAU_df
    #, CM_P_S_up_BAU_df
    )
    CM_all_dn_BAU_df = vcat(CM_P_dn_BAU_df, CM_P_R_dn_BAU_df)
    
    CM_twoweek_RDup_BAU_df = subset(CM_all_up_BAU_df, :time => x -> 7681 .<= x .<= 8016)
    CM_twoweek_RDdn_BAU_df = subset(CM_all_dn_BAU_df, :time => x -> 7681 .<= x .<= 8016)
    CM_twoweek_LL_BAU_df = subset(CM_LL_BAU_df, :time => x -> 7681 .<= x .<= 8016)

    sum(subset(CM_P_R_dn_BAU_df, :time => x -> 7681 .<= x .<= 8016).output)

    TwoWeekVol_RD_BAU = DataFrame(SummeUp = sum(CM_twoweek_RDup_BAU_df.output), SummeLL = sum(CM_twoweek_LL_BAU_df.output), SummeDn = sum(CM_twoweek_RDdn_BAU_df.output))
    TotalVol_RD_BAU = DataFrame(SummeUp = sum(CM_all_up_BAU_df.output), SummeLL = sum(CM_LL_BAU_df.output), SummeDn = sum(CM_all_dn_BAU_df.output))

    CSV.write(joinpath(processpath, BAU, "BAU_RD_Gesamtvolumen.csv"), TotalVol_RD_BAU)
    sum(CM_all_up_BAU_df.output)+sum(CM_LL_BAU_df.output)
    sum(CM_all_dn_BAU_df.output)

    g = unique(CM_twoweek_RDup_BAU_df.fuel)

    RDprofile_BAU_twoweeks_df = DataFrame(hour=Int64[], fuel=String[], RD_up=Float64[], RD_down=Float64[], RD_net=Float64[])
    begin 
            prog = Progress(14*24)
        for h in 7681:8016
            for z in g
                up = sum(subset(CM_twoweek_RDup_BAU_df, :fuel => x -> x .== z, :time => y -> y .== h).output)
                dn = sum(subset(CM_twoweek_RDdn_BAU_df, :fuel => x -> x .== z, :time => y -> y .== h).output)
                net = up + dn
            push!(RDprofile_BAU_twoweeks_df, [h, string(z), up, dn, net])
            end
            ll = sum(subset(CM_twoweek_LL_BAU_df, :time => y -> y .== h).output)
            push!(RDprofile_BAU_twoweeks_df, [h "LostLoad" ll 0 ll])
            next!(prog)
        end
    end

    RDprofile_BAU_twoweeks_df = vcat(rename!(RDprofile_BAU_twoweeks_df[!, [:hour, :fuel, :RD_up]], :RD_up => :RD), rename!(RDprofile_BAU_twoweeks_df[!, [:hour, :fuel, :RD_down]], :RD_down => :RD))

    CM_P_up_TSOM_df = CSV.read(joinpath(TSOMpath, "CM_P_up.csv"), DataFrame)
    CM_P_dn_TSOM_df = CSV.read(joinpath(TSOMpath, "CM_P_dn.csv"), DataFrame)
    CM_P_R_up_TSOM_df = CSV.read(joinpath(TSOMpath, "CM_P_R_up.csv"), DataFrame)
    CM_P_R_dn_TSOM_df = CSV.read(joinpath(TSOMpath, "CM_P_R_dn.csv"), DataFrame)
    # CM_P_S_up_TSOM_df = CSV.read(joinpath(TSOMpath, "CM_P_S_up.csv"), DataFrame)
    CM_LL_TSOM_df = CSV.read(joinpath(TSOMpath, "CM_P_load_lost.csv"), DataFrame)
    CM_LS_TSOM_df = CSV.read(joinpath(TSOMpath, "CM_Loadshift.csv"), DataFrame)
    CM_all_up_TSOM_df = vcat(CM_P_up_TSOM_df, CM_P_R_up_TSOM_df
    #, CM_P_S_up_TSOM_df
    )
    CM_all_dn_TSOM_df = vcat(CM_P_dn_TSOM_df, CM_P_R_dn_TSOM_df)

    CM_LS_TSOM_df.ShiftUp = zeros(nrow(CM_LS_TSOM_df))
    CM_LS_TSOM_df.ShiftDn = zeros(nrow(CM_LS_TSOM_df))
    for i in 1:nrow(CM_LS_TSOM_df)
        if CM_LS_TSOM_df[i, 5] > 0
                CM_LS_TSOM_df[i, 6] = CM_LS_TSOM_df[i, 5]
        elseif CM_LS_TSOM_df[i, 5] < 0
                CM_LS_TSOM_df[i, 7] = -CM_LS_TSOM_df[i, 5]
        end
    end

    sum(CM_LS_TSOM_df.ShiftDn)

    # TWO WEEK SUM OF CURTAILMENT
    sum(subset(CM_P_R_dn_TSOM_df, :time => x -> 7681 .<= x .<= 8016).output)

    CM_twoweek_RDup_TSOM_df = subset(CM_all_up_TSOM_df, :time => x -> 7681 .<= x .<= 8016)
    CM_twoweek_RDdn_TSOM_df = subset(CM_all_dn_TSOM_df, :time => x -> 7681 .<= x .<= 8016)
    CM_twoweek_LL_TSOM_df = subset(CM_LL_TSOM_df, :time => x -> 7681 .<= x .<= 8016)
    CM_twoweek_LS_TSOM_df = subset(CM_LS_TSOM_df, :time => x -> 7681 .<= x .<= 8016)


    TwoWeekVol_RD_TSOM = DataFrame(SummeUp = sum(CM_twoweek_RDup_TSOM_df.output), SummeLL = sum(CM_twoweek_LL_TSOM_df.output), SummeDn = sum(CM_twoweek_RDdn_TSOM_df.output))

    TotalVol_RD_TSOM = DataFrame(SummeUp = sum(CM_all_up_TSOM_df.output), SummeLL = sum(CM_LL_TSOM_df.output), SummeDn = sum(CM_all_dn_TSOM_df.output))

    CSV.write(joinpath(processpath, TSOM, "TSOM_RD_Gesamtvolumen.csv"), TotalVol_RD_TSOM)
    sum(CM_all_up_TSOM_df.output)+sum(CM_LL_TSOM_df.output)
    sum(CM_all_dn_TSOM_df.output)

    RDprofile_TSOM_twoweeks_df = DataFrame(hour=Int64[], fuel=String[], RD_up=Float64[], RD_down=Float64[], RD_net=Float64[])

    g = unique(CM_twoweek_RDup_TSOM_df.fuel)

    begin
        prog = Progress(14*24)
        for h in 7681:8016
            for z in g
                up = sum(subset(CM_twoweek_RDup_TSOM_df, :fuel => x -> x .== z, :time => y -> y .== h).output)
                dn = sum(subset(CM_twoweek_RDdn_TSOM_df, :fuel => x -> x .== z, :time => y -> y .== h).output)
                net = up + dn
                push!(RDprofile_TSOM_twoweeks_df, [h, string(z), up, dn, net])
            end
            ll = sum(subset(CM_twoweek_LL_TSOM_df, :time => y -> y .== h).output)
            lsup = -sum(subset(CM_twoweek_LS_TSOM_df, :time => y -> y .== h).ShiftUp)
            lsdn = sum(subset(CM_twoweek_LS_TSOM_df, :time => y -> y .== h).ShiftDn)
            push!(RDprofile_TSOM_twoweeks_df, [h "LostLoad" ll 0 ll])
            push!(RDprofile_TSOM_twoweeks_df, [h "LoadShift" lsup lsdn lsup+lsdn])
            next!(prog)
        end
    end


    RDprofile_TSOM_twoweeks_df = vcat(rename!(RDprofile_TSOM_twoweeks_df[!, [:hour, :fuel, :RD_up]], :RD_up => :RD_tot), rename!(RDprofile_TSOM_twoweeks_df[!, [:hour, :fuel, :RD_down]], :RD_down => :RD_tot))

    CM_P_up_DSOM_df = CSV.read(joinpath(DSOMpath, "CM_P_up.csv"), DataFrame)
    CM_P_dn_DSOM_df = CSV.read(joinpath(DSOMpath, "CM_P_dn.csv"), DataFrame)
    CM_P_R_up_DSOM_df = CSV.read(joinpath(DSOMpath, "CM_P_R_up.csv"), DataFrame)
    CM_P_R_dn_DSOM_df = CSV.read(joinpath(DSOMpath, "CM_P_R_dn.csv"), DataFrame)
    # CM_P_S_up_DSOM_df = CSV.read(joinpath(DSOMpath, "CM_P_S_up.csv"), DataFrame)
    CM_LL_DSOM_df = CSV.read(joinpath(DSOMpath, "CM_P_load_lost.csv"), DataFrame)
    CM_LS_DSOM_df = CSV.read(joinpath(LoadShiftpath, "Loadshift.csv"), DataFrame)
    CM_all_up_DSOM_df = vcat(CM_P_up_DSOM_df, CM_P_R_up_DSOM_df 
    #, CM_P_S_up_DSOM_df
    )
    CM_all_dn_DSOM_df = vcat(CM_P_dn_DSOM_df, CM_P_R_dn_DSOM_df)


    # TWO WEEK SUM OF CURTAILMENT
    sum(subset(CM_P_R_dn_DSOM_df, :time => x -> 7681 .<= x .<= 8016).output)    

    CM_twoweek_RDup_DSOM_df = subset(CM_all_up_DSOM_df, :time => x -> 7681 .<= x .<= 8016)
    CM_twoweek_RDdn_DSOM_df = subset(CM_all_dn_DSOM_df, :time => x -> 7681 .<= x .<= 8016)
    CM_twoweek_LL_DSOM_df = subset(CM_LL_DSOM_df, :time => x -> 7681 .<= x .<= 8016)
    CM_twoweek_LS_DSOM_df = subset(CM_LS_DSOM_df, :hour => x -> 7681 .<= x .<= 8016)

    TwoWeekVol_RD_DSOM = DataFrame(SummeUp = sum(CM_twoweek_RDup_DSOM_df.output), SummeLL = sum(CM_twoweek_LL_DSOM_df.output), SummeDn = sum(CM_twoweek_RDdn_DSOM_df.output))
    TotalVol_RD_DSOM = DataFrame(SummeUp = sum(CM_all_up_DSOM_df.output), SummeLL = sum(CM_LL_DSOM_df.output), SummeDn = sum(CM_all_dn_DSOM_df.output))
    
    CSV.write(joinpath(processpath, DSOM, "DSOM_RD_Gesamtvolumen.csv"), TotalVol_RD_DSOM)

    g = unique(CM_twoweek_RDup_DSOM_df.fuel)

    RDprofile_DSOM_twoweeks_df = DataFrame(hour=Int64[], fuel=String[], RD_up=Float64[], RD_down=Float64[], RD_net=Float64[])

    for h in 7681:8016
        for z in g
            up = sum(subset(CM_twoweek_RDup_DSOM_df, :fuel => x -> x .== z, :time => y -> y .== h).output)
            dn = sum(subset(CM_twoweek_RDdn_DSOM_df, :fuel => x -> x .== z, :time => y -> y .== h).output)
            net = up + dn
        push!(RDprofile_DSOM_twoweeks_df, [h, string(z), up, dn, net])
        end
        ll = sum(subset(CM_twoweek_LL_DSOM_df, :time => y -> y .== h).output)
        lsup = -sum(subset(CM_twoweek_LS_DSOM_df, :hour => y -> y .== h).ShiftUp)
        lsdn = sum(subset(CM_twoweek_LS_DSOM_df, :hour => y -> y .== h).ShiftDn)
        push!(RDprofile_DSOM_twoweeks_df, [h "LostLoad" ll 0 ll])
        push!(RDprofile_DSOM_twoweeks_df, [h "LoadShift" lsup lsdn lsup+lsdn])
    end

    RDprofile_DSOM_twoweeks_df = vcat(rename!(RDprofile_DSOM_twoweeks_df[!, [:hour, :fuel, :RD_up]], :RD_up => :RD_tot), rename!(RDprofile_DSOM_twoweeks_df[!, [:hour, :fuel, :RD_down]], :RD_down => :RD_tot))

    TotalVol_RD_Global_yearlong_BYCASE_df = insertcols!(vcat(TotalVol_RD_BAU, TotalVol_RD_DSOM, TotalVol_RD_TSOM), 1, :Case => ["BAU","DSOM","TSOM"])

    CSV.write(joinpath(processpath, Global, "Total_RD_Vol_yearlong_BYCASE.csv"), TotalVol_RD_Global_yearlong_BYCASE_df)    
    CSV.write(joinpath(processpath, BAU, "RDprofile_BAU_twoweeks.csv"), RDprofile_BAU_twoweeks_df)
    CSV.write(joinpath(processpath, "Sens_quart", "RDprofile_DSOM_twoweeks.csv"), RDprofile_DSOM_twoweeks_df)
    CSV.write(joinpath(processpath, "Sens_quart", "RDprofile_TSOM_twoweeks.csv"), RDprofile_TSOM_twoweeks_df)

    (sum(CM_all_dn_BAU_df.output) - sum(CM_all_dn_DSOM_df.output))/sum(CM_all_dn_BAU_df.output)
    (sum(CM_all_dn_BAU_df.output) - sum(CM_all_dn_TSOM_df.output))

    BAU_info_df = CSV.read(joinpath(BAUpath, "CM_info.csv"), DataFrame)
    TSOM_info_df = CSV.read(joinpath(TSOMpath, "CM_info.csv"), DataFrame)
    DSOM_info_df = CSV.read(joinpath(DSOMpath, "CM_info.csv"), DataFrame)

    """CALCULATE SAVE PER MWh OF LOAD SHIFTING"""

    Save_per_MWh_BYCASE_df = DataFrame(Case=["TSOM", "DSOM"], Save=
        [(sum(BAU_info_df.objective_value)-sum(TSOM_info_df.objective_value))/sum(CM_LS_TSOM_df.ShiftUp), (sum(BAU_info_df.objective_value)-sum(DSOM_info_df.objective_value))/sum(CM_LS_DSOM_df.ShiftUp)])
    Save_per_MWh_BYCASE_twoweeks_df = DataFrame(Case=["TSOM", "DSOM"], Save=
        [(sum(BAU_info_df.objective_value[321:334])-sum(TSOM_info_df.objective_value[321:334]))/sum(CM_LS_TSOM_df.ShiftUp[320*numberofnodes*24+1:334*numberofnodes*24]), 
        (sum(BAU_info_df.objective_value[321:334])-sum(DSOM_info_df.objective_value[321:334]))/sum(CM_LS_DSOM_df.ShiftUp[320*numberofnodes*24+1:334*numberofnodes*24])])
    CSV.write(joinpath(processpath, "Sens_quart", "Save_per_MWh_twoweeks_BYCASE.csv"), Save_per_MWh_BYCASE_twoweeks_df)

nod_vol_dict = Dict((nodal_dailyvolume_df[i, 1] => Dict((j => nodal_dailyvolume_df[i, j+1]) for j in 1:365)) for i in 1:nrow(nodal_dailyvolume_df))
nod_cost_dict = Dict((nodal_dailycost_df[i, 1] => Dict((j => nodal_dailycost_df[i, j+1]) for j in 1:365)) for i in 1:nrow(nodal_dailycost_df))

### Calculate Daily and hourly global Cost and Volume of RD, sum for year and twoweeks

        working_df = CSV.read(joinpath(BAUpath, "C_V_RD_N__1-8760.csv"), DataFrame)
        daily_global_BAU_df = DataFrame(day=1:365, 
            Global_Cost=((sum((working_df[s, :CostTotal]) for s in (endhour(i-1)*numberofnodes+1):endhour(i)*numberofnodes) for i in 1:365) |> collect),
            Global_Vol_Res=((sum((working_df[s, :VolumeRes]) for s in (endhour(i-1)*numberofnodes+1):endhour(i)*numberofnodes) for i in 1:365) |> collect),
            Global_Volume=((sum((working_df[s, :VolumeTotal]) for s in (endhour(i-1)*numberofnodes+1):endhour(i)*numberofnodes) for i in 1:365) |> collect),
            )
        hourly_global_BAU_df = DataFrame(hour=1:8760, 
            Total_RD_Cost=((sum((working_df[s, :CostTotal]) for s in ((i-1)*numberofnodes+1:i*numberofnodes)) for i in 1:8760) |> collect),
            Global_Vol_Res=((sum((working_df[s, :VolumeRes]) for s in ((i-1)*numberofnodes+1:i*numberofnodes)) for i in 1:8760) |> collect),
            Global_Vol_PP=((sum((working_df[s, :VolumePP]) for s in ((i-1)*numberofnodes+1:i*numberofnodes)) for i in 1:8760) |> collect),
            Global_Volume=((sum((working_df[s, :VolumeTotal]) for s in ((i-1)*numberofnodes+1:i*numberofnodes)) for i in 1:8760) |> collect),
            )
            highcostweeks=highestcost2weeks(daily_global_BAU_df)
            highcoststartday = Int64(highcostweeks[2])
            period_string = string(highcoststartday) * "_to_" * string(highcoststartday+13)
            print("the highest cost over two weeks is " * string(highcostweeks[1]) * " in the two weeks starting with day " * string(highcoststartday))
        working_df = CSV.read(joinpath(TSOMpath, "C_V_RD_N.csv"), DataFrame)
        daily_global_TSOM_df = DataFrame(day=1:365, 
            Global_Cost=((sum((working_df[s, :CostTotal]) for s in (endhour(i-1)*numberofnodes+1):endhour(i)*numberofnodes) for i in 1:365) |> collect),
            Global_Vol_Res=((sum((working_df[s, :VolumeRes]) for s in (endhour(i-1)*numberofnodes+1):endhour(i)*numberofnodes) for i in 1:365) |> collect),
            Global_Volume=((sum((working_df[s, :VolumeTotal]) for s in (endhour(i-1)*numberofnodes+1):endhour(i)*numberofnodes) for i in 1:365) |> collect),
            )
        hourly_global_TSOM_df = DataFrame(hour=1:8760, 
            Total_RD_Cost=((sum((working_df[s, :CostTotal]) for s in ((i-1)*numberofnodes+1:i*numberofnodes)) for i in 1:8760) |> collect),
            Global_Vol_Res=((sum((working_df[s, :VolumeRes]) for s in ((i-1)*numberofnodes+1:i*numberofnodes)) for i in 1:8760) |> collect),
            Global_Vol_PP=((sum((working_df[s, :VolumePP]) for s in ((i-1)*numberofnodes+1:i*numberofnodes)) for i in 1:8760) |> collect),
            Global_Volume=((sum((working_df[s, :VolumeTotal]) for s in ((i-1)*numberofnodes+1:i*numberofnodes)) for i in 1:8760) |> collect),
            )
        working_df = CSV.read(joinpath(DSOMpath, "C_V_RD_N__1-8760.csv"), DataFrame)
        daily_global_DSOM_df = DataFrame(day=1:365, 
            Global_Cost=((sum((working_df[s, :CostTotal]) for s in (endhour(i-1)*numberofnodes+1):endhour(i)*numberofnodes) for i in 1:365) |> collect),
            Global_Vol_Res=((sum((working_df[s, :VolumeRes]) for s in (endhour(i-1)*numberofnodes+1):endhour(i)*numberofnodes) for i in 1:365) |> collect),
            Global_Volume=((sum((working_df[s, :VolumeTotal]) for s in (endhour(i-1)*numberofnodes+1):endhour(i)*numberofnodes) for i in 1:365) |> collect),
            )
        hourly_global_DSOM_df = DataFrame(hour=1:8760, 
            Total_RD_Cost=((sum((working_df[s, :CostTotal]) for s in ((i-1)*numberofnodes+1:i*numberofnodes)) for i in 1:8760) |> collect),
            Global_Vol_Res=((sum((working_df[s, :VolumeRes]) for s in ((i-1)*numberofnodes+1:i*numberofnodes)) for i in 1:8760) |> collect),
            Global_Vol_PP=((sum((working_df[s, :VolumePP]) for s in ((i-1)*numberofnodes+1:i*numberofnodes)) for i in 1:8760) |> collect),
            Global_Volume=((sum((working_df[s, :VolumeTotal]) for s in ((i-1)*numberofnodes+1:i*numberofnodes)) for i in 1:8760) |> collect),
            )

        GlobalVol_Curtailment_BAU_hourly_twoweeks_df = subset(hourly_global_BAU_df[!, [:hour, :Global_Vol_Res]], :hour => x -> (321-1)*24+1 .<= x .<= 334*24)
        GlobalVol_Curtailment_BAU_hourly_twoweeks_df[!, 2] = GlobalVol_Curtailment_BAU_hourly_twoweeks_df[!, 2]*-1
        GlobalVol_Curtailment_DSOM_hourly_twoweeks_df = subset(hourly_global_DSOM_df[!, [:hour, :Global_Vol_Res]], :hour => x -> (321-1)*24+1 .<= x .<= 334*24)
        GlobalVol_Curtailment_DSOM_hourly_twoweeks_df[!, 2] = GlobalVol_Curtailment_DSOM_hourly_twoweeks_df[!, 2]*-1
        GlobalVol_Curtailment_TSOM_hourly_twoweeks_df = subset(hourly_global_TSOM_df[!, [:hour, :Global_Vol_Res]], :hour => x -> (321-1)*24+1 .<= x .<= 334*24)
        GlobalVol_Curtailment_TSOM_hourly_twoweeks_df[!, 2] = GlobalVol_Curtailment_TSOM_hourly_twoweeks_df[!, 2]*-1

        sum(GlobalVol_Curtailment_BAU_hourly_twoweeks_df.Global_Vol_Res)
        sum(GlobalVol_Curtailment_DSOM_hourly_twoweeks_df.Global_Vol_Res)
        sum(GlobalVol_Curtailment_TSOM_hourly_twoweeks_df.Global_Vol_Res)

        GlobalCost_BAU_hourly_twoweeks = subset(hourly_global_BAU_df[!, [:hour, :Total_RD_Cost]], :hour => x -> (321-1)*24+1 .<= x .<= 334*24)
        GlobalCost_DSOM_hourly_twoweeks = subset(hourly_global_DSOM_df[!, [:hour, :Total_RD_Cost]], :hour => x -> (321-1)*24+1 .<= x .<= 334*24)
        GlobalCost_TSOM_hourly_twoweeks = subset(hourly_global_TSOM_df[!, [:hour, :Total_RD_Cost]], :hour => x -> (321-1)*24+1 .<= x .<= 334*24)

        GlobalVol_BAU_Curtailment_yearlong = sum(daily_global_BAU_df.Global_Vol_Res)
        GlobalVol_DSOM_Curtailment_yearlong = sum(daily_global_DSOM_df.Global_Vol_Res)
        GlobalVol_TSOM_Curtailment_yearlong = sum(daily_global_TSOM_df.Global_Vol_Res)
        Annual_Curtailment_BYCASE_df = DataFrame(BAU=[GlobalVol_BAU_Curtailment_yearlong], DSOM=[GlobalVol_DSOM_Curtailment_yearlong], TSOM=[GlobalVol_TSOM_Curtailment_yearlong],)

        CSV.write(joinpath(processpath, Global, "Annual_Curtailment_BYCASE.csv"), Annual_Curtailment_BYCASE_df)

        CSV.write(joinpath(processpath, BAU, "GlobalVol_Curtailment_BAU_hourly_twoweeks.csv"), GlobalVol_Curtailment_BAU_hourly_twoweeks_df)
        CSV.write(joinpath(processpath, DSOM, "GlobalVol_Curtailment_DSOM_hourly_twoweeks.csv"), GlobalVol_Curtailment_DSOM_hourly_twoweeks_df)
        CSV.write(joinpath(processpath, TSOM, "GlobalVol_Curtailment_TSOM_hourly_twoweeks.csv"), GlobalVol_Curtailment_TSOM_hourly_twoweeks_df)

        CSV.write(joinpath(processpath, BAU, "GlobalCost_BAU_hourly_twoweeks.csv"), GlobalCost_BAU_hourly_twoweeks)
        CSV.write(joinpath(processpath, DSOM, "GlobalCost_DSOM_hourly_twoweeks.csv"), GlobalCost_DSOM_hourly_twoweeks)
        CSV.write(joinpath(processpath, TSOM, "GlobalCost_TSOM_hourly_twoweeks.csv"), GlobalCost_TSOM_hourly_twoweeks)

### Calculate total RD cost over two weeks

    Global_RD_Cost_twoweeks_df = DataFrame()
    working_df = CSV.read(joinpath(BAUpath, "C_V_RD_N__1-8760.csv"), DataFrame)
    Global_RD_Cost_twoweeks_df[!, :BAU] = [sum(subset(working_df, :time => x -> (321-1)*24+1 .<= x .<= 334*24).CostTotal)]
    working_df = CSV.read(joinpath(DSOMpath, "C_V_RD_N__1-8760.csv"), DataFrame)
    Global_RD_Cost_twoweeks_df[!, :DSOM] = [sum(subset(working_df, :time => x -> (321-1)*24+1 .<= x .<= 334*24).CostTotal)]
    working_df = CSV.read(joinpath(TSOMpath, "C_V_RD_N.csv"), DataFrame)
    Global_RD_Cost_twoweeks_df[!, :TSOM] = [sum(subset(working_df, :time => x -> (321-1)*24+1 .<= x .<= 334*24).CostTotal)]
    working_df=nothing
    (Global_RD_Cost_twoweeks_df.BAU[1]-Global_RD_Cost_twoweeks_df.DSOM[1])/Global_RD_Cost_twoweeks_df.BAU[1]
    (Global_RD_Cost_twoweeks_df.BAU[1]-Global_RD_Cost_twoweeks_df.TSOM[1])/Global_RD_Cost_twoweeks_df.BAU[1]


    CSV.write(joinpath(processpath, Global, "GlobalCost_compare_twoweeks.csv"), Global_RD_Cost_twoweeks_df)

    (Global_RD_Cost_twoweeks_df.BAU[1]-Global_RD_Cost_twoweeks_df.TSOM[1])/Global_RD_Cost_twoweeks_df.BAU[1]
    (Global_RD_Cost_twoweeks_df.BAU[1]-Global_RD_Cost_twoweeks_df.DSOM[1])/Global_RD_Cost_twoweeks_df.BAU[1]


daily_global_BAU_df = DataFrame(day=1:365, 
    Global_Cost=((sum((working_df[s, :CostTotal]) for s in (endhour(i-1)*numberofnodes+1):endhour(i)*numberofnodes) for i in 1:365) |> collect),
    Global_Vol_Res=((sum((working_df[s, :VolumeRes]) for s in (endhour(i-1)*numberofnodes+1):endhour(i)*numberofnodes) for i in 1:365) |> collect),
    Global_Volume=((sum((working_df[s, :VolumeTotal]) for s in (endhour(i-1)*numberofnodes+1):endhour(i)*numberofnodes) for i in 1:365) |> collect),
    )

sum(subset(working_df, :time => x -> (321-1)*24+1 .<= x .<= 334*24).CostTotal)


highcost_hourly_df = copy(working_df[numberofnodes*24*(highcoststartday-1)+1:(numberofnodes*24*(highcoststartday+13)), [:node, :time, :CostTotal]])
highcost_hourly_df = unstack(highcost_hourly_df, :time, :CostTotal)
highcost_hourly_df[!, :total] = ((sum(highcost_hourly_df[i, 2:end]) for i in 1:numberofnodes) |> collect)
sum(highcost_hourly_df[!, :total])

nodes_df[!, :totalcost] = highcost_hourly_df[!, :total]
nodes_df
nodes_sorted_df = sort(nodes_df, [:totalcost], rev=true)

sum(nodes_df[!, :totalcost])
sum(nodes_sorted_df[!, :totalcost])

daily_global_df


GlobalVol_RD_twoweeks_BYCASE_df = vcat(TwoWeekVol_RD_BAU, TwoWeekVol_RD_DSOM, TwoWeekVol_RD_TSOM)
GlobalVol_RD_twoweeks_BYCASE_df.Case = ["BAU","DSOM","TSOM"]
GlobalVol_RD_twoweeks_BYCASE_df
CSV.write(joinpath(processpath, Global, "GlobalVol_RD_twoweeks_BYCASE.csv"), GlobalVol_RD_twoweeks_BYCASE_df)

sum(subset(CM_LL_BAU_df, :output => x -> x .> 0).output)

))