# Copyright 2019 Fabio Mazza 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

using StatsBase
using DataFrames
using CSV

function linspace(startN::Number,endN::Number,nel::Int)
    step::Float64 = (endN-startN)/(nel-1)
    collect(startN:step:endN)
end

function cycleindx(idx::Int,maxind::Int)
    if(idx>0 && idx<=maxind)
        return idx
    elseif(idx >maxind)
        return idx-maxind
    else
        return idx+maxind
    end
end

function meanVfromVStat(vstat::Array{Int64,2},lane::Int)
    z  = size(vstat,2)
    sum(vstat[lane,:].*collect(0:z-1))/sum(vstat[lane,:])
end
function calc_headway(bigpos::T,smallpos::T,L::T) where T<:Real
    ## Subtract after, so that if the two cars are superimposed the distance becomes negative
    d::T = bigpos-smallpos
    if(d<0)
        return d+L-1
    else
        return d-1
    end
end
function calc_dist(pos::Array{Int64,1},L::Int64;debug=false)
    # function to compute the distances between cars
    # GOOD FOR ONE LANE ONLY
    if(length(pos)>1)
        dist = calc_headway.(circshift(pos,-1),pos,L)
    else
        dist = pos
    end
    if(debug==true)
        println("calc_dist: ",dist)
    end
    return dist
end

##Generate the initial configuration
function gen_cars(L::Integer,prob::AbstractFloat)
    x=0
    cars_pos = zeros(Int,0)
    while true
        increm = trunc(Int,log(1-rand())/log(1-prob))
        while (increm) <= 0
            increm = trunc(Int,log(1-rand())/log(1-prob))
        end
        x += increm
        if x>L
            break
        end
        append!(cars_pos,x)
    end
    return cars_pos
end

function newVel(vcurr::Integer,vmax::Integer,head::Integer,pslow::Float64)
    max(0,min(vmax,head,vcurr+1)-Int(rand()<pslow))
end
function prop_rule(v,vmax,dist)
    v<=vmax & v>=dist
end
## Generate occupancy
function add_occup!(x::Int,lane::Int,vector::Array{Int,2})
    vector[lane,x] += 1
    return
end
function add_occup!(stat::Array{T,1},record::Array{T,2}) where T <: Real
    lane::T = stat[3]
    pos::T = stat[2]
    record[lane,pos] += 1
    return
end
function find_nearest(status::Array{Int,2},idxcent::Int,targetlane::Int)
    prec_found::Bool=false
    succ_found::Bool=false;
    N::Int = size(status,1)
    if(N<=1)
        return idxcent,idxcent
    elseif(N<3)
        prec_idx=cycleindx(idxcent-1,N);
        succ_idx = prec_idx;
        return prec_idx,succ_idx
    end
    prec_idx::Int=idxcent
    succ_idx::Int=idxcent
    for d=1:floor(Int,N/2)
        if(!prec_found)
            prec_idx = cycleindx(idxcent-d,N)
            if(status[prec_idx,3]==targetlane)
                prec_found=true
            end
        end
        if(!succ_found)
            succ_idx = cycleindx(idxcent+d,N)
            if(status[succ_idx,3]==targetlane)
                succ_found=true
            end
        end
        if(prec_found && succ_found)
            break
        end
    end
    prec_idx,succ_idx
end

function move_car_if_possible!(idx::Int,toLane::Int,status::Array{Int,2},
            changeCounter::Array{Int,1},L::Int)

    pos::Int = status[idx,2]
            ## Since cars are ordered by position, check the position of the neighboring cars on target lane
    beh_idx,for_idx = find_nearest(status,idx,toLane)
    d_prec::Int = calc_headway(pos,status[beh_idx,2],L)
    d_succ::Int = calc_headway(status[for_idx,2],pos,L)
            ## The distance is negative if the cars are on the same position (but different lane)
    if(status[beh_idx,4]<d_prec && status[idx,4]<d_succ)
        status[idx,3] = toLane
        changeCounter[toLane] += 1
    end
end
#=
Two LANES!!!
=#
function run_sim_double(p0::Float64,L::Int64,T::Int64,pslow::Float64,
        vmax::Int64,what::String,startw2::Bool=true)
    T_START_STATS=200
    
    
    if(startw2)
        #cars_px = gen_cars(2*L,p0)
        #sidx= findall(cars_px.>L)[1]
        cars_posl1 =  gen_cars(L,p0)#cars_px[1:sidx-1]
        #cars_posl2 = cars_px[sidx:end].-L
        cars_posl2 = gen_cars(L,p0)
    else
        cars_posl2 = Array{Int,1}[]
        cars_posl1 = gen_cars(L,p0)
    end
    N1::Int64 = size(cars_posl1,1)
    N2::Int64 = size(cars_posl2,1)
    #N = size(cars_pos,1)
    rho::Float64 = 0;
    
    if(startw2)
        rho = (N1+N2)/(L*2)
    else
        rho = (N1)/L
    end
    
    cars_mark_pos = zeros(Int,2,L)
    jumps = zeros(Int,2,T)

   
    N::Int64 = N1+N2
    status = zeros(Int,N,4)
    status[:,1] = collect(1:(N1+N2)) #Cars_ID
    status[:,2] = vcat(cars_posl1,cars_posl2) # Position
    status[:,3] = vcat(fill(1,N1),fill(2,N2)) #lane
    status[:,4] = fill(0,(N1+N2)) #speed
    
    
    statusrec = zeros(Int,T+1,N,4)
    
    ##REORDER ARRAY BY POSITION
    status = sortslices(status,dims=1,by=x->x[2])
    
    #println("Cars: ",N,"\t","Actual density: ",rho);
    statusrec[1,:,:] = status;
    alert = true;
    ngoback = zeros(Int,2)
    debug::Bool = false;verbose=false;
    vstat = zeros(Int,2,vmax+1)
    for t=1:T
        verbose ==true ? println("time $t") : 0;
        #divide by lane
        carsmoved = zeros(Int,2)
        ##temporary indices
        index_l1 = findall(status[:,3] .== 1)
        index_l2 = findall(status[:,3] .== 2)
        #CHECK no lost vehicules
        if(length(index_l1)+length(index_l2) < N)
            println("REALLY BAD ERROR")
        end
        if(debug == true)
            println("Status")
            println("\t",status[index_l1,1],"\t\t",status[index_l2,1])
            println("\t",status[index_l1,2],"\t\t",status[index_l2,2])
            println("\t",status[index_l1,4],"\t\t",status[index_l2,4])
        end
        ##MOVE FROM LANE 2 TO LANE 1
        #find distances USELESS
        #disttot = calc_dist(status[:,2],L)
        N_lane2::Int = size(index_l2,1)
        for idx2 in index_l2
            move_car_if_possible!(idx2,1,status,carsmoved,L)
           
        end
        ##MOVE FROM LANE 1 TO LANE 2
        #REBUILD indices
        #index_l1 = findall(status[:,3] .== 1)
        #index_l2 = findall(status[:,3] .== 2)
        distl1 = calc_dist(status[index_l1,2],L)
        if(debug==true)
            println("\tl1->l2");
            println(status[index_l1,4])
            println(distl1)
        end
        mask = prop_rule.(status[index_l1,4],vmax,distl1)
        checkarr = BitArray{1}(undef,length(mask))
        for i=1:length(mask)
            if(mask[i])
                caridx = index_l1[i]
                checkarr[i] = !any(status[index_l2,2].== status[caridx,2])
            end
        end
        prop_car = index_l1[checkarr]
         #this contains the indices of the cars that want to be brought to lane 2
        if(length(prop_car)>0 && debug == true)
            println("stpr",status[prop_car,:])
        end
        # check if the proposed car have enough space in lane 2, 
        # also considering the ones that will move to lane 1
        for idx in prop_car
            move_car_if_possible!(idx,2,status,carsmoved,L)
        end
        #=
        selcarst = sortslices(vcat(status[index_l2,:],status[prop_car,:]),dims=1,by=x->x[2])
        # calc distances
        if(size(selcarst,1)>0 && size(status[index_l2,:])==0)
            #no cars in the lane 2 -> put all the proposed car from lane 1 to lane 2
            status[prop_car,3] = 2
            carsmoved[2]+=length(prop_car)
        elseif(size(selcarst,1)>0)
            disttot = calc_dist(selcarst[:,2],L)
            lane1change = BitArray{1}(undef,length(prop_car))
            
            for i=1:length(prop_car)
                ##find again the prop_car
                newidx = findfirst(selcarst[:,1] .== status[prop_car[i],1])
                lessidx = cycleindx(newidx-1,size(disttot,1))
                #println(newidx,"  ",lessidx)
                if(selcarst[lessidx,4]<disttot[lessidx] && selcarst[newidx,4]<disttot[newidx])
                    ##move the car to lane 2
                    status[prop_car[i],3] = 2
                    carsmoved[2]+=1
                end
            end
            
        else
            debug == true ? println("timestep $t: no cars selected") : 1;
        end
        =#
        if(verbose == true)
            println("Finished lane changes;\n moved ",carsmoved[1]," cars from lane 2 to 1 \n and ",carsmoved[2]," from lane 1 to 2")
        end

        #REBUILD indices
        index_l1 = findall(status[:,3] .== 1)
        index_l2 = findall(status[:,3] .== 2)
        #CHECK no lost vehicules
        if(length(index_l1)+length(index_l2) < N)
            println("REALLY BAD ERROR")
        end
        jumps[:,t]=carsmoved
        if(debug == true)
            println("Status")
            println("\t",status[index_l1,1],"\t\t",status[index_l2,1])
            println("\t",status[index_l1,2],"\t\t",status[index_l2,2])
            println("\t",status[index_l1,4],"\t\t",status[index_l2,4])
        end
        ##MAIN MOVES OF Na-Sch
        ##compute distances

        for lid=1:2 ##for each lane
            lanemask = status[:,3] .== lid
            #lanearr
            dist = calc_dist(status[lanemask,2],L)
            N_lane = count(lanemask)
            #println(dist)
            #=
            # Do not save headway...
            =#
            new_vs = newVel.(status[lanemask,4],vmax,dist,pslow)
            debug == true ? println("lid: $lid, new vels = ",new_vs) : 1;
            if(any(new_vs.<0))
                println("THIS SHOULD NOT HAVE HAPPENED, time ",t)
            end
            status[lanemask,2] = status[lanemask,2] .+new_vs
            status[lanemask,4] = new_vs
            ngoback[lid] += count(status[lanemask,2].>L)
            #= DEBUGGING OVERLAPPING VEHICLES
            if(!checkOverlaps(status,L,t,true))
                println("ERROR AFTER NASCH LANE ",lid)
                println("VEHICLES: ",status[lanemask,:])
                println(oldstatus)
                throw(status)
            end
            oldstatus = copy(status)
            =#
            if(t>T_START_STATS)
                for v=0:vmax
                    vstat[lid,v+1] += sum(new_vs.== v)
                end
            end
        end
        
        indec = status[:,2].>L
        status[indec,2] .-=L
        if(debug == true)
            println("Status")
            println("\t",status[index_l1,1],"\t\t",status[index_l2,1])
            println("\t",status[index_l1,2],"\t\t",status[index_l2,2])
            println("\t",status[index_l1,4],"\t\t",status[index_l2,4])
        end
        
        status = sortslices(status,dims=1,by=x->x[2])
        if(t>T_START_STATS)
            for i=1:N
                add_occup!(status[i,:],cars_mark_pos)
            end
        end
        statusrec[t+1,:,:] = status
    end
    #println("$pslow -> $rho")
    if(what=="dataframe")
        return (rho,vmax,pslow,L,T,ngoback[1]/T,ngoback[2]/T,true,mean(jumps[1,:])/(L*T),
            mean(jumps[2,:])/(L*T),meanVfromVStat(vstat,1),meanVfromVStat(vstat,2))
    elseif(what=="stat")
        return cars_mark_pos,vstat,jumps
    else
        return statusrec
    end
end


meas2lane = DataFrame(rho=AbstractFloat[],vmax=Int[],pslow=AbstractFloat[],L=Int[],T=Int[],fluxl1=AbstractFloat[],fluxl2=AbstractFloat[],
    twolanes=Bool[],jumpL1=AbstractFloat[],jumpL2=AbstractFloat[],vmeanl1=AbstractFloat[],vmeanl2=AbstractFloat[])
L = 1000
T = 4000
pslow=0.1;
vmax = 5;
file = "doubleRuns_no2Start_renewd.csv"
#oldmisure = CSV.read(file)
oldmisure = copy(meas2lane)


#densarr = vcat(linspace(0.01,0.3,200),linspace(0.3,0.97,800));
densarr = linspace(0.01,0.99,1100)
#println(densarr)

println("STARTING")
Threads.@threads for i=1:size(densarr,1)
    for v=1:5
        #print("step ",i,"  ",v,"\r")
        #put!(ch2,(i,v))
        push!(meas2lane,run_sim_double(densarr[i],L,T,pslow,v,"dataframe",false))
        #yield()
    end
end

CSV.write(file,vcat(meas2lane,oldmisure))

