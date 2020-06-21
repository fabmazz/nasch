# Copyright 2019 Fabio Mazza 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

using StatsBase
using DataFrames
using CSV

##substitute for numpy's linspace
function linspace(startN::Number,endN::Number,nel::Int)
    step::Float64 = (endN-startN)/(nel-1)
    collect(startN:step:endN)
end

function cycleindx(idx::Int,maxind::Int)
    ##Keep index in the same range
    if(idx>0 && idx<=maxind)
        return idx
    elseif(idx >maxind)
        return idx-maxind
    else
        return idx+maxind
    end
end
function calc_dist(pos::Array{Int64,1},L::Int64;debug=false)
    # function to compute the distances between cars
    if(length(pos)>1)
        dist = circshift(pos,-1)-pos.-1
        dist[findall(dist.<0)] .+=L
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

#=
Two LANES!!!
=#
function run_sim_double(p0::Float64,L::Int64,T::Int64,pslow::Float64,
        vmax::Int64)

    cars_pos = gen_cars(L,p0)
    N = size(cars_pos,1)
    rho::Float64 = N/L

    cars_lane = fill(1,N)


    status = zeros(Int,N,4)
    status[:,1] = collect(1:N) #Cars_ID
    status[:,2] = cars_pos # Position
    status[:,3] = fill(1,N) #lane
    status[:,4] = fill(0,N) #speed


    statusrec = zeros(Int,T+1,N,4)

    #println("Cars: ",N,"\t","Actual density: ",rho);
    statusrec[1,:,:] = status;
    alert = true;
    ngoback = zeros(Int,2)
    debug::Bool = false;verbose=false;
    for t=1:T
        verbose ==true ? println("time $t") : 0;
        #divide by lane
        carsmoved = zeros(Int,2)
        ##temporary indeces
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
        #find distances, putting all cars in the same lane
        disttot = calc_dist(status[:,2],L)
        for idx2 in index_l2
            ##check speeds
            if (idx2-1>0) #idx2-1 is the preceding car
                if(status[idx2-1,4]<disttot[idx2-1] && status[idx2,4]<disttot[idx2])
                    #change lane
                    status[idx2,3] = 1
                    #lane2change[idx2] = true
                    carsmoved[1] += 1
                end
            else
                if(status[idx2-1+N,4]<disttot[idx2-1+N] && status[idx2,4]<disttot[idx2])
                    #change lane
                    status[idx2,3] = 1
                    carsmoved[1] += 1
                end
            end
        end
        ##MOVE FROM LANE 1 TO LANE 2
        distl1 = calc_dist(status[index_l1,2],L)
        if(debug==true)
            println("\tl1->l2");
            println(status[index_l1,4])
            println(distl1)
        end
        prop_car = index_l1[status[index_l1,4].<=vmax .& status[index_l1,4].>=distl1]
         #this contains the indices of the cars that want to go to lane 2
        if(length(prop_car)>0 && debug == true)
            println("stpr",status[prop_car,:])
        end
        # check if the proposed car have enough space in lane 2,
        # also considering the ones that will move to lane 1
        selcarst = sortslices(vcat(status[index_l2,:],status[prop_car,:]),dims=1,by=x->x[2])
        # calc distances
        if(size(selcarst,1)>0 && size(status[index_l2,:])==0)
            #no cars in the lane 2 -> put all the proposed car from lane 1 to lane 2
            status[prop_car,3] = 2
            carsmoved[2]+=length(prop_car)
        elseif(size(selcarst,1)>0)
            ##there are cars in lane 2
            # this represents the distance in the new lane configuration
            disttot = calc_dist(selcarst[:,2],L)
            lane1change = BitArray{1}(undef,length(prop_car))

            for i=1:length(prop_car)
                ##find again the prop_car (idx)
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
            #lstat = status[lanemask,:]
            dist = calc_dist(status[lanemask,2],L)
            N_lane = count(lanemask)
            #println(dist)
            #=
            # Do not save headway...
            #for x=1:L
            #    headwayc[x] += count(dist.==(x-1)) ##correct for non catching zero headaway
            #end
            =#
            new_vs = newVel.(status[lanemask,4],vmax,dist,pslow)
            debug == true ? println("lid: $lid, new vels = ",new_vs) : 1;
            if(any(new_vs.<0))
                print(t)
                println("THIS SHOULD NOT HAVE HAPPENED \n VTEMP:")
            end
            status[lanemask,2] = status[lanemask,2] .+new_vs
            status[lanemask,4] = new_vs
            ngoback[lid] += count(status[lanemask,2].>L)
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
        statusrec[t+1,:,:] = status
    end
    #println("$pslow -> $rho")
    return (rho,vmax,pslow,ngoback[1]/T,ngoback[2]/T,true)
    #statusrec
end

meas2lane = DataFrame(rho=AbstractFloat[],vmax=Int[],
    pslow=AbstractFloat[], fluxl1=AbstractFloat[],fluxl2=AbstractFloat[],twolanes=Bool[])
L = 1000
T = 3000
pslow=0.1;
vmax = 5;
file = "twoLanesandFlux.csv"
#oldmisure = CSV.read(file)
oldmisure = DataFrame(meas2lane)


#densarr = vcat(linspace(0.01,0.3,200),linspace(0.3,0.97,800));

densarr = linspace(0.01,0.999,200);

for i=1:size(densarr,1)
    for v=1:5
        print("step ",i," v ",v,"\r")
        push!(meas2lane,run_sim_double(densarr[i],L,T,pslow,v))
    end
end

CSV.write(file,vcat(meas2lane,oldmisure))
