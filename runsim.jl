# Copyright 2019 Fabio Mazza 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

using PyPlot
using StatsBase
using DataFrames
using CSV

##substitute for numpy's linspace
function linspace(startN::Number,endN::Number,nel::Int)
    step::Float64 = (endN-startN)/(nel-1)
    collect(startN:step:endN)
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

#=
Make the simulation with normal cars and trucks
=#

function run_sim_truck(p0::Float64,L::Int64,T::Int64,pslow::Float64,
        vmax::Array{Int64,1},ptruck::Float64)
    
    cars_pos = gen_cars(L,p0)
    N = size(cars_pos,1)
    rho::Float64 = N/L
    
    truckval::BitArray{1} = rand(N).<ptruck
    ntrucks = count(truckval)
    deltav::Float64=maximum(vmax)-minimum(vmax)    
    vmaxarr::Array{Int64,1} = Int.(truckval).*(-1)*deltav .+maximum(vmax);
   
    vels = zeros(Int,N)
    veh_pos = zeros(Int,T,N)
    velrec = zeros(Int,T,N)
    #distory = zeros(Int,T,N)  ##do not memorize DeltaX
    headwayc = zeros(Integer,L) ##to store the headaway distribution
    pos0 = cars_pos;
    #println("Cars: ",N,"\t","Actual density: ",rho);
    alert = true;
    ngoback::Int = 0
    for t=1:T
        ##compute distances
        dist = circshift(cars_pos,-1)-cars_pos.-1
        dist[findall(dist.<0)] .+=L
        #println(dist)
        for x=1:L
            headwayc[x] += count(dist.==(x-1)) ##correct for non catching zero headaway
        end
        #distory[t,:] = dist
        etas = Int.(rand(N).<pslow);
        #etas = [draw_ps([1-pslow,pslow]) for i=1:N].-1;
        #etas .= 0
        vtemp = min.(vmaxarr,dist,vels.+1)

        #print(size(vtemp),size(etas))

        new_vs = max.(zeros(N),vtemp-etas)
        if(alert && any(new_vs.<0))
            print(t)
            println("THIS SHOULD NOT HAVE HAPPENED \n VTEMP:")
            #println(vtemp)
            #println("ETAS")
            #println(etas)
            #println("new vs:")
            #println(new_vs)
        end
        cars_pos = cars_pos.+new_vs
        vels= new_vs
        ngoback += count(cars_pos.>L)
        cars_pos[findall(cars_pos.>L)] .-=L
        veh_pos[t,:] = cars_pos
        velrec[t,:] = vels
    end
    #println("$pslow -> $rho")
    return (rho,maximum(vmaxarr),minimum(vmaxarr),pslow,ptruck,ngoback/T)#,veh_pos,velss)
    #veh_pos,velrec
end

misure = DataFrame(rho=AbstractFloat[],vmaxup=Int[],vmaxdown=Int[],
    pslow=AbstractFloat[],ptruck=AbstractFloat[], flux=AbstractFloat[])

p0 = 0.1
ptruck = 0.1;
L = 300
T = 500
pslow=0.1;
file = "trucksim.dat"
ptruck_arr =collect(0.1:0.1:0.6)
densarr = linspace(0.1,0.7,100);
oldmisure = CSV.read(file)

for i=1:size(densarr,1)
    print("step $i\r")
    for j=1:size(ptruck_arr,1)
        push!(misure,run_sim_truck(densarr[i],L,T,pslow,[2,4],ptruck_arr[j]))
    end
end

CSV.write(file,vcat(misure,oldmisure))
