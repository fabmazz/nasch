
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
Normal simulation
=#
function newVel(vcurr::Integer,vmax::Integer,head::Integer,pslow::Float64)
    max(0,min(vmax,head,vcurr+1)-Int(rand()<pslow))
end

function run_NS(p0::Float64,L::Int64,T::Int64,pslow::Float64,
        vmax::Int64,outCSV::Bool)

    cars_pos = gen_cars(L,p0)
    N::Int64 = size(cars_pos,1)
    rho::Float64 = N/L

    vmaxarr = fill(vmax,N)

    vels = zeros(Int64,N)
    veh_pos = zeros(Int64,T,N)
    velrec = zeros(Int64,T,N)
    #distory = zeros(Int,T,N)  ##do not memorize DeltaX
    headwayc = zeros(Int64,L) ##to store the headaway distribution
    pos0 = cars_pos;
    #println("Cars: ",N,"\t","Actual density: ",rho);
    alert::Bool = true;
    ngoback::Int64 = 0
    for t=1:T
        ##compute distances
        dist = circshift(cars_pos,-1)-cars_pos.-1
        dist[findall(dist.<0)] .+=L
        #println(dist)
        if(outCSV == false)
            for x=1:L
                headwayc[x] += count(dist.==(x-1)) ##correct for non catching zero headaway
            end
        end
        #distory[t,:] = dist
       
        new_vs = newVel.(vels,vmaxarr,dist,pslow)
        #=
        
            etas = Int.(rand(N).<pslow);
            vtemp = min.(vmaxarr,dist,vels.+1)
            new_vs = max.(zeros(N),vtemp-etas)
        
        =#
        if(alert && any(new_vs.<0))
            print(t)
            println("THIS SHOULD NOT HAVE HAPPENED \n VTEMP:")
        end
        cars_pos = cars_pos.+new_vs
        vels= new_vs
        ngoback += count(cars_pos.>L)
        cars_pos[findall(cars_pos.>L)] .-=L
        veh_pos[t,:] = cars_pos
        velrec[t,:] = vels
    end
    #println("$pslow -> $rho")
    if(outCSV)
        return (rho,vmax,pslow,ngoback/T)#,veh_pos,velss)
    else 
        return veh_pos,velrec,headwayc
    end
end

misure = DataFrame(rho=AbstractFloat[],vmax=Int[],
    pslow=AbstractFloat[],flux=AbstractFloat[])


L = 1000
T = 4000
pslow=0.25;
file = "normal_p0.25.dat"
densarr = linspace(0.01,0.15,200);

oldmisure = CSV.read(file)
#oldmisure =DataFrame(misure)

for i=1:size(densarr,1)
    print("step $i\r")
    for v=1:5
        push!(misure,run_NS(densarr[i],L,T,pslow,v,true))
    end
end

CSV.write(file,vcat(misure,oldmisure))
