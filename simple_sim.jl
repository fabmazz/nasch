
using PyPlot
using StatsBase
#using DataFrames
#using CSV

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

function run_simple(p0::Float64,L::Int64,T::Int64,pslow::Float64,vmax::Int64)
    
    cars_pos = gen_cars(L,p0)
    N = size(cars_pos,1)
    rho::Float64 = N/L
    
    vmaxarr = fill(vmax,N)
   
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
        end
        cars_pos = cars_pos.+new_vs
        vels= new_vs
        ngoback += count(cars_pos.>L)
        cars_pos[findall(cars_pos.>L)] .-=L
        veh_pos[t,:] = cars_pos
        velrec[t,:] = vels
    end
    #println("$pslow -> $rho")
    return veh_pos,velrec#,veh_pos,velss)
    #veh_pos,velrec
end

p = 0.4
L = 200
T = 400
pslow = 0.1
vmax = 5
xrec,vrec = run_simple(p,L,T,pslow,vmax)

pcolor(vrec)
colorbar()
show()

