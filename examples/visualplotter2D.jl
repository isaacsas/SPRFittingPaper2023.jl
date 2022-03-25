using SPRFitting
using Plots
using DataFrames, Printf, CSV
using UnPack


reaches = [2.0,10.0,30.0]
kon     = 10^(-4)   # per concentration per time 
koff    = 0.05      # per time
konb    = 1.0       # per time
CP      = 1.0        

# A -> B rate in simulation is (kon * antibodyconcens[i])
antibodyconcens = [9.375,18.75,37.5,75.0,150.0,300.0]  
antigenconcen   = 0.0025   # L = sqrt(N/antigenconcen) with N=250 by default

# saving and plotting controls
plotsixpanelplot = false
savesixpanelplot = false
plottoprowplot   = true
savetoprowplot   = true
saveCSVofcurves  = false

# this saves to a subfolder named "figures and data" of the folder containing this script
savefolder       = joinpath(@__DIR__, "figures and data")
sixplotname      = joinpath(savefolder, "fullfig.pdf")
toprowplotname   = joinpath(savefolder, "toprow.pdf")
datacsvname      = joinpath(savefolder, "curves_for_bottom_row.csv")

############### plotting code 
function circleshape(h,k,r)
    theta = LinRange(0,2*pi,500)
    h .+ r*sin.(theta), k.+ r*cos.(theta)
end

function makeplots(bioparams, numparams, reach, antibodyconcens)
    @unpack initlocs,L,N,tstop = numparams
    bioparams.reach = reach

    SD = []    
    for (i,antibodyconcen) in enumerate(antibodyconcens)
        println("Running reach = $reach, [antibody] = $antibodyconcen ($i/$(length(antibodyconcens)))")
        bioparams.antibodyconcen = antibodyconcen
        bindcnt = run_spr_sim(bioparams, numparams)
        push!(SD, bindcnt)
    end

    # plot circles at initial locations
    MonoPhasic = []
    BiPhasic   = []    
    for i in 1:N
        mono=1
        for j in 1:N
            if i==j
            else
                if SPRFitting.sqdist_periodic(initlocs[:,i],initlocs[:,j],L) < reach^2
                    mono=0
                end
            end
        end
        if mono==1
            push!(MonoPhasic,i)
        else
            push!(BiPhasic,i)
        end
    end
    X = [initlocs[1,MonoPhasic],initlocs[1,BiPhasic]]
    Y = [initlocs[2,MonoPhasic],initlocs[2,BiPhasic]]
    circles = [circleshape(initlocs[1,i],initlocs[2,i],reach/2) for i in 1:N]
    Colors = ["blue" for i in 1:N]
    Colors = reshape(Colors,1,N) # casts to a matrix which is what plot() needs
    for a in BiPhasic
        Colors[a] = "red"
    end
    p1 = plot(circles,fillalpha=.5,seriestype=[:shape,],labels=:none,color=Colors,xlims=(-50,50),ylims=(-50,50),title="reach = $reach nm")

    # plots curves
    X2 = 0:tstop
    Y2 = [SD[i] for i in eachindex(antibodyconcens)]
    p2 = plot(X2,Y2,labels=:none,color="black",xlims=(0,tstop),ylims=(0,1))

    return p1,p2,X2,Y2
end

########### ACTUAL SCRIPT TO RUN SIMS AND MAKE FIGURES

# this stores the biological parameters and is passed around in the code
bioparams = BioPhysParams(; kon, koff, konb, reach=reaches[1], CP, antibodyconcen=antibodyconcens[1], antigenconcen)

# you can pass non-default numerical parameters as keywords to SimParams
# an example would be SimParams(bioparams.antigenconcen; dt=1.0, N=1000)
# to change dt and N from the defaults, see definition of SimParams above.
numparams = SimParams(bioparams.antigenconcen)

# run simulations and plot/save data
p1v = []; p2v = []; Xv = []; Yv = [];
idx = 1
for (i,reach) in pairs(reaches)
    p1,p2,X,Y = makeplots(bioparams, numparams, reaches[i], antibodyconcens)
    push!(p1v,p1); push!(p2v,p2); push!(Xv,X); push!(Yv,Y)
end
# create the output directory if it doesn't exist
if savesixpanelplot || saveCSVofcurves || savetoprowplot
    mkpath(savefolder)
end

# plots and optionally saves six panel plot
fullplot = plot(vcat(p1v,p2v)...,layout=(2,length(p1v)))
savesixpanelplot && savefig(fullplot, sixplotname)   

# plots and optionally saves top row of the preceding plot
pmerge = plot(p1v...,layout=(1,length(p1v)), size=(900,300),xtickfontsize=12,ytickfontsize=12)
savetoprowplot && savefig(pmerge, toprowplotname)   

# saves data from bottom row curves in CSV
if saveCSVofcurves
    data = zeros(length(Xv[1]),length(antibodyconcens)*length(reaches)+1)
    data[:,1] .= Xv[1]
    idx = 2
    paramnames = ["times"]
    for (i,reach) in pairs(reaches)
        kon = bioparams.kon
        for (j,abconcen) in pairs(antibodyconcens)
            data[:,idx] = Yv[i][j]
            idx += 1
            kon *= abconcen
            konstr = @sprintf "reach=%2.0f, kon=%2.5f" reach kon
            push!(paramnames, konstr)
        end
    end
    df = DataFrame(data,paramnames)
    CSV.write(datacsvname,df)
end

if plotsixpanelplot
    display(fullplot)
elseif plottoprowplot
    display(pmerge)
end