using SPRFitting
using Roots, Interpolations
using Plots, DataFrames, XLSX

########### PARMETERS TO MODIFY ###########

antigenconcens  = [0.2*(1/5)^(i) for i in 9:-1:1]
antibodyconcens = 10.0 .^ range(-2,6,step=8/29)

#FD-11A
kon = 1.027E-04
koff = 0.042
konb =  0.764
reach = 32.766

#FD-5D
# kon = 4.93E-05
# koff = 4.22E-03
# konb = 2.04E+00
# reach = 3.48E+01

#REGN10987
# kon = 8.50E-04
# koff = 3.79E-02
# konb = 3.06E+00
# reach = 4.02E+01

plotfigure  = true
savefigure  = true
savecsv     = true
savefolder  = joinpath(@__DIR__, "figures and data")
figname     = joinpath(savefolder, "DoseResponse.pdf")
xlsxname    = joinpath(savefolder, "DoseResponseData.xlsx")

tstop  = 15.0*60.0       # time in seconds
dt     = tstop/1000      # times at which we save the data
N      = 1000            # number of particles (tethers)
nsims  = 20

########### END USER PARAMETERS ###########

########### INTERNAL PARAMETER STRUCTURES ###########
# this stores the biological parameters and is passed around in the code
# we update this in the simulator for the current antibody and antigen concentrations
biopars = BioPhysParams(; kon, koff, konb, reach)

# this will need to get updated in the code as the domain length changes
# as the antigenconcentration changes
numpars = SimParams(biopars.antigenconcen[1]; tstop, dt, N, nsims)

########### END INTERNAL PARAMETER STRUCTURES ###########

function getIC50(freeantigen, antibodyconcens)

    itp     = interpolate(freeantigen[:,1], BSpline(Linear()))
    abrange = range(log10(antibodyconcens[1]),log10(antibodyconcens[end]); length=length(antibodyconcens))
    itp_x   = interpolate(abrange, BSpline(Linear()))

    f(x) = itp(x)-0.5

    z = find_zero(f,(1,30),Bisection())

    return 10^(itp_x(z))
    
end

# note that this modifies biopars and numpars!!!
function runpotencysims!(biopars, numpars, antigenconcen, antibodyconcens)
    # we just save the value at the final time point
    freeantigensims = zeros(Float64, length(antibodyconcens))
    for (ax,antibodyconcen) in enumerate(antibodyconcens)
        numpars.L              = (nsims/antigenconcen)^(1/2) 
        biopars.antigenconcen  = antigenconcen
        biopars.antibodyconcen = antibodyconcen    
        freeantigen            = run_spr_sim(biopars, numpars)
        freeantigensims[ax]    = freeantigen[end]
    end
    return freeantigensims
end

function varyantigenconcen(biopars, numpars, antigenconcens, antibodyconcens)
    freeantigenstore = zeros(Float64,30,length(antigenconcens))
    for (ax,antigenconcen) in enumerate(antigenconcens)
        @time freeantigensims   = runpotencysims!(biopars, numpars, antigenconcen, antibodyconcens)
        freeantigenstore[:,ax]  = freeantigensims
    end
    return freeantigenstore
end

################### ACTUAL SCRIPT TO RUnsimsSIMS AND PLOT %%%%%%%%%%%%%%%%%%%%%%%%

@time freeantigenstore = varyantigenconcen(biopars, numpars, antigenconcens, antibodyconcens)

if plotfigure
    X=[antibodyconcens for antigenconcen in antigenconcens]
    Y=[freeantigenstore[:,ax] for ax in 1:length(antigenconcens)]
    Labels=["$(antigenconcen)" for antigenconcen in antigenconcens]
    plt = plot(X,Y,seriestype=:line,label=:none,xscale=:log10,
               xticks=[1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5],
               yticks=[0,.25,.5,.75,1.0].
               xlabel = "antibody concentration nM",
               ylabel = "Free antigen")

    display(plt)

    savefigure && savefig(plt,figname)
end

if savecsv
    headers = [Symbol("Antigen Concentration $(antigenconcen) nM") for antigenconcen in antigenconcens]
    df = DataFrame(freeantigenstore,headers)

    # Write the DataFrame to an xlsx file
    XLSX.writetable(xlsxname,collect(eachcol(df)),names(df))
end