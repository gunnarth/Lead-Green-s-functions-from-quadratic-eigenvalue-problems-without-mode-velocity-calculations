#===============Installing/Loading external packages==========#
using Pkg
try
    using DataFrames
    using CSV
    using IterativeSolvers
    println("Dataframes, CSV, and IterativeSolvers installed")
catch

    Pkg.add("DataFrames")
    Pkg.add("CSV")
    Pkg.add("IterativeSolvers")
    Pkg.precompile()
    println("Dataframes, CSV, and IterativeSolvers not installed")
    using DataFrames
    using CSV
    using IterativeSolvers
    println("Dataframes, CSV, and IterativeSolvers installed")
end

#=====Loading code=============================================#
DIR=@__DIR__
include(joinpath(DIR, "src", "hamiltonian_functions", "grapheneMatrices.jl"))
include(joinpath(DIR, "src", "hamiltonian_functions", "rashbaMatrices.jl"))
include(joinpath(DIR, "src", "selfenergy_functions", "selfenergy.jl"))
include(joinpath(DIR, "src", "selfenergy_functions",  "LopesSanchos.jl"))
include(joinpath(DIR, "src", "helper_functions",      "Timer.jl"))

#=======Setting Parameters=====================================#
etaList=[10.0^(-p) for p in [16,15,14,13]] #Imaginary shift to energy

#--list of energy points to calculate at
Energies=collect(0:0.05:10) #Rashba
Energies_ribbon=collect(-5:0.05:5) #Nanoribbon


repIt=10 #how many time to repeat calculation
N=24     #size parameter for matrices, corresponds to size 98

#--Hamiltonian matrices calculated and put in a dictionary
Hr,Vr=rashbaMatrices((2*N+1),1);
Hg,Vg=grapheneMatrices(N);

SystemType=Dict("Nanoribbon" => (Hg,Vg), "Rashba" => (Hr,Vr))

#========Dataframe for results created==========================#
dfE=DataFrame(E = Float64[],
              S = String[],
              eta = Float64[],
              GESpert_Min = Union{Float64,Missing}[],
              GESpert_FE = Union{Float64,Missing}[],
              GESstandard_Min = Union{Float64,Missing}[],
              GESstandard_FE = Union{Float64,Missing}[],
              LS_Min = Union{Float64,Missing}[],
              LS_FE = Union{Float64,Missing}[]);

#=======Iteration over system, eta and energy===================#
for key in keys(SystemType)
    for eta in etaList
        for (i,e) in enumerate(Energies)

            if key=="Nanoribbon"
                e=Energies_ribbon[i]
            end
            H=SystemType[key][1]
            V=SystemType[key][2]
            push!(dfE[!,:E],e)
            push!(dfE[!,:S],key)
            push!(dfE[!,:eta],eta)

            push!(dfE[!,[:GESpert_Min,:GESpert_FE]],
                        Timer(selfenergy_pert,e,H,V,eta,repIt,""))

            push!(dfE[!,[:GESstandard_Min,:GESstandard_FE]],
                        Timer(selfenergy_standard,e,H,V,eta,repIt,""))
 
            push!(dfE[!,[:LS_Min,:LS_FE]],
                        Timer(LS,e,H,V,eta,repIt,"LS"))  
                    
        end
    end
    CSV.write(joinpath(DIR,"Data","TimingData","Data_for_figs_5_and_6.csv"), dfE)
end