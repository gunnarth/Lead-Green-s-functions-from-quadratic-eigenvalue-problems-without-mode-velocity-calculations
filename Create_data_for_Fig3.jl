#===============Installing/Loading external packages==========#
using Pkg
try
    using DataFrames
    using CSV
    using IterativeSolvers
    println("Dataframes,CSV, and IterativeSolvers installed")
catch

    Pkg.add("DataFrames")
    Pkg.add("CSV")
    Pkg.add("IterativeSolvers")
    Pkg.precompile()
    println("Dataframes,CSV, and IterativeSolvers installed")
    using DataFrames
    using CSV
    using IterativeSolvers
end

#=====Loading code=============================================#
DIR=@__DIR__
include(joinpath(DIR, "src", "hamiltonian_functions", "grapheneMatrices.jl"))
include(joinpath(DIR, "src", "hamiltonian_functions", "rashbaMatrices.jl"))
include(joinpath(DIR, "src", "selfenergy_functions", "selfenergy.jl"))
include(joinpath(DIR, "src", "selfenergy_functions",  "LopesSanchos.jl"))
include(joinpath(DIR, "src", "helper_functions",      "Timer.jl"))
include(joinpath(DIR, "src", "helper_functions",      "CountBands.jl"))

#=======Setting Parameters=====================================#
eta=1e-14
N=24
Energies=collect(0:0.001:1);
repIt=10
Hg,Vg=grapheneMatrices(N);

#========Dataframe for results created==========================#
dfE=DataFrame(E = Float64[],
              N = Int64[],
              GESpert_Min = Union{Float64,Missing}[],GESpert_FE = Union{Float64,Missing}[],
              GESstandard_Min = Union{Float64,Missing}[],GESstandard_FE = Union{Float64,Missing}[],
              nBands = Union{Float64,Missing}[]);

#=======Iteration over energy===================#
for (i,e) in enumerate(Energies)
     push!(dfE[!,:E],e)
     push!(dfE[!,:N],2*(2*N+1)) 
             
     nBands=CountBands(e,Hg,Vg,eta)
     push!(dfE[!,:nBands],nBands)
           #
           push!(dfE[!,[:GESpert_Min,:GESpert_FE]],
                     Timer(selfenergy_pert,e,Hg,Vg,eta,repIt,""))
              
           #
           push!(dfE[!,[:GESstandard_Min,:GESstandard_FE]],
                     Timer(selfenergy_standard,e,Hg,Vg,eta,repIt,""))            
end
CSV.write(joinpath("Data","TimingData","Data_for_fig3_test.csv"),dfE)