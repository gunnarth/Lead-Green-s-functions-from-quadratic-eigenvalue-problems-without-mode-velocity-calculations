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
#using ProgressMeter
using CSV
DIR=@__DIR__
println(joinpath(DIR, "src","helper_functions","CodeLoader.jl"))
include(joinpath(DIR, "src","helper_functions","CodeLoader.jl"))



etaList=[10.0^(-p) for p in [14,]] #[6,10,12,14]]
#Energies=collect(0:0.5:100);
Energies=collect(0:1:100);

#Energies=collect(0:100:2000);
#sort!(push!(Energies,5,10,25,50))
#Energies_ribbon=collect(-5:0.5:5);
Energies_ribbon=collect(0:0.001:1);


println(length(Energies),", ",length(Energies_ribbon))



#repIt=20
repIt=10
for N in [12,24,48]
    println("======== N=",N,"==============")
    for ML in  [100.0,]# 1e2]#[0.03,0.3,1.0,1.0e2]#[1.0e2, 1.0]
        #Hg,Vg=grapheneMatrices_1Dinner(2*N);
        Hg,Vg=grapheneMatrices(N);
        println("Hg:",size(Hg),", Vg:",size(Vg))

        #Hr,Vr=rashbaMatrices((2*N+1),1,ML);
        #println("Hr:",size(Hr),", Vr:",size(Vr))
        

        println("================================")
        SystemType=Dict("Nanoribbon" => (Hg,Vg))#, "Rashba" => (Hr,Vr))
        #SystemType=Dict("Rashba" => (Hr,Vr))


        for key in keys(SystemType)
            println("====================================",key,"====================================")
            dfE=DataFrame(E = Float64[],
            N = Int64[],
            S = String[],
            eta = Float64[],
            GESpert_Min = Union{Float64,Missing}[], GESpert_Mean = Union{Float64,Missing}[], GESpert_Std = Union{Float64,Missing}[], GESpert_FE = Union{Float64,Missing}[],
            GESpower_Min = Union{Float64,Missing}[], GESpower_Mean = Union{Float64,Missing}[], GESpower_Std = Union{Float64,Missing}[], GESpower_FE = Union{Float64,Missing}[],
            GESpower1_Min = Union{Float64,Missing}[], GESpower1_Mean = Union{Float64,Missing}[], GESpower1_Std = Union{Float64,Missing}[], GESpower1_FE = Union{Float64,Missing}[],
            GESpower2_Min = Union{Float64,Missing}[], GESpower2_Mean = Union{Float64,Missing}[], GESpower2_Std = Union{Float64,Missing}[], GESpower2_FE = Union{Float64,Missing}[],
            GESnaive_Min = Union{Float64,Missing}[], GESnaive_Mean = Union{Float64,Missing}[], GESnaive_Std = Union{Float64,Missing}[], GESnaive_FE = Union{Float64,Missing}[],
            GEScombined_Min = Union{Float64,Missing}[], GEScombined_Mean = Union{Float64,Missing}[], GEScombined_Std = Union{Float64,Missing}[], GEScombined_FE = Union{Float64,Missing}[],
           LS_Min = Union{Float64,Missing}[], LS_Mean = Union{Float64,Missing}[], LS_Std = Union{Float64,Missing}[],LS_FE = Union{Float64,Missing}[], 
            MagLen =Float64[],
            BE = Union{Float64,Missing}[],
            WE = Union{Float64,Missing}[],
            nBands = Union{Float64,Missing}[]);
            


            for eta in etaList
                for (i,e) in enumerate(Energies_ribbon)
                  #if key=="Nanoribbon"
                  #      e=Energies_ribbon[i]
                  #end
                  H=SystemType[key][1]
                  V=SystemType[key][2]
                  push!(dfE[!,:E],e)
                  push!(dfE[!,:N],2*(2*N+1))
                  push!(dfE[!,:S],key)
                  push!(dfE[!,:eta],eta)
                  BE=BackwardError(e,H,V)    
                    
                  nBands=CountBands(e,H,V,eta)
                  push!(dfE[!,:nBands],nBands)

                  #
                  push!(dfE[!,[:GESpert_Min,:GESpert_Mean,:GESpert_Std,:GESpert_FE]],
                            Timer(selfenergy_pert,e,H,V,eta,repIt,""))
                     
                  #
                  push!(dfE[!,[:GESpower_Min,:GESpower_Mean,:GESpower_Std,:GESpower_FE]],
                            Timer(selfenergy_power,e,H,V,eta,repIt,""))
                  #
                  push!(dfE[!,[:GESpower1_Min,:GESpower1_Mean,:GESpower1_Std,:GESpower1_FE]],
                            Timer(selfenergy_power1,e,H,V,eta,repIt,""))
                  #
                  push!(dfE[!,[:GESpower2_Min,:GESpower2_Mean,:GESpower2_Std,:GESpower2_FE]],
                            Timer(selfenergy_power2,e,H,V,eta,repIt,""))
                  #
                  push!(dfE[!,[:GESnaive_Min,:GESnaive_Mean,:GESnaive_Std,:GESnaive_FE]],
                            Timer(selfenergy_naive,e,H,V,eta,repIt,""))
                    
                  #
                  push!(dfE[!,[:GEScombined_Min,:GEScombined_Mean,:GEScombined_Std,:GEScombined_FE]],
                            Timer(selfenergy_combined,e,H,V,eta,repIt,""))
                
                  push!(dfE[!,[:LS_Min,:LS_Mean,:LS_Std,:LS_FE]],
                            Timer(leadGF,e,H,V,eta,repIt,"LS"))  
                    
                    
                    
                  #Backward skekkja fyrir eigingildi nálægt einingahring
                  push!(dfE[!,:BE],BE)
                  
                  #Skekkja skv. Wimmer
                  #push!(dfE[!,:WE],WimmerError(e,H,V,eta))
                  push!(dfE[!,:WE],0.0)
                    
                  #Segullengd
                  push!(dfE[!,:MagLen],ML)
                    
                end
            end

            println("writeout")
            println("N$(2*(2*N+1))_$(key)_MagLen$(ML)")

            CSV.write("Tímatökuniðurstöður-2/N$(2*(2*N+1))_$(key)_MagLen$(ML).csv", dfE)
        end        
    end
end
