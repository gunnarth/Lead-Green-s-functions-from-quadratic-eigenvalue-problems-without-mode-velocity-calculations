using LinearAlgebra

function GenEigenSelfAlaSchur_SIEmethod(E,H,V,etaR=1e-14,printEig=false)
    #------ Constructing the eigenvalue problem -------#
    n=convert(Int,size(H)[1])
    O=zeros(ComplexF64,n,n)
    I0=Matrix{ComplexF64}(I, n, n) 

    A=[O I0;-V' (E*I0-H)]
    B=[I0 O; O V]

    #---Generlized Schur factoring performed-------#
    F0=schur(A,B)
    eigenvalues=F0.α ./ F0.β # Eigenvalues ordered is ascending order

# Þessi hluti er óþarfi    
#    #---Greens function created with naive ordering--#
#    P=Array{Bool}([(!isnan(x)) && (norm(x)<1.0) for x in eigenvalues ])
#    F=ordschur(F0,P)
#    enew=F.α./F.β
#    
# Þessi hluti er óþarfi
    
    I22=[O O; O I0]
  
#    T=Q'*A*Z  
    Tnn=[F0.T[i,i] for i in 1:2*n]    
    dT=im*(F0.Q)'*I22*F0.Z
    dTnn=[dT[i,i] for i in 1:2*n]
  
#    etaNum=1.0e-11
#    P=Array{Bool}(undef,2*n)
#    for i=1:2*n
#        la=(F0.α[i]+etaNum*dTnn[i])/F0.β[i]  # Ekki góð leið því dTnn getur verið lítil stærð.
#        if(!isnan(la) && (abs(la)< 1.0))
#            P[i]=true
#        else
#            P[i]=false
#        end
#        println("dT=",abs(eigenvalues[i])," , ",abs(dT[i,i]))
#    end
    
    dla=Float64(1.0e-8)
    P=Array{Bool}(undef,2*n)
    for i=1:2*n
        if(isnan(eigenvalues[i]) || (abs(eigenvalues[i])>(Float64(1.0)+dla)))
            P[i]=false
        elseif(abs(eigenvalues[i])<(Float64(1.0)-dla))
            P[i]=true
        else
           if(real(conj(eigenvalues[i])*(dTnn[i]/F0.β[i]))>Float64(0.0))
                P[i]=false
                #println("|la0|,real(dT/S0),P[i]=",abs(eigenvalues[i])," , ",real(dTnn[i]/F0.β[i])," , ",P[i])
                #println(dTnn[i]," , ",F0.β[i])
            else
                P[i]=true
                #println("|la0|,real(dT/S0),P[i]=",abs(eigenvalues[i])," , ",real(dTnn[i]/F0.β[i])," , ",P[i])
                #println(dTnn[i]," , ",F0.β[i])
            end
        end
    end
      
    F=ordschur(F0,P)
#    println("neweig=",abs.(F.α./F.β))
    Z11=F.Z[1:n,1:n]
    Z21=F.Z[n+1:2*n,1:n]
    Gout=Z21*inv(Z11)
 
end