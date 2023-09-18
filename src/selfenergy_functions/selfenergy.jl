using LinearAlgebra
using IterativeSolvers


#=======================================================================================#
"""Function that calculated the backward error for a solution of a Green's function
matrix GV for a lead.
Input: GV Green's function matrix of a lead
       H Hamilton matrix of a slice in lead
       V Coupling matrix between slices in lead
       E Energy
       eta imaginary shift to energy
Output: Forward error according to Eq. 6 

"""
function GreenCheck(GV,H,V,E,eta)

    try
        n=size(GV,1)
        I0=Matrix{ComplexF64}(I, n, n)
        return norm(V*GV*GV-((E+im*eta)*I0-H)*GV+V')
    catch error
        return missing
    end
end

#=======================================================================================#
"""Creating matrices A and B (see Eq. 14 and 15) and identity matrix I0 and
zero matrix O.
Input:  n size of matrices
        E Energy
        H Hamilton matrix of slice in lead
        V Coupling matrix between slices in lead
        etaR imaginary shift of energy
Output: Matrices A,B,I0,O
"""
function setup_problem(n,E,H,V,etaR) 


    #-----zero matrix and identity matrix constructed--#
    O=zeros(ComplexF64,n,n)
    I0=Matrix{ComplexF64}(I, n, n)
    
    #-----A and B matrices of the eigenvalue problem contructed--#
    A=[O I0;-V' ((E+etaR)*I0-H)]
    B=[I0 O; O V]
    
    return A,B,I0,O
end

#=======================================================================================#
"""Eigenvalues within 1e-10 offset of the unit circle the unit circle sorted into
inside and outside eigenvalues with the pertubation method.
Input:  n size of matrices
        F Output from the Schur factorization of the A and B matrices
        I0 Identity matrix
        O zero matrix
Output: P Boolean array containing information on if eigenvalue is inside or outside
        of unit circel
"""
function eigenorder_pert(n,F,I0,O)

    #------size of hamilton slice matrices -------#
    n=convert(Int,size(I0)[1])
    
    #-----Matrices for perturbation calculation--#
    I22=[O O; O I0]
    Tnn=[F.T[i,i] for i in 1:2*n]    
    dT=im*(F.Q)'*I22*F.Z
    dTnn=[dT[i,i] for i in 1:2*n]

    #-----All eigenvalues and offset from unit circle-#
    la=F.α./F.β    
    dla=Float64(1.0e-10)
    
    #-----Boolian mask for if eigenvalues are inside unit circle--#
    P=Array{Bool}(undef,2*n)
    for i=1:2*n
        #---Check if eigenvalue is NaN or outside of offset from unit circle--#
        if(isnan(la[i]) || (abs(la[i]) > (Float64(1.0)+dla)))
            P[i]=false
        #--Check if eigenvalue is inside of offset from unit circle--#
        elseif(abs(la[i])<(Float64(1.0)-dla))
            P[i]=true
            
        #--Here eigenvalue are within offset of unit circle--#
        else
           #--sign of perturbation is positive--#
           if(real(conj(la[i])*(dTnn[i]/F.β[i])) > Float64(0.0)) 
                P[i]=false
            else
                P[i]=true
            end
        end
    end
    return P
end
#=======================================================================================#
"""Eigenvalues within 1e-10 offset of the unit circle the unit circle sorted into
inside and outside eigenvalues with the standard method.
Input:  n size of matrices
        F Output from the Schur factorization of the A and B matrices
        A matrix (see Eq. 14)
        B matrix (see Eq. 15)
        K negative imaginary unit matrix
        Hsl Hamilton matrix of slice in lead
        Vsl Coupling matrix between slices in lead
        I0 Identity matrix
        O zero matrix
Output: P Boolean array containing information on if eigenvalue is inside or outside
        of unit circel
"""
function eigenorder_standard(n,F,A,B,K,Hsl,Vsl)


    #-----All eigenvalues and offset from unit circle-#
    eigenvalues=F.α./F.β 
    dla=Float64(1.0e-10)
    BT=transpose(B)
    AT=transpose(A)
    
    #----Check if any eigenvalue is within the offset from unit circle-------------#
    dla=Float64(1.0e-10)
    if any(abs.(norm.(eigenvalues).-1.0).<dla)##
        
        #----Boolean for eigenvalues within the offset of unit circle-----------------------------------------#
        P=Array{Bool}([(!isnan(x)) && (norm(x)<(1.0-dla)) for x in eigenvalues])
   
        #----Boolean for eigenvalues within offset of unit circle-----------------------------------------#
        Pb=Array{Bool}([((1.0-dla)<=norm(x)<=(1.0+dla)) for x in eigenvalues])

   
        #---The inverse Power method for generalized eigenvalue problems performed for eiginvalues on the unit circle----#
        la=eigenvalues[Pb]
        m=count(Pb)
        mask=fill(true,m) 
        qr=ones(ComplexF64,n)

        for i in collect(1:m)
            mu=la[i]
            H=Hsl+Vsl*mu+Vsl'*conj(mu)
            fill!(qr,1.0)
            _,qr=powm!(H,qr;shift=mu,inverse=true)
            
            mup=(qr'*K*qr) 
            if real(mu'*mup+mu*mup')>0 
                mask[i]=false
            end  
        end  
        #----Complete the Bool mask for for which eigenvalue are within the unit circle--# 
        Pb[Pb.==true]=mask
        P=P .| Pb  

    else  
        
        #----Bolean for eigenvalues within unit circle--------------------------------------#
        P=Array{Bool}([(!isnan(x)) && (norm(x)<1.0) for x in eigenvalues ])
    end  
end

#=======================================================================================#
""" Using Schur decomposition and eigenvalur ordering to calculate Green's function,
see Eq. 6 and Eq. 40.
Input: n size of matrice
       F0 output of schur decompostion of A and B matrices
       P Boolean array containing information on if eigenvalue is inside or outside
        of unit circel
Output: GV^+
"""
function calc_self_en(n,F0,P)

    F=ordschur(F0,P)
    enew=F.α./F.β

    Z11=F.Z[1:n,1:n]
    Z21=F.Z[n+1:2*n,1:n]
    
    return Z21*inv(Z11)
end
#=======================================================================================#
"""Calculates the solution to Eq. 6 using the perturbation method.
Input: E Energy
       H Hamilton matrix of a slice in lead
       V Coupling matrix between slices in lead
       etaR imaginary shift to energy
Output: GV^+, see Eq. 40.
"""
function selfenergy_pert(E,H,V,etaR)

    #------size of hamilton slice matrices -------#
    n=convert(Int,size(H)[1])
    A,B,I0,O=setup_problem(n,E,H,V,etaR)
    
    #------schur factorization--------------------#
    F=schur(A,B)
    
    #-----eigenvalues diveded into within and outside unit circle---#
    P=eigenorder_pert(n,F,I0,O)
    
    return calc_self_en(n,F,P)
end
#=======================================================================================#
"""Calculates the solution to Eq. 6 using the standard method.
Input: E Energy
       H Hamilton matrix of a slice in lead
       V Coupling matrix between slices in lead
       etaR imaginary shift to energy
Output: GV^+, see Eq. 40.
"""
function selfenergy_standard(E,H,V,etaR)

    #------size of hamilton slice matrices -------#
    n=convert(Int,size(H)[1])
    A,B,I0,O=setup_problem(n,E,H,V,etaR)
    
    #------schur factorization--------------------#
    F=schur(A,B) 
    
    #-----eigenvalues diveded into within and outside unit circle---#
    K=-im*I0
    P=eigenorder_standard(n,F,A,B,K,H,V)
    
    return calc_self_en(n,F,P)
end

#=======================================================================================#
function all_eigenvalues(E,H,V,etaR)

    #------size of hamilton slice matrices -------#
    n=convert(Int,size(H)[1])
    A,B,_,_=setup_problem(n,E,H,V,etaR)
    
    #------schur factorization--------------------#
    F=schur(A,B) 

    return F.α./F.β 
end

#=======================================================================================#

