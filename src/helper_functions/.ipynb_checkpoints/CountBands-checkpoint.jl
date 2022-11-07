using LinearAlgebra

function CountBands(E,H,V,etaR=1e-14)
    #------ Constructing the eigenvalue problem -------#
    n=convert(Int,size(H)[1])
    O=zeros(ComplexF64,n,n)
    I0=Matrix{ComplexF64}(I, n, n) 
    #C=(-(E+5.0e-16)*I0+H)
    C=(-(E+etaR)*I0+H)

    K=V'
    M=V
    A=[O I0;-K -C];

    B=[I0 O; O M];
    AT=transpose(A)
    BT=transpose(B)
    K=[O  O; O -im*I0];

    #---Generlized Schur factoring performed-------#
    F=schur(A,B)
    #eigenvalues=F[:alpha]./F[:beta]
    eigenvalues=F.α./F.β 
    de=1.0e-10
    if any(abs.(norm.(eigenvalues).-1.0).<de)
        #----Eigenvalues on the unit circle-----------------------------------------#
        Pb=Array{Bool}([((1.0-de)<=norm(x)<=(1.0+de)) for x in eigenvalues])
        m=count(Pb)
    else
        m=0
    end
        
    return m
end