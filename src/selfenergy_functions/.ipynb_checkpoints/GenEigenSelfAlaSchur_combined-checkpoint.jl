using LinearAlgebra

function GenEigenSelfAlaSchur_combined(E,H,V,etaR=1e-14,printEig=false)
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
#        println("==== Inni ====")
#        println(((norm.(eigenvalues)-1.0)))
#        println("==============")
        #----Eigenvalues within unit circle-----------------------------------------#
        P=Array{Bool}([(!isnan(x)) && (norm(x)<(1.0-de)) for x in eigenvalues])
   
        #----Eigenvalues on the unit circle-----------------------------------------#
        Pb=Array{Bool}([((1.0-de)<=norm(x)<=(1.0+de)) for x in eigenvalues])

        la=eigenvalues[Pb]
   
        #---The inverse Power method for generalized eigenvalue problems performed for eiginvalues on the unit circle----#
        m=count(Pb)
        mask=fill(true,m) 
        qr=ones(ComplexF64,2*n)
        qld=ones(ComplexF64,2*n)
        for i in collect(1:m)
            mu=la[i]
            fill!(qr,1.0)
            fill!(qld,1.0)
            countit=0
            error=1
            while error>5e-5 && (countit<10)
                qr=(A-mu.*B)\(B*qr);
                qr=qr/norm(qr)
                qld=(AT-mu.*BT)\(BT*qld)
                qld=qld/norm(qld)
                countit+=1
                error=norm(A*qr-mu*B*qr)
            end
            dB=(transpose(qld)*B*qr)
            mup=(transpose(qld)*K*qr)/dB

            if real(mu'*mup+mu*mup')>0 
                mask[i]=false
            end  
        end  
        #----Complete the Bool mask for for which eigenvalue are within the unit circle--# 
        Pb[Pb.==true]=mask
        P=P .| Pb  

    else  

        P=Array{Bool}([(!isnan(x)) && (norm(x)<1.0) for x in eigenvalues ])
    end
    F=ordschur(F,P)
    enew=F.α./F.β

    Z11=F.Z[1:n,1:n]
    Z21=F.Z[n+1:2*n,1:n]
    if printEig
        return Z21*inv(Z11),enew
    else
        return Z21*inv(Z11)
    end
end