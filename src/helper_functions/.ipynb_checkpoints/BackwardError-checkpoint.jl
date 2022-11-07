using LinearAlgebra

function BackwardError(E,H,V)
    try
        #------ Constructing the eigenvalue problem -------#
        n=convert(Int,size(H)[1])
        O=zeros(ComplexF64,n,n)
        I0=Matrix{ComplexF64}(I, n, n) 
        C=(-(E+0.0)*I0+H)

        K=V'
        M=V
        A=[O I0;-K -C];

        B=[I0 O; O M];
        AT=transpose(A)
        BT=transpose(B)
        K=[O  O; O -im*I0];

        #---Generlized Schur factoring performed-------#
        F=schur(A,B)

        eigenvalues=F.alpha./F.beta

        de=1.0e-6
        de=1e3
        de=0.9
        if any(abs.(norm.(eigenvalues).-1.0).<de)
            #----Eigenvalues on the unit circle-----------------------------------------#
            Pb=Array{Bool}([((1.0-de)<=norm(x)<=(1.0+de)) for x in eigenvalues])

            la=eigenvalues[Pb]
            #---The inverse Power method for generalized eigenvalue problems performed for eiginvalues on the unit circle----#
            m=count(Pb)
            be=zeros(Float64,m) 
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
                be[i]=eps(1.0)/sqrt.(abs.(transpose(qld)*A*qr)^2+abs.(transpose(qld)*B*qr)^2)
            end  
        else
                #return -1
                return missing
        end
        return maximum(be)
    catch error
        println(error)
        #return -1
        return missing
    end
end