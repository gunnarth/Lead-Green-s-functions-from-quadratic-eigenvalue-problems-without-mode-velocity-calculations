using LinearAlgebra

function GenEigenSelfAlaSchur(E,H,V,etaR=1e-10,printEig=false)
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
  de=1.0e-14
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

  #----Correct eigenvalues (within unit circle) calculated and returned-------------#
  F=ordschur(F,P)
  enew=F.alpha./F.beta
  Z11=F.Z[1:n,1:n]
  Z21=F.Z[n+1:2*n,1:n]

  if printEig
    return Z21*inv(Z11),enew
  else
    return Z21*inv(Z11)
  end
 end