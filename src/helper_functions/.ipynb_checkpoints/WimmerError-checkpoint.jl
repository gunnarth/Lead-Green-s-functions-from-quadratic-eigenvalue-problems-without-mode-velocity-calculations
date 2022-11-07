using LinearAlgebra

function WimmerError(E,H,V,etaR=1e-10)
  #------ Constructing the eigenvalue problem -------#
  eta=im*etaR
  n=convert(Int,size(H)[1])
  O=zeros(ComplexF64,n,n)
  I0=Matrix{ComplexF64}(I, n, n)
  C=(-(E+eta)*I0+H)

  K=V'
  M=V
  A=[O I0;-K -C];

  B=[I0 O; O M];
  AT=transpose(A)
  BT=transpose(B)
  K=[O  O; O -im*I0];

  #---Generlized Schur factoring performed-------#
  local F
  try 
      F=schur(A,B)
  catch
      return -1
  end
  eigenvalues=F.alpha./F.beta
  de=1.0e-10
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
  wErr=zeros(m)
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
        
      A1=(transpose(qld[1:n])*V*qr[1:n])
      A2=-(transpose(qld[1:n])*C*qr[1:n])
      wErr[i]=2*imag(1.0/(2*A1-A2))
  end  
 if length(wErr) == 0
        #return -1
        return missing
 else
    return etaR*maximum(wErr)
 end
 end