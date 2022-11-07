using LinearAlgebra
#===============GenEigenSelf============#
# V can be uninvertable
function GenEigenSelf_new(E,H,V,etaR=1e-14)
  #------ Constructing the eigenvalue problem -------#
  n=size(H)[1]
  etaR=1e-10
  #eta=etaR*im;
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
    
  #-----Eigenvalues an eigenvectors calculated------#
  e,b=eigen(A,B);
   
    #----Eigenvalues within unit circle-----------------------------------------#
  P=Array{Bool}([(!isnan(x)) && (norm(x)<(1.0-etaR)) for x in e])

  #----Eigenvalues on the unit circle-----------------------------------------#
  Pb=Array{Bool}([((1.0-etaR)<=norm(x)<=(1.0+etaR)) for x in e])

  la=e[Pb]


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

  sb=b[:,P][1:n,1:n]
  se=e[P]
   

  #--eigenvectors normalized in the new nXn subspace----#
  for i in 1:n
       sb[:,i]=sb[:,i]/norm(sb[:,i])
  end


  #---SelfEnergy constructed and returned-------------#
  return  sb*Matrix(Diagonal(se))*inv(sb)#,eold
 end