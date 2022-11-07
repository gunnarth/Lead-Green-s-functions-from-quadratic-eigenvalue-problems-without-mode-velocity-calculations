using LinearAlgebra
#===============GenEigenSelf============#
# V can be uninvertable
function GenEigenSelf(E,H,V,etaR=1e-14)
  #------ Constructing the eigenvalue problem -------#
  n=size(H)[1]
  eta=etaR*im;
  O=zeros(ComplexF64,n,n)
  I0=Matrix{ComplexF64}(I,n, n)
  C=(-(E+eta)*I0+H)
  K=V'
  M=V
  A=[O I0;-K -C];
  B=[I0 O; O M];


  #-----Eigenvalues an eigenvectors calculated------#
  e,b=eigen(A,B);
  eold=e
  #--Nan eigenvalues and negative imaginary part eigenvalues and their vectors excluded--#
  b=b[:,.!(isnan.(e))]
  e=e[.!(isnan.(e))]
  b=b[:,abs.(e).<1.0]
  e=e[abs.(e).<1.0]

  #--einvalues and their vectors sorted by size, only subspace of n smallest used-----#
  p=sortperm(abs.(e))
  sb=b[:,p][1:n,1:n]
  se=e[p][1:n]

  #--eigenvectors normalized in the new nXn subspace----#
  for i in 1:n
       sb[:,i]=sb[:,i]/norm(sb[:,i])
  end


  #---SelfEnergy constructed and returned-------------#
  return  sb*Matrix(Diagonal(se))*inv(sb)#,eold
 end
