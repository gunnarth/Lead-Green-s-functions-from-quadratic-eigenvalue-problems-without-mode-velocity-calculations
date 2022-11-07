using LinearAlgebra
const la=LinearAlgebra

#===============GenEigenSelf============#
# V can be uninvertable
function GenEigenSelfAlaSchur_simple(E,H,V,etaR=1e-14,printEig=false)
  #------ Constructing the eigenvalue problem -------#
  n=convert(Int,size(H)[1])
  eta=etaR*im;
  O=zeros(n,n)
  I0=Matrix{ComplexF64}(la.I, n, n)
  C=(-(E+eta)*I0+H)
  K=V'
  M=V
  A=[O I0;-K -C];
  B=[I0 O; O M];

  #---Generlized Schur factoring performed-------#
  F=schur(A,B)
  e=F.alpha./F.beta
  logEta=-log10(etaR)
  P=Array{Bool}([(!isnan(x)) && (norm(x)<1.0) for x in e ])
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
