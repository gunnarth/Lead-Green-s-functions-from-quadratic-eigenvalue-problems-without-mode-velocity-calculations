function LS(E,Hsl,Vsl,eta)
#  This routine is written for a right lead with coupling
#  matrix V connecting slice n and n+1.
#

#Iteration control parameters
maxIter=5000
eps=1e-15;

Ny2=size(Hsl)[1];  # Dimension of matrices
I0=Matrix{ComplexF64}(I, Ny2, Ny2)
z=I0*(E+im*eta);  # Imaginary energy
  
#  Iteration routine
i=1;
Err=1;

H=Hsl;
HP=Hsl;
H0=Hsl;
H0P=Hsl;

V=Vsl;
VD=Vsl';
while(Err > eps && i< maxIter)

  invGz=inv(z-H);
  H0P=H0+V*invGz*VD;
  HP=H+V*invGz*VD+VD*invGz*V;
  V=V*invGz*V;
  VD=VD*invGz*VD;

  H=HP;
  H0=H0P;

  i=i+1;

  Err=maximum(maximum(abs.(V*invGz*VD)))

end
#
# Uncomment to see number of iterations needed to reach eps=10^-15
#
#println("iter=",i," norm(V)=",norm(V*invGz*VD))
return inv(z-H0)#,i,Err

end
