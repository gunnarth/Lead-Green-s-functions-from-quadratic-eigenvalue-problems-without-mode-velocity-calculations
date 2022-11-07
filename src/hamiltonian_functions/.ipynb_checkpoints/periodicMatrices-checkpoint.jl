using LinearAlgebra
const la=LinearAlgebra
#=========================H and V matrices for a graphene==========#
function periodicMatrices(Nx,Ny,V1c,V2c,V3c,lBLy=1e5)
#
# Nx: Number of point in each (periodic) cell
#
# Ny: Number of points in transverse direction 

# Periodic potential V(x)=V1c*cos(1*2 \pi x /L)+V2c*cos(2*2 \pi x /L)+V3c*cos(3*2 \pi x /L)
#
# The matrix Hsl contains V(x) evaluated at x_1=0, ..., x_N=L*N/(L+1)
#
# lBLy: The ratio of magnetic length lB and Ly [Default value B \approx 0]
#
#  The energy is scaled in units of \hbar^2Q^2/2m where Q=2 \pi/L_x 
#

# The dimension of the matrices
Nsl=Nx*Ny
INy=Matrix{ComplexF64}(la.I, Ny, Ny) #julia 1.0
#
#  Calculate the periodic potential
#
xVec=collect(0:Nx-1)/float(Nx)
VP=V1c*cos.(2.0*pi*xVec)+V2c*cos.(4.0*pi*xVec)+V2c*cos.(6.0*pi*xVec)  # Measured relative to \hbar^2 Q^2/2m where Q=\pi/L

invaQ2=float(Nx)^2/(2.0*pi)^2

##
##  Start with small slice of size Ny*Ny
##
##  Works for Ny=0, which gives the correct 1D version.
##
Hsl=zeros(ComplexF64,Ny,Ny)
#Hsl=invaQ2*((2.0+((Ny>1) ? 2.0 : 0.0))*diagm(vec(ones(Ny,1)),0)-1.0*diagm(vec(ones(Ny-1,1)),+1)-1.0*diagm(vec(ones(Ny-1,1)),-1)) #julia 0.6
Hsl=invaQ2*(
            (2.0+((Ny>1) ? 2.0 : 0.0))*
             la.diagm(-1=>-ones(Ny-1),
                       0=> ones(Ny  ),
                       1=>-ones(Ny-1)
                     )      
            )
Vsl=zeros(ComplexF64,Ny,Ny)
for iy=1:Ny
    Vsl[iy,iy]=invaQ2*(-1.0*exp(im*lBLy^(-2)*float(iy)/float(Ny+1)^2))
end
##
## Creating the super-slice matrix HHsl
##
HHsl=zeros(ComplexF64,Nsl,Nsl)
VVsl=zeros(ComplexF64,Nsl,Nsl)
for ix=1:Nx-1
    #HHsl[Ny*(ix-1)+1:Ny*ix,Ny*(ix-1)+1:Ny*ix]=Hsl+VP[ix]*eye(Ny,Ny) #julia 0.6
    HHsl[Ny*(ix-1)+1:Ny*ix,Ny*(ix-1)+1:Ny*ix]=Hsl+VP[ix]*INy
    HHsl[Ny*(ix-1)+1:Ny*ix,Ny*ix+1:Ny*(ix+1)]=Vsl
    HHsl[Ny*ix+1:Ny*(ix+1),Ny*(ix-1)+1:Ny*ix]=Vsl'
end
#HHsl[Ny*(Nx-1)+1:Ny*Nx,Ny*(Nx-1)+1:Ny*Nx]=Hsl+VP[Nx]*eye(Ny,Ny) #julia 0.6
HHsl[Ny*(Nx-1)+1:Ny*Nx,Ny*(Nx-1)+1:Ny*Nx]=Hsl+VP[Nx]*INy
VVsl[Ny*(Nx-1)+1:Ny*Nx,1:Ny]=Vsl

return (2*float(Ny+1)/float(Nx))^2*HHsl,(2*float(Ny+1)/float(Nx))^2*VVsl

end
