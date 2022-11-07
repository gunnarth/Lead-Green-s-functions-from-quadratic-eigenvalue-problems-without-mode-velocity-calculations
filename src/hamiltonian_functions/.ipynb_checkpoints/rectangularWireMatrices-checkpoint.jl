using LinearAlgebra
const la=LinearAlgebra
#=========================H and V matrices for a square ==========#
function rectangularWireMatrices(Nx,Ny,lBLx=1e5)
#
#  NEW VERSION WITH TWO LAYERS TO GET CORRECT CORNER LOCALIZATION
#
# Nx: Number of points on the top surfaces.
#
# Ny: Number of points on side surfaces.
#
# lBLx: The ratio of magnetic length lB and Lx [Default value B \approx 0]
#
#  The energy is scaled in units of \hbar^2Q^2/2m where Q=2 \pi/L_x 
#

#
# Creating an effective model with two outer layers
#
# Outer layer, layer 1
#
N1=2*Nx+2*(Ny-2)  # Corners belong to Nx
H1=zeros(ComplexF32,N1,N1)
H1=diagm(0=>vec(fill(6.0,(N1,1))))+diagm(+1=>vec(fill(-1.0,(N1-1,1))))+diagm(-1=>vec(fill(-1.0,(N1-1,1))))
H1[1,N1]=-1
H1[N1,1]=-1

#
# Inner layer, layer 2
#
N2=2*(Nx-2)+2*((Ny-2)-2)  # Corners belong to (Nx-2)
H2=zeros(ComplexF32,N2,N2)
H2=diagm(0=>vec(fill(6.0,(N2,1))))+diagm(+1=>vec(fill(-1.0,(N2-1,1))))+diagm(-1=>vec(fill(-1.0,(N2-1,1))))
H2[1,N2]=-1
H2[N2,1]=-1

#
#  Couling between layers 1 and 2
#
V21=zeros(ComplexF32,N2,N1)
for is=1:(Nx-2)
   V21[is,is+1]=-1.0 
end
V21[1,N1]=-1.0
V21[Nx-2,Nx+1]=-1.0

for is=(Nx-2)+1:(Nx-2)+(Ny-2)-2
   V21[is,is+3]=-1.0  
end

for is=(Nx-2)+(Ny-2)-1:2*(Nx-2)+(Ny-2)-2
   V21[is,is+5]=-1.0  
end
V21[(Nx-2)+(Ny-2)-1,Nx+Ny-2]=-1
V21[2*(Nx-2)+(Ny-2)-2,2*Nx+(Ny-2)+1]=-1

for is=2*(Nx-2)+(Ny-2)-1:2*(Nx-2)+2*((Ny-2)-2)
   V21[is,is+7]=-1.0  
end

#
#  Create the effective slice model combining both layers
#
Nsl=N1+N2
Hsl=zeros(ComplexF32,Nsl,Nsl)
Hsl[1:N1,1:N1]=H1
Hsl[N1+1:N1+N2,N1+1:N1+N2]=H2
Hsl[N1+1:N1+N2,1:N1]=V21
Hsl[1:N1,N1+1:N1+N2]=V21'
    
#
#  Interslice coupling Vsl
#
Vsl1=zeros(ComplexF32,N1,N1)
u1Vec=vec(range(1.0/(Nx+1),step=1.0/(Nx+1),length=Nx))
Vsl2=zeros(ComplexF32,N2,N2) 
u2Vec=u1Vec[2:Nx-1]
Vsl=Matrix{ComplexF64}(la.I, Nsl, Nsl) #julia 1.0
    
for ix=1:Nx
    Vsl1[ix,ix]=-exp(im*(lBLx)^(-2)*(u1Vec[ix]-0.5)/(Nx-1))
end
for ix=1:(Nx-2)
    Vsl2[ix,ix]=-exp(im*(lBLx)^(-2)*(u2Vec[ix]-0.5)/(Nx-1))
end
for iy=1:Ny-2
    Vsl1[iy+Nx,iy+Nx]=-exp(im*(lBLx)^(-2)*(u1Vec[Nx]-0.5)/(Nx-1))  
end
for iy=1:(Ny-2)-2
    Vsl2[iy+(Nx-2),iy+(Nx-2)]=-exp(im*(lBLx)^(-2)*(u2Vec[Nx-2]-0.5)/(Nx-1))  
end
for ix=1:Nx
    Vsl1[ix+Nx+(Ny-2),ix+Nx+(Ny-2)]=-exp(im*(lBLx)^(-2)*(u1Vec[Nx+1-ix]-0.5)/(Nx-1))
end
for ix=1:(Nx-2)
    Vsl2[ix+(Nx-2)+(Ny-2)-2,ix+(Nx-2)+(Ny-2)-2]=-exp(im*(lBLx)^(-2)*(u2Vec[(Nx-2)+1-ix]-0.5)/(Nx-1))
end
for iy=1:Ny-2
    Vsl1[iy+2*Nx+(Ny-2),iy+2*Nx+(Ny-2)]=-exp(im*(lBLx)^(-2)*(u1Vec[1]-0.5)/(Nx-1))  
end
for iy=1:(Ny-2)-2
    Vsl2[iy+2*(Nx-2)+(Ny-2)-2,iy+2*(Nx-2)+(Ny-2)-2]=-exp(im*(lBLx)^(-2)*(u2Vec[1]-0.5)/(Nx-1))  
end
Vsl[1:N1,1:N1]=Vsl1
Vsl[N1+1:N1+N2,N1+1:N1+N2]=Vsl2
    
#return Hsl*(Nx+1)^2,Vsl*(Nx+1)^2
return Hsl,Vsl

end
