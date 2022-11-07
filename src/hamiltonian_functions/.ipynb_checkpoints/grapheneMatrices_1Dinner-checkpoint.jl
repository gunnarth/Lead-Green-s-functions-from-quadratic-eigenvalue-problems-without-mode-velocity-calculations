using LinearAlgebra
const la=LinearAlgebra
#=========================H and V matrices for a graphene==========#
function grapheneMatrices_1Dinner(J)
# Number of stacked hexgons
#
#J=4
#
# Number of vertically staked atoms
#
M=2*J+1
#
# Total number of atoms in slice
#
Nsl=2*M

##
## Creating Hsl
##

# Hsl=zeros(Complex128,Nsl,Nsl) # julia 0.6
Hsl=zeros(ComplexF64,Nsl,Nsl)   # julia 1.0
t=-1.0
HJ=[0.0 t t 0.0;t 0.0 0.0 t;t 0.0 0.0 t;0.0 t t 0.0]
HV=[0.0 0.0 0.0 0.0;0.0 0.0 0.0 0.0;t 0.0 0.0 0.0;0.0 t 0.0 0.0]
#
# First fill the "cups" of the hexagons
#
for i=1:(J-1)

  Hsl[4*(i-1)+1:4*i,4*(i-1)+1:4*i]=HJ
  Hsl[4*(i-1)+1:4*i,4*i+1:4*(i+1)]=HV
  Hsl[4*i+1:4*(i+1),4*(i-1)+1:4*i]=HV'
end
Hsl[4*(J-1)+1:4*J,4*(J-1)+1:4*J]=HJ

#
#  Finalize the top
#
Hsl[Nsl-1:Nsl,Nsl-1:Nsl]=[0.0 t;t 0.0]
Hsl[Nsl-3,Nsl-1]=t
Hsl[Nsl-2,Nsl  ]=t
Hsl[Nsl-1,Nsl-3]=t
Hsl[Nsl  ,Nsl-2]=t

##
## Creating Vsl
##
# Vsl=zeros(Complex128,Nsl,Nsl)  # julia 0.6
Vsl=zeros(ComplexF64,Nsl,Nsl)    # julia 1.0
for i=1:J
  Vsl[4*i,4*i-1]=t
end
    return Hsl[1:2:end,1:2:end],Vsl[1:2:end,1:2:end]
end
