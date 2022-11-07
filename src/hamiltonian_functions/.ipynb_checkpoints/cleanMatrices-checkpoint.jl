using LinearAlgebra
#=========================H and V matrices for a graphene==========#
function cleanMatrices(Ny)::Tuple{Array{Complex{Float64},2},Array{Complex{Float64},2}}
# Number of lattice points
#
# Ny
#
Nsl=2*Ny

Hsl=zeros(Complex{Float64},Nsl,Nsl)
Vsl=zeros(Complex{Float64},Nsl,Nsl)

for i=1:Ny-1
    Hsl[2*(i-1)+1:2*i,2*(i-1)+1:2*i]=[4.0 0.0 ; 0.0 4.0]
    Hsl[2*(i-1)+1:2*i,2*(i)+1:2*(i+1)]=[-1.0 0.0 ; 0.0 -1.0]   
    Hsl[2*(i)+1:2*(i+1),2*(i-1)+1:2*i]=[-1.0 0.0 ; 0.0 -1.0]   
        
    Vsl[2*(i-1)+1:2*i,2*(i-1)+1:2*i]=[-1.0 0.0 ; 0.0 -1.0]
    
end
i=Ny
Hsl[2*(i-1)+1:2*i,2*(i-1)+1:2*i]=[4.0 0.0 ; 0.0 4.0]
Vsl[2*(i-1)+1:2*i,2*(i-1)+1:2*i]=[-1.0 0.0 ; 0.0 -1.0]
    
return float(Ny+1)^2/pi^2*Hsl,float(Ny+1)^2/pi^2*Vsl

end