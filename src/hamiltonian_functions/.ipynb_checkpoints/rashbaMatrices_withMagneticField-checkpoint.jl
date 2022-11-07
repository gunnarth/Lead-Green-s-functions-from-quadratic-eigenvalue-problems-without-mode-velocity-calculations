using LinearAlgebra
#=========================H and V matrices for a graphene==========#
function rashbaMatrices_withMagneticField(Ny,kRLy,lBLy)::Tuple{Array{Complex{Float64},2},Array{Complex{Float64},2}}
    # Number of lattice points
    #
    # Ny
    #
    # Ratio of magnetic length and Ly
    #
    #lBLy
    #
    # Ratio of Rasha length and Ly
    #
    #lRLy
    #

    Nsl=2*Ny
    Hsl=zeros(Complex{Float64},Nsl,Nsl)
    Vsl=zeros(Complex{Float64},Nsl,Nsl)
    I2=Matrix{Complex{Float64}}(I, 2, 2) #julia 1.0

    sx=[0.0 1.0; 1.0 0.0]
    sy=[0.0 -im; im 0.0]
    sz=[1.0 0.0; 0.0 -1.0]

    for iy=1:Ny-1
        Hsl[2*(iy-1)+1:2*iy,2*(iy-1)+1:2*iy]=4.0*I2
        Hsl[2*(iy-1)+1:2*iy,2*iy+1:2*(iy+1)]=-1.0*I2+kRLy*im*sy/(Ny+1)
        Hsl[2*iy+1:2*(iy+1),2*(iy-1)+1:2*iy]=-1.0*I2-kRLy*im*sy/(Ny+1)
    end

    Hsl[2*(Ny-1)+1:2*Ny,2*(Ny-1)+1:2*Ny]=4.0*I2

    for iy=1:Ny
        Vsl[2*(iy-1)+1:2*iy,2*(iy-1)+1:2*iy]=(-1.0*I2-kRLy*im*sx/(Ny+1))*exp(im*lBLy^(-2)*float(iy)/float(Ny+1)^2)
    end

    return float(Ny+1)^2/pi^2*Hsl,float(Ny+1)^2/pi^2*Vsl
end