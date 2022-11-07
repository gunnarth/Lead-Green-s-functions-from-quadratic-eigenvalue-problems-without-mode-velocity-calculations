using LinearAlgebra
const la=LinearAlgebra
#=========================H and V matrices for a core/shell ==========#
function snakingMatrices(M,lcR,lER)
# Highest value of angular momentum 
#
# M -> Dimension of matrices Nsl=2M+1
#
# Ratio of magnetic length and R
#
#lcR
#
# Ratio of electric field confinement length and R
#
#lER  # Note that lER-> \inft : E->0  and E becomes relevant when lER \approx lcR
#
    #
    # aR=a/R  # ratio of lattice constant a and radius R
    #  
    aR=1.0/(float(M)*sqrt(1.0+lcR^(-4)))
    ##
    ## Matrices scaled in E0
    ##
    Nsl=2*M+1
    H=zeros(ComplexF64,Nsl,Nsl)
    H=diagm(0=>range(-M,step=1,length=2*M+1).^2)+lcR^(-4)*(0.5*diagm(0=>vec(fill(1.0,Nsl,1)))-0.25*diagm(-2=>vec(fill(1.0,Nsl-2,1)))-0.25*diagm(2=>vec(fill(1.0,Nsl-2,1)))) +0.5*lER^(-3)*(diagm(-1=>vec(fill(1.0,Nsl-1,1)))+diagm(1=>vec(fill(1.0,Nsl-1,1))))   
    V=zeros(ComplexF64,Nsl,Nsl)
    V=2.0*lcR^(-2)*(-0.5*im)*(diagm(-1=>vec(fill(1.0,Nsl-1,1)))-diagm(1=>vec(fill(1.0,Nsl-1,1))))

    Hsl=zeros(ComplexF64,Nsl,Nsl)
    Hsl=(H+2.0*aR^(-2)*Matrix{ComplexF64}(1.0I, Nsl, Nsl))
    Vsl=zeros(ComplexF64,Nsl,Nsl)
    Vsl=(-1.0*aR^(-2)*Matrix{ComplexF64}(1.0I, Nsl, Nsl)-im*0.5*aR^(-1)*V)
    
    return Hsl,Vsl
    
end