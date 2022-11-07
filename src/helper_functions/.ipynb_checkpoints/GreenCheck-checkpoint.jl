using LinearAlgebra
function GreenCheck(G,H,V,E,eta)
    try
        n=size(G,1)
        I0=Matrix{ComplexF64}(I, n, n)
        return norm(V*G*G-((E+im*eta)*I0-H)*G+V')
    catch error
       # return -1
        return missing
    end

end
