using Statistics

"""
Timer is a function that times selfenergy functions

fun: function that calculates selfenergy
e: Energy
H: Hamilton matrix
V: coupling matrix
small: eta in the greensfunction method
nmin: how many times fun is timed for an average
funName: if equal LS, return adjusted for quadratic equation (Eq. 40)

returns the minimum of the nmin times and the forward error of the
Quadtratic Matrix equation (Eq. 6) for the selfenergy function.
"""
function Timer(fun,e,H,V,small,nmin,funName)
    tmin=zeros(nmin);
    for r in 1:nmin
        try
          tmin[r]= @elapsed fun(e,H,V,small);
        catch BoundsError
            return missing,missing,missing,missing
        end
    end
    if(funName=="LS")
        G=fun(e,H,V,small)
        #The Lopes-Sanchos method only returns G, while eigenvalue methods
        #return GV^+
        G=G*V'
        forwardError=GreenCheck(G,H,V,e,0.0)
        return minimum(tmin),forwardError
    else
        G=fun(e,H,V,small)
        forwardError=GreenCheck(G,H,V,e,0.0)

        return minimum(tmin),forwardError
    end
end

