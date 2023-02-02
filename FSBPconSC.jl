using LinearAlgebra


# Function finds a mimetic boundary operator by calling
# FindSurfCubSVDONB
# with different values for the basis cutoff, i.e. the tolerance under which 
# a singular value under which respective functions are declared linearly dependend
# in the inner product defined throught the collocation points
function FindSurfCubIter(Op, N)
        minind = 3                                                      # start at 10^-3
        maxind = 16                                                     #  end trying at 10^-16
        
        errar = zeros(maxind)                                           # save the errors in an array
        errar .= Inf                                                    # and initialize them as infinite
        bestMx, bestMy = 0, 0                                           # memorize the best boundary operators.
        for k = minind:maxind                                           # try out the respective tolerances
                println("Trying out ", 10.0^(-k))                       
                Mx, My, err = FindSurfCubSVDONB(Op, N, 10.0^(-k))       # build boundary operators with this tolerance
                
                if err < minimum(errar)                                 # if their error is smaller than the previous best
                        bestMx, bestMy = Mx, My                         # they are saved as new best candidates
                end
                errar[k] = err                                          # and their error memorized
        end

        k = argmin(errar)                                               
        println("Returning for ", 10.0^(-k), " with total error ", errar[k])
        return bestMx, bestMy
        
end


# Returns a set of surface operators and their respective errors
# when a cutoff tolerance of tol1 is used to sort out basis functions 
# that are not relevant
# General overview of the algorithm:
#       - Build a list of all possible products of basis functions, as we would like to
#               design a quadrature rule exact for all products of basis functions.
#       - Assemble matrices Ax, Ay and vectors bx, by defining exactness conditions
#       - Decompose both matrices into USV = Ax and USV = Ay
#       - Use U and S to calculate a new basis
#       - Recalculate Ax, Ay and bx, by with respect to the new basis
#       - Use the POCS algorithm to find a vector of quadrature weights
#       - Calculate the errors (||Aw - b||^2 + || Aw - b ||^2)^1/2

function FindSurfCubSVDONB(Op, N=1000, tol1=1.0E-8, tol2=1.0E-3, debug = 0)
        nfunc = length(Op.funcs)                                                        # number of basis functions
        M = zeros(Op.No)                                                                # weights for surface quadrature 
              
                                                                                        # building the square function list
        funcsq = Array{Any}(undef, nfunc*nfunc)                                         # A list of the squaref functions        
        for k=1:nfunc                                                                   # iterating over all combinatins of
                for l=1:nfunc                                                           # basis functions.
                        funcsq[k + nfunc*(l-1)] = x->Op.funcs[k](x)*Op.funcs[l](x)
                end
        end
        
        nfuncsq = nfunc*nfunc                                                           # the number of new basis functions

        Axprim = zeros(nfuncsq, Op.No)                                                  # Ax and Ay map a quadrature rule to its results
        Ayprim = zeros(nfuncsq, Op.No)                                                  # when applied to all functions and their products

        bxprim, byprim = zeros(nfuncsq), zeros(nfuncsq)                                 # holds boundary integrals of basis functions
         
                                                                                        # initialisation of Ax, Ay, bx, by      
        for k=1:nfuncsq
                for l=1:Op.No
                        pt, normal = Op.cpoints[l, :], Op.normals[l, :]                 # using the saved node and normal values
                        Axprim[k, l] = funcsq[k](pt)*normal[1]                          # every row of Ax multiplies the weight, if multiplied
                        Ayprim[k, l] = funcsq[k](pt)*normal[2]                          # from the right, with Ax, with the normal and the function value
                end                                                                     # and vice versa for Ay
                print(".")
                integ1(x) = funcsq[k](Op.Omega.param(x))*Op.Omega.nparam(x)[1]          # integrand for the surface integral using domain and normal parametrisation 
                bxprim[k], err = quadgk(integ1, Op.Omega.discs..., maxevals=10^3)       # surface integral, calculated using numerical quadrature
                integ2(x) = funcsq[k](Op.Omega.param(x))*Op.Omega.nparam(x)[2]          
                byprim[k], err = quadgk(integ2, Op.Omega.discs..., maxevals=10^3)
        end
        println("")

        if debug == -2                                                                  # The original Ax and Ay can be returned for debug reasons
               return Ax, Ay
        end

        sAx, sAy = svd(Axprim, full=true), svd(Ayprim, full=true)                       # Calculation of the SVDs of the original Ax and Ay matrices.
                                                                                        
        println("Detected ", sum(sAx.S .> tol1), " nondegenerate x basis functions")    # Every singular value above tol1 designates a relevant basis function
        println("Detected ", sum(sAy.S .> tol1), " nondegenerate y basis functions")    

                                                                                        # Calculate scaling/normalization factors for the new basis functions
        sfx, sfy = zeros(nfuncsq), zeros(nfuncsq)                                       # these are saved in sfx and sfy
        for k = 1:min(nfuncsq, Op.No)                                                   # If the function is discarded the corresponding scaling value is 
                sfx[k] = sAx.S[k] > tol1 ? 1.0/(sAx.S[k]) : 0.0                         # set to zero, otherwise to the inverse singular value to enforce
                sfy[k] = sAy.S[k] > tol1 ? 1.0/(sAy.S[k]) : 0.0                         # a normed basis vector in the half-norm induces by the two-norm
        end                                                                             # of the vector formed by the basis function evaluated at the nodes.

        funcsqx, funcsqy = Array{Any}(undef, nfuncsq), Array{Any}(undef, nfuncsq)       # These function lists will collect the new base functions. 
        for k=1:nfuncsq                                                                 # every new base functions is the result of the linear combination. 
                funcsqx[k] = function (x)                                               # The coefficients are given by the columns of U
                        val = 0.0                                                       # and scaling with sfx results in a discard or scaling to one
                        for l=1:nfuncsq                                                 # Warning:  We are multiplying with the adjoint - therefore swap the indices
                                val = val + sfx[k]*sAx.U[l, k]*funcsq[l](x)
                        end
                        return val
                end
                funcsqy[k] = function (x)                                               # vice versa for the y basis functions with the corresponding U stemming
                        val = 0.0                                                       # from Ay
                        for l=1:nfuncsq # See above - indices are swaped
                                val = val + sfy[k]*sAy.U[l, k]*funcsq[l](x)
                        end
                        return val
                end
        end

                                                                                        # Now: Build new Ax and Ay matrices encoding the exactness conditions
        Ax = zeros(nfuncsq, Op.No)                                                      # This time with the rhs bx, by and matrix elements of Ax and Ay
        Ay = zeros(nfuncsq, Op.No)                                                      # encoded with respect to the new bases in saved in funcsqx, funcsqy
        bx, by = zeros(nfuncsq), zeros(nfuncsq)
        for k=1:nfuncsq
                for l=1:Op.No
                        pt, normal = Op.cpoints[l, :], Op.normals[l, :]                 # Note: same as above, but funcsqx and funcsqy are evaluated. 
                        Ax[k, l] = funcsqx[k](pt)*normal[1]
                        Ay[k, l] = funcsqy[k](pt)*normal[2]
                end
                print(".")
                integ1(x) = funcsqx[k](Op.Omega.param(x))*Op.Omega.nparam(x)[1]
                bx[k], err = quadgk(integ1, Op.Omega.discs..., maxevals=10^3)
                integ2(x) = funcsqy[k](Op.Omega.param(x))*Op.Omega.nparam(x)[2]  
                by[k], err = quadgk(integ2, Op.Omega.discs..., maxevals=10^3)
        end
        println("")

                                                                                        # As sfx and sfy can become large when a function vanishes on the collocation points
        if max(maximum(abs.(bx)), maximum(abs.(by))) > 10.0                             # the integrals in bx and by can also become large. Therefore, a warning is in place.
                println("Warning: High quad value")                                     # A firing of this warning is a symptom of (far) to few collocation points on the boundary.
        end
        Axi = pinv(Ax, atol=tol2)                                                       # recalculate the penrose inverses
        Ayi = pinv(Ay, atol=tol2)


        for k=1:N                                                                       # This loop implements the POCS algorithm
                rx = bx - Ax*M                                                          # These lines project onto the two exactness conditions
                ry = by - Ay*M                                                          # for the quadrature
                corrx = Axi*rx
                corry = Ayi*ry
                M = M + 0.5.*(corrx + corry)                    
                M = max.(0.0, M)                                                        # This line projects the weights onto the positive orthant.
        end
                                                                                        # Errors are calculates using the original exactness conditions!
                                                                                        # Otherwise conditions could be removed in error and the results for
                                                                                        # different tol1 would be not comparable.
        println("Errors of x Boundary norm matrix exactness: ", norm(bxprim - Axprim * M)) 
        println("Errors of y Boundary norm matrix exactness: ", norm(byprim - Ayprim * M))

        toterr = sqrt(norm(bxprim - Axprim * M)^2 + norm(byprim - Ayprim * M)^2)
        Mx = M .* Op.normals[:, 1]                                                      # The usable quadrature rule for surve quadratures are the result
        My = M .* Op.normals[:, 2]                                                      # of the product between the quadrature weights and the normal components.

        if debug==-1                                                                    # For debug reasons one can return the exactness matrices
                return Ax, Ay                                                           # otherwise the resulting quadratures and the error is returned.
        else
                return Mx, My, toterr
        end
end

