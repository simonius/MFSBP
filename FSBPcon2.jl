using LinearAlgebra                                                                             # Used for Matrix inversion and penrose inverses, etc.

include("structs.jl")                                                                           # include definitions of structs used globaly
include("FSBPconSC.jl")                                                                          # include stabilized mimetic surface quadrature construction


# Uses the POCS algorithm to find a surface quadrature for the operator described in P
# N iterations of the pocs algorithm are done.
# Consult also "FSBPconSC.jl" for a more mature version.
function FindSurfCubQR(P, N)
        nfunc = length(P.funcs)
        Mx = zeros(P.No)
        My = zeros(P.No)

        Ax = zeros(nfunc*nfunc, P.No) #Ax and Ay map a quadrature rule to its results
        Ay = zeros(nfunc*nfunc, P.No) #when applied to all functions and their products

        bx, by = zeros(nfunc*nfunc), zeros(nfunc*nfunc) # holds boundary integral of basis functions
	for k=1:nfunc
        	for l=1:nfunc
                        for m=1:P.No               
                                pt, normal = P.cpoints[m, :], P.normals[m, :]
                                Ax[k + nfunc*(l-1), m] = P.funcs[k](pt)*P.funcs[l](pt)*normal[1]
                                Ay[k + nfunc*(l-1), m] = P.funcs[k](pt)*P.funcs[l](pt)*normal[2]
                        end
                        bx[k + nfunc*(l-1)] = P.BxEl(k, l)
                        by[k + nfunc*(l-1)] = P.ByEl(k, l)
                end
        end
        
        Axi = pinv(Ax, atol=1.0E-10)			# The penrose inverses are needed to calculate the orthogonal projection onto the set of
        Ayi = pinv(Ay, atol=1.0E-10)			# exact quadratures

        for k=1:N					# This loop implements the POCS algorithm to find a solution to the positive quadrature problem.
                rx = bx - Ax*Mx				# rx is the residual of the quadrature exactness conditions
                ry = by - Ay*My
                corrx = Axi*rx				# And the Penrose inverse of this residual is the smalles distance possible
                corry = Ayi*ry				# to translate into the set of exact quadratures
                Mx = Mx + corrx
                My = My + corry
                Mx = max.(0, Mx)			# The projection onto the positive quadratures is given by cuting of negative values.
                My = max.(0, My)
        end
        println("Errors of Boundary norm matrix exactness: " ,norm(bx - Ax*Mx[1:P.No]) + norm(by - Ay*My[1:P.No]))
        println("Minimum surface quadrature entry: ", min(minimum(Mx), minimum(My)))
        Mx = Mx .* P.normals[:, 1]			# As the vectors Mx and My only holds weights without the normals the normals
        My = My .* P.normals[:, 2]			# are multiplied with the positive weights.
        return Mx, My
end

# Returns a Vandermonde matrix for the basis and collocation points in P
function MakeVander(P)
        nfunc = length(P.funcs)
        V = zeros(P.No+P.Ni, nfunc)
        for k=1:(P.No+P.Ni)
                for l=1:nfunc
                        V[k, l] = P.funcs[l](P.cpoints[k, :])
                end
        end
        return V
end

# Returns the Vandermonde matrix for the derivatives in x direction
# of the basis and collocation points in P.
function MakeDxVander(P)
        nfunc = length(P.funcs)
        dV = zeros(P.No+P.Ni, nfunc)
        for k=1:(P.No + P.Ni)
                for l=1:nfunc
                        dV[k, l] = P.dxfuncs[l](P.cpoints[k, :])
                end
        end
        return dV
end

# See function above, vice versa for the y direction.
function MakeDyVander(P)
        nfunc = length(P.funcs)
        dV = zeros(P.No + P.Ni, nfunc)
        for k=1:(P.No + P.Ni)
                for l=1:nfunc
                        dV[k, l] = P.dyfuncs[l](P.cpoints[k, :])
                end
        end
        return dV
end


# Returns a volume quadrature, found using the POCS algorithm.
# P is the operator definition and N the maximum number of POCS iterations
# dfac factor how much the minimum quadrature value is allowed to lie 
# below equidistributed weights.
function FindQuadQR(P, N, dfac = 0.1)
        nfunc = length(P.funcs)
        M = zeros(P.No+P.Ni)
                                                                                                # A maps a quadrature rule onto its values
                                                                                                # for all combinations of $dx(fg)$ and $dy(fg)$
        A = zeros(nfunc*(nfunc + nfunc), P.No + P.Ni)                                           # b saves the integrals $\int dx(fg) dV$ 
        b = zeros(nfunc*(nfunc + nfunc))                                                        # even rows are for x direction and uneven rows for
        for k=1:nfunc                                                                           # y direction derivatives
                for l=1:nfunc
                        for m=1:(P.No + P.Ni)
                                pt = P.cpoints[m, :]                                            # pt is the node coordinate
                                A[2*(k + nfunc*(l-1)), m] = P.funcs[k](pt)*P.dxfuncs[l](pt) + P.dxfuncs[k](pt)*P.funcs[l](pt)
                                A[2*(k + nfunc*(l-1)) - 1, m] = P.funcs[k](pt)*P.dyfuncs[l](pt) + P.dyfuncs[k](pt)*P.funcs[l](pt)
                        end
                        b[2*(k + nfunc*(l-1))] = P.SxEl(k, l)  
                        b[2*(k + nfunc*(l-1)) - 1] = P.SyEl(k, l)
                end
        end
                                                                                                # A will in have linearly dependend rows, and bad condition numbers,
        Ai = pinv(A)                                                                            # because the fg product is commutative and even disjoint $f_1g_1$ and $f_2g_2$ 
        MinVal = dfac/(P.No + P.Ni)                                                             # can be nearly colinear. Therefore, a Penrose inverse is used
        for k=1:N                                                                               # Every loop run does one complete run through all projections in the POCS algorithm.
                r = b - A*M                                                                     # Projection onto the 
                corr = Ai*r                                                                     # Exactness condition
                M = M + corr                                                                    # using the Penrose inverse 
                M = max.(MinVal, M)                                                             # Projection onto the set of Vectors with minimum value MinVal
        end
        println("Quadrature Error: ", norm(b-A*M, Inf))
        return M, A, b                                                                          # Return the quadrature weights, the exactness condition matrix A and the RHS b.
end


# Constructs a MFSBP Operator with mimetic boundary matrices Bx and By
# that is based on a single quadrature rule, multiplied with the norma directions.
# Using a maximum number of N POCS iterations. 
function MakeFSBPmimeticB(Op, N)
        nfunc = length(Op.funcs)
        V = MakeVander(Op)                                                                      # Initialize V, Vx and Vy as the Vandermonde matrices for
        Vx = MakeDxVander(Op)                                                                   # the given nodes and functions or derivatives of functions
        Vy = MakeDyVander(Op)                                                                   # 
        Bx, By = zeros(Op.No + Op.Ni, Op.No+ Op.Ni), zeros(Op.No + Op.Ni, Op.No+ Op.Ni)         # Define the matrices for the mimetic boundary operators
        diagBx, diagBy = FindSurfCubIter(Op, 10000)                                             # Calculate the surface quadratures for the diagonals
        Bx[1:Op.No, 1:Op.No] = diagm(diagBx)                                                    # and set them
        By[1:Op.No, 1:Op.No] = diagm(diagBy)
        M, A, b = FindQuadQR(Op, 30000)                                                         # Calculate quadrature rule for the mass matrix.
        M = diagm(M)                                                                            # and write it onto the diagonal.
        Rx = M*Vx - Bx*V/2                                                                      # precalculate values used for the orth. projection onto
        Ry = M*Vy - By*V/2                                                                      # the set of exact differentiation matrices.
        SxM = zeros(Op.No + Op.Ni, Op.No + Op.Ni)                                               # Initialize space for the Skew symmetric parts of the
        SyM = zeros(Op.No + Op.Ni, Op.No + Op.Ni)                                               # differentiation matrices.

        Vi = pinv(V)                                                                            # precalculate the penrose inverse of the Vandermone matrix.
        for k=1:N                                                                               # Loop for the POCS algorithm.
                rx = Rx - SxM*V                                                                 # Projection onto the exactness conditions
                ry = Ry - SyM*V                                                                 # for the Skew symmetric part

                corrx = rx*Vi
                corry = ry*Vi

                SxM = SxM + corrx                                                               # Addition of the correction calculated above.
                SyM = SyM + corry

                SxM = (SxM - transpose(SxM))/2                                                  # Projection of the Skewsymmetric part onto the skewsymmetric matrices
                SyM = (SyM - transpose(SyM))/2                                                  # to enforce the Skewsymmetry.
        end
        println("S Error: ", norm(Rx - SxM*V) + norm(Ry - SyM*V))
        
        Mi = pinv(M)                                                                            # Invert the mass matrix for later use.
        Dx, Dy =  Mi*(SxM + Bx/2), Mi*(SyM + By/2)                                              # Calculate the stiffness matrices by assembling the symmetric
        Op.P = M                                                                                # part, defined by the boundary operator, and the skewsymmetric part.
        Op.Pi = Mi                                                                              # Save the results in the struct.
        Op.Dx, Op.Dy = Dx, Dy
        Op.Bx, Op.By = Bx, By

	return Op.P, Op.Dx, Op.Dy, Op.Bx, Op.By, sqrt(norm(Vx - Dx*V)^2 + norm(Vy - Dy*V)^2)
end


# Constructs a MFSBP Operator
# Using a maximum number of N POCS iterations.
# The boundary matrices Bx and By are constructed using two separate quadrature rules.
# Therefore not advised for multiblock constructions
function MakeFSBPbB(Op, N)
	nfunc = length(Op.funcs)
        V = MakeVander(Op)
        Vx = MakeDxVander(Op)
        Vy = MakeDyVander(Op)
        Bx, By = zeros(Op.No + Op.Ni, Op.No+ Op.Ni), zeros(Op.No + Op.Ni, Op.No+ Op.Ni)
        diagBx, diagBy = FindSurfCubQR(Op, 30000)
        Bx[1:Op.No, 1:Op.No] = diagm(diagBx)
        By[1:Op.No, 1:Op.No] = diagm(diagBy)
        M, A, b = FindQuadQR(Op, 30000)
        M = diagm(M)
        Rx = M*Vx - Bx*V/2
        Ry = M*Vy - By*V/2
        SxM = zeros(Op.No + Op.Ni, Op.No + Op.Ni)
        SyM = zeros(Op.No + Op.Ni, Op.No + Op.Ni)

        Vi = pinv(V)
        for k=1:N
                SxMold = SxM
                SyMold = SyM

                rx = Rx - SxM*V
                ry = Ry - SyM*V

                corrx = rx*Vi
                corry = ry*Vi

                SxM = SxM + corrx
                SyM = SyM + corry

                SxM = (SxM - transpose(SxM))/2
                SyM = (SyM - transpose(SyM))/2
        end
        println("S Error: ", norm(Rx - SxM*V) + norm(Ry - SyM*V))

        Mi = pinv(M)
        Dx, Dy =  Mi*(SxM + Bx/2), Mi*(SyM + By/2)
        Op.P = M
        Op.Pi = Mi
        Op.Dx, Op.Dy = Dx, Dy
        Op.Bx, Op.By = Bx, By

        return Op.P, Op.Dx, Op.Dy, Op.Bx, Op.By
end


# Moves operator around, first a rotation around phi radians and second a translation 
# of (a, b). 
function moveOp(Op, a, b, phi)
        Rmat = [cos(phi) sin(phi);                                                              # Matrix describing the rotation. 
               -sin(phi) cos(phi)]
        mDom = Domain(Op.Omega.maxv, x->Op.Omega.CharFun(inv(Rmat)*x-[a, b]),                   # A new struct Domain. The characteristic function and parametrization
                t->Rmat*(Op.Opmega.param(t)+[a, b]), t->Rmat*Op.Omega.dparam(t),                # are transformed accordingly
                t->Rmat*Op.nparam(t), Op.Omega.discs)
        rpoints = zeros(Op.No + Op.Ni, 2)                                                       # The nodes and their normals are also transformed
        rnormals = zeros(Op.No, 2)
        for k=1:(Op.No + Op.Ni)
                rpoints[k, :] = Rmat*(Op.cpoints[k, :]) + [a, b]                                
        end
        for k=1:(Op.No)
                rnormals[k, :] = Rmat*Op.normals[k, :]
        end

        Dx = cos(phi)*Op.Dx + sin(phi)*Op.Dy                                                    # The same holds for the differentiation matrices
        Dy = -sin(phi)*Op.Dx + cos(phi)*Op.Dy                                                   # as any directional derivative can be expressed as
        Bx = cos(phi)*Op.Bx + sin(phi)*Op.By                                                    # a linear combination of directional derivatives in the 
        By = -sin(phi)*Op.Bx + cos(phi)*Op.By                                                   # original coordinate system.

                                                                                                # Warning: the function basis are not rotated.
        mOp = FSBP(mDom, Op.No, Op.Ni, rpoints, rnormals, Op.funcs, Op.dxfuncs, Op.dyfuncs,     # Therefore should the matrices P, Bx, By, Dx and Dy be constructed first
                        Op.SxEl, Op.SyEl, Op.BxEl, Op.ByEl, Op.P, Op.Pi,  # before the operator is rotated.
                        Dx, Dy, Bx, By)
        return mOp                                                                              # Return the new operator struct.
end

