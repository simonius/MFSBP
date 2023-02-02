# File implements an inhomogeneuous steady state problem on
# triangular blocks. The domain is formed by a grid of rectangular blocks 
# Every block is split via a diagonal from the top left to the bottom right corner to produce two triangles. Boundary conditions are
# implemented via a Lax-Friedrichs numerical flux and the strong form. The PDE is solved in the reference element element scale.

# inclusion of the needed operators. 
include("polyTriOp.jl")					# Exact for polynomials with a degree up to 3
include("steadyTriOp.jl")                               # Exact for the steady state and linear polynomials
include("pfs.jl")					# plotting routines 
using DifferentialEquations				# time integration methods

TL = 10							# the typical length of of the solution.
u0func(x) = 0.0						# zero as initial condition.
amp = 1.0						# amplitude of the steady state
tau = 1.0
ub(x, t) = amp*sin(2*pi*x[1]/TL)*sin(2*pi*x[2]/TL)	# boundary condition - the steady state at the boundary 
rhs(u, x, y) = amp*(dxp1([x, y]) + dyp1([x, y])) + (amp*p1([x, y]) - u)/tau # the needed right hand side
kufunc(x, y, t) = ub([x, y], t)				# used for error calculation


# Performs a numerical test for the steady state problem. Op is the used operator.
# N and K are the number of elements in x and y direction. v is the advection direction and dt the 
# timestep.
function TriangStNT(Op; N = 10, K = 10, v = [1.0, 1.0], dt = 0.01, tend=20.0)
        tspan = (0.0, tend)							
        u0ar = zeros(N, K, 2, Op.Ni + Op.No)					# The tensor initialized with the initial condition.
        Grid = zeros(N, K, 2, Op.Ni + Op.No, 2)					# A tensor with coordinates of all nodes.
        									# rotate one Triangle by pi and move it to the right upper corner
        turnedOp = moveOp(Op, 1.0, 1.0, pi) 

        #initialize grid
        for n=1:N
                for k=1:K
                        for i=1:(Op.Ni + Op.No)
                                Grid[n, k, 1, i, :] = [n-1, k-1] + Op.cpoints[i, :]
                                Grid[n, k, 2, i, :] = [n-1, k-1] + turnedOp.cpoints[i, :]
                                u0ar[n, k, 1, i] = u0func(Grid[n, k, 1, i, :])
                                u0ar[n, k, 2, i] = u0func(Grid[n, k, 2, i, :])
                        end
                end
        end

        p = TriangMB(N, K, dt, ub, Grid, Op, turnedOp, v)
        									# Solve the problem
        println("Solving the Problem")
        prob = ODEProblem(LinAdvLLF!, u0ar, tspan, p)
        sol = solve(prob, SSPRK33(), dt=dt)
        return p, sol								# and return the result
end


# The analytic fluxes for the linear advection
fLin(u, v) = v[1]*u
gLin(u, v) = v[2]*u
# The local Lax-Friedrichs flux for a linear advection with speed v = (v_1, v_2),
# internal value ui, external value uo and outer surface normal n.
fLinLLF(ui, uo, v, n) = 0.5*(fLin(ui, v) + fLin(uo, v) + abs(v[1])*(ui - uo)*sign(n[1]))
gLinLLF(ui, uo, v, n) = 0.5*(gLin(ui, v) + gLin(uo, v) + abs(v[2])*(ui - uo)*sign(n[2]))

# Discretisation of the steady state problem using a local Lax-Friedrichs flux for the coupling of a multiblock discretisation. 
function LinAdvLLF!(du, u, p, t)
        									# Calculation of the inner derivative
        numbound = floor(Int, p.SBPOp.No/3)
        for n=1:p.N
                for k=1:p.K
                        							# lower left triangle
			du[n, k, 1, :] = - p.vvec[1]*p.SBPOp.Dx*u[n, k, 1, :] - p.vvec[2]*p.SBPOp.Dy*u[n, k, 1, :]
										# upper right triangle
			du[n, k, 2, :] = - p.vvec[1]*p.MSBPOp.Dx*u[n, k, 2, :] - p.vvec[2]*p.MSBPOp.Dy*u[n, k, 2, :]
                end
        end
        									# Calculation of the boundary terms
        ffluxdiff = zeros(p.SBPOp.No + p.SBPOp.Ni)				# Flux differences in the f (x-direction)
        gfluxdiff = zeros(p.SBPOp.No + p.SBPOp.Ni)				# and g (y-direction) fluxes.
        for n=1:p.N								# iteration over all triangles
                for k=1:p.K
                        # First for the lower triangle. Sides are named
			# 1 - bottom,  2 - diagonal, 3 - left 
                        # first side: neighbor is upper triangle of lower cell, side 1
                        for l=1:numbound
                                if k!=1						# Use neighbor on the interrior
                                        ubs = u[n, k-1, 2, numbound-l+1]
                                else						# or boundary value at the boundary
                                        ubs = p.ub(p.Grid[n, k, 1, l, :], t)
                                end

                                ffluxdiff[l] = fLin(u[n, k, 1, l], p.vvec) - fLinLLF(u[n,k, 1, l], ubs, p.vvec, p.SBPOp.normals[l, :])
                                gfluxdiff[l] = gLin(u[n, k, 1, l], p.vvec) - gLinLLF(u[n,k, 1, l], ubs, p.vvec, p.SBPOp.normals[l, :])

                        end
                        							# second side: upper triangle of same cell, side 2
                        for l=1:numbound
                                ubs = u[n, k, 2, 2*numbound-l+1]
                                ffluxdiff[l+numbound] = fLin(u[n, k, 1, l+numbound], p.vvec) - fLinLLF(u[n,k, 1, l+numbound], ubs, p.vvec, p.SBPOp.normals[l+numbound, :])
                                gfluxdiff[l+numbound] = gLin(u[n, k, 1, l+numbound], p.vvec) - gLinLLF(u[n,k, 1, l+numbound], ubs, p.vvec, p.SBPOp.normals[l+numbound, :])

                        end
                        							# third side: upper triangle of left cell, side 3
                        for l=1:numbound
                                if n!=1
                                        ubs = u[n-1, k, 2, 3*numbound-l+1]
                                else
                                        ubs = p.ub(p.Grid[n, k, 1, l+2*numbound, :], t)
                                end
                                ffluxdiff[l+2*numbound] = fLin(u[n, k, 1, l+2*numbound], p.vvec) - fLinLLF(u[n,k, 1, l+2*numbound], ubs, p.vvec, p.SBPOp.normals[l+2*numbound, :])
                                gfluxdiff[l+2*numbound] = gLin(u[n, k, 1, l+2*numbound], p.vvec) - gLinLLF(u[n,k, 1, l+2*numbound], ubs, p.vvec, p.SBPOp.normals[l+2*numbound, :])
                        end

                        # add in boundary contributions
                        du[n, k, 1, :] += p.SBPOp.Pi*(p.SBPOp.Bx*ffluxdiff + p.SBPOp.By*gfluxdiff)



                        # Now for the upper triangle. Sides are numbered
			# 1 - upper side, 2 - diagonal, 3 - right side
                        # first side: neighbor is lower triangle of upper cell, side 1
                        for l=1:numbound
                                if k!=p.K
                                        ubs = u[n, k+1, 1, numbound-l+1]
                                else
                                         ubs = p.ub(p.Grid[n, k, 2, l, :], t)
                                end

                                ffluxdiff[l] = fLin(u[n, k, 2, l], p.vvec) - fLinLLF(u[n,k, 2, l], ubs, p.vvec, p.MSBPOp.normals[l, :])
                                gfluxdiff[l] = gLin(u[n, k, 2, l], p.vvec) - gLinLLF(u[n,k, 2, l], ubs, p.vvec, p.MSBPOp.normals[l, :])

                        end
                        # second side: neighbor is lower triangle of same cell, side 2
                        for l=1:numbound
                                ubs = u[n, k, 1, 2*numbound-l+1]
                                ffluxdiff[l+numbound] = fLin(u[n, k, 2, l+numbound], p.vvec) - fLinLLF(u[n,k, 2, l+numbound], ubs, p.vvec, p.MSBPOp.normals[l+numbound, :])
                                gfluxdiff[l+numbound] = gLin(u[n, k, 2, l+numbound], p.vvec) - gLinLLF(u[n,k, 2, l+numbound], ubs, p.vvec, p.MSBPOp.normals[l+numbound, :])

                        end
                        # third side: neighbor is lower triangle of right cell, side 3
                        for l=1:numbound
                                if n!=p.N 
                                        ubs = u[n+1, k, 1, 3*numbound-l+1]
                                else
                                        ubs = p.ub(p.Grid[n, k, 2, l+2*numbound, :], t)
                                end

                                ffluxdiff[l + 2*numbound] = fLin(u[n, k, 2, l+2*numbound], p.vvec) - fLinLLF(u[n,k, 2, l+2*numbound], ubs, p.vvec, p.MSBPOp.normals[l+2*numbound, :])
                                gfluxdiff[l + 2*numbound] = gLin(u[n, k, 2, l+2*numbound], p.vvec) - gLinLLF(u[n,k, 2, l+2*numbound], ubs, p.vvec, p.MSBPOp.normals[l+2*numbound, :])

                        end

                         # add in boundary contributions
                        du[n, k, 2, :] += p.MSBPOp.Pi*(p.MSBPOp.Bx*ffluxdiff + p.MSBPOp.By*gfluxdiff)

                end
        end
        # add in rhs
        du[:, :, :, :] += rhs.(u[:, :, :, :], p.Grid[:, :, :, :, 1], p.Grid[:, :, :, :, 2])
end






