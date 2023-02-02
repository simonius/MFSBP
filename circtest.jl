# File implements a linear Advection test using a FSBP operator on a circular disc. Obviously only a single block is (and can) be used.
# Boundary conditions are implemented using a local Lax-Friedrichs flux.

# inclusion of needed libraries
include("circOp.jl")						# inclusion of a SBP operator for a circular disk, exact for polynomials of degree smaller than 4
include("fcircOp.jl")						# inclusion of a FSBP operator for a circlar disk, exact for polynomials of degree 1 and travelling waves
								# with direction v and wavelength TL/||v||
include("pfs.jl")						# needed plotting functions
using DifferentialEquations					# time integration routines


TL = 2.0							# Make sure the typical length is suitable for the size of the single block. 
u0func(x) = sin(pi*dot(v, x)/TL)				# initial condition
ufunc(x, t) = u0func(x - v*t)					# solution
kufunc(x, y, t) = ufunc([x, y], t)

# driver routine for a numerical test using operator Op
# v is the advection direction and dt the timestep, tmax the end time.
function CircAdvNT(Op, v = [1.0, 1.0]; dt = 0.001, tmax = 1.0)		
	tspan = (0.0, tmax)						
	u0ar = zeros(Op.Ni + Op.No)				# arrays holding the initial condition
        Grid = zeros(Op.Ni + Op.No, 2)				# and the Grid coordinates for plotting 
	for i=1:(Op.Ni + Op.No)					# Initialization of both arrays.
        	Grid[i, :] = Op.cpoints[i, :]
        	u0ar[i] = u0func(Grid[i, :])
        end
	
	p = (Op, v, ufunc)

	println("Solving the Problem")				
	prob = ODEProblem(LinAdvCirc!, u0ar, tspan, p)
	sol = solve(prob, SSPRK33(), dt=dt, adaptive=false)	# Solve the problem
	return p, sol						# And return the result.
end

# The needed analytic fluxes
fLin(u, v) = v[1]*u
gLin(u, v) = v[2]*u
# And the corresponding (local) Lax-Friedrichs fluxes.
fLinLLF(ui, uo, v, n) = 0.5*(fLin(ui, v) + fLin(uo, v) + abs(v[1])*(ui - uo)*sign(n[1]))
gLinLLF(ui, uo, v, n) = 0.5*(gLin(ui, v) + gLin(uo, v) + abs(v[2])*(ui - uo)*sign(n[2]))

# Single block semidisretisation of a linear Advection equation using the LF flux above.
function LinAdvCirc!(du, u, p, t)
	du[:] = -p[2][1]*p[1].Dx*u - p[2][2]*p[1].Dy*u		# Calculate the volume term
	ffluxdiff = zeros(p[1].No + p[1].Ni)			# Reserve space for the flux over the surface
	gfluxdiff = zeros(p[1].No + p[1].Ni)
	for l=1:p[1].No						# Calculate the boundary terms.
		x = p[1].cpoints[l, :]
		ffluxdiff[l] = fLin(u[l], p[2]) - fLinLLF(u[l], p[3](x, t), p[2], p[1].normals[l, :])
		gfluxdiff[l] = gLin(u[l], p[2]) - gLinLLF(u[l], p[3](x, t), p[2], p[1].normals[l, :])
	end
	du[:] += p[1].Pi*(p[1].Bx*ffluxdiff + p[1].By*gfluxdiff) # And add the boundary terms.
end

# First, the testcase is solved using a polynomial SBP operator.
clf()
p, sol = CircAdvNT(Op)
# And the results are shown as a heatmap
tricontourf(Op.cpoints[:, 1], Op.cpoints[:, 2], sol(1.0), cmap="Reds")
cbar = colorbar()
cbar.formatter.set_powerlimits((0, 0))
xlabel("x")
ylabel("y")
PyPlot.tight_layout()

PyPlot.savefig("CircSol.jpg")

# Errors are calculated using the analytic solution
clf()
refvec = kufunc.(Op.cpoints[:, 1], Op.cpoints[:, 2], 1.0)
tricontourf(Op.cpoints[:, 1], Op.cpoints[:, 2], sol(1.0) .- refvec, cmap="Reds")
cbar = colorbar()
cbar.formatter.set_powerlimits((0, 0))
xlabel("x")
ylabel("y")
PyPlot.tight_layout()

PyPlot.savefig("CircSolDist.jpg")

# And the analyic solution is plotted using a high number of nodes
Np = 10000
clf()
ppoints = EnterHalton(2, 3,xf, [1.0, 1.0], Np)
refvec = kufunc.(ppoints[:, 1], ppoints[:, 2], 1.0)
tricontourf(ppoints[:, 1], ppoints[:, 2], refvec, cmap="Reds")
cbar = colorbar()
cbar.formatter.set_powerlimits((0, 0))
xlabel("x")
ylabel("y")
PyPlot.tight_layout()

PyPlot.savefig("CircRef.jpg")

# Then the same is carried out once more for the MFSBP operator
clf()
p, sol = CircAdvNT(fcOp)
tricontourf(fcOp.cpoints[:, 1], fcOp.cpoints[:, 2], sol(1.0), cmap="Reds")
cbar = colorbar()
cbar.formatter.set_powerlimits((0, 0))
xlabel("x")
ylabel("y")
PyPlot.tight_layout()

PyPlot.savefig("fCircSol.jpg")

clf()
refvec = kufunc.(fcOp.cpoints[:, 1], fcOp.cpoints[:, 2], 1.0)
tricontourf(fcOp.cpoints[:, 1], fcOp.cpoints[:, 2], sol(1.0) .- refvec, cmap="Reds")
cbar = colorbar()
cbar.formatter.set_powerlimits((0, 0))
xlabel("x")
ylabel("y")
PyPlot.tight_layout()

PyPlot.savefig("fCircSolDist.jpg")

