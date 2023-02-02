# File implements a convergence analysis for an inhomogeneuous steady state problem on
# triangular blocks. The domain is formed by a grid of rectangular blocks 
# Every block is split via a diagonal from the top left to the bottom right corner to produce two triangles. Boundary conditions are
# implemented via a Lax-Friedrichs numerical flux and the strong form. The PDE is solved in the reference element element scale.


include("triangSteadyTest.jl")				# include the problem definition and discretisation from the linear advection test.
using Makie						# Plotting libraries
using GLMakie
using CairoMakie


minind = 6						# The minimum number of blocks per direction tested
maxind = 10						# The maximum number of blocks per direction tested
errors = zeros(maxind)					# the error norms of the solution by the FSBP operator. 
residuals = zeros(maxind)				# the residuals of the operator construction
perrors = zeros(maxind)					# the errors of the solution by the polynomial SBP operator

for k=minind:maxind
        						# set global variable for the function space
        global TL = k
        tend = 2*k					# and the corresponding end time in the time scale of the reference element. 
        global tau = k/10				# and the corresponding relaxation time 
	                                                # First a solution using a FSBP operator is calculated. 
        	                                        # As this operator depends on the size of its block in the grid we have to
                	                                # rebuild the operator for this particular grid, i.e. this TL.

        M, Dx, Dy, Bx, By, res = MakeFSBPmimeticB(TriOp, 40000)
        residuals[k] = res				# The residual is saved for later usage
        						# Solve the system
        p, sol = TriangStNT(TriOp, N =k, K = k, tend=tend, dt=0.0001*k)
        Gx = p.Grid[:, :, :, p.SBPOp.No+1:end, 1]
        Gy = p.Grid[:, :, :, p.SBPOp.No+1:end, 2]
        u = sol(tend)[:, :, :, p.SBPOp.No+1:end]
        uref = kufunc.(Gx, Gy, tend)
						        # For scaling we have to divide the error through the total node number.
        Ntot = length(uref)
        errors[k] = norm((u-uref))/sqrt(Ntot)		# Please note that we are by purpose calculating the per point error and are not using the inner product induced by the FSBP operator.
	                                                # This is done as otherwise the error for the second operator would be calculated with another norm, as the other operator uses another volume quadrature.

        	                                        # The same as before for the FSBP oeprator is done 
                	                                # for the poly. SBP operator.
                        	                        # Solve the problem
        pp, psol = TriangStNT(PolyTriOp, N =k, K = k, tend=tend, dt=0.0001*k)
        Gx = pp.Grid[:, :, :, pp.SBPOp.No+1:end, 1]
        Gy = pp.Grid[:, :, :, pp.SBPOp.No+1:end, 2]
        u = psol(tend)[:, :, :, pp.SBPOp.No+1:end]
        uref = kufunc.(Gx, Gy, tend)
        						# For scaling we have to divide the error through the total node number.
        Ntot = length(uref)
        perrors[k] = norm((u-uref))/sqrt(Ntot)		#  As we calculate an l^2 norm we only scale by the square root of the number of points
end

Narr = collect(minind:maxind)				# For convenience, the unused array positions are cut of.
errarr = errors[minind:maxind]				# the index where the order lines are "mounted" to the error markers

mi = 1
CairoMakie.activate!()
fontsize_theme = Theme(fontsize = 28, linewidth=3.0, markersize = 15.0) # overwrite the fontsize, linewidth and markersize for readable figures

set_theme!(fontsize_theme)
laxis = (xlabel = "blocks per direction", ylabel = "error", yscale=log10, xscale=log10)
fig = Makie.scatter(Narr, perrors[minind:maxind], axis = laxis, color=:black, label = "poly. MSBP", marker=:cross)	# plot the errors of the polynomial operators solution
f2 = Makie.scatter!(Narr, errarr, color=:black, label="MFSBP")	# plot the errors of the MFSBP operators solution

errarr = perrors[minind:maxind]				# add order 3 and 4 lines over the poly errors
f3 = lines!([Narr[1], Narr[end]], [errarr[mi]*(Narr[mi]/Narr[1])^3, errarr[mi]*(Narr[mi]/Narr[end])^3], label = "order 3", linestyle=:dot)

f4 = lines!([Narr[1], Narr[end]], [errarr[mi]*(Narr[mi]/Narr[1])^4, errarr[mi]*(Narr[mi]/Narr[end])^4], label = "order 4", linestyle=:dot)

errarr = errors[minind:maxind]				# add order 1 line over the FSBP operator error
f5 = lines!([Narr[1], Narr[end]], [errarr[mi]*(Narr[mi]/Narr[1])^1, errarr[mi]*(Narr[mi]/Narr[end])^1], label = "order 1")

							# add a legend.
Legend(fig.figure[1, 2], [fig.plot, f2, f5, f3, f4], ["poly. MSBP", "MFSBP", "order 1", "order 3", "order 4"])


save(string("stconvana.pdf"), fig)


