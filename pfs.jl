# This file holds all functions used for plotting solutions.
# A mix of PyPlot and Plots with the pyplot backend is used.

# Inclusion of needed libraries
using Plots
pyplot()
using PyPlot
using PyCall

# Enlarging the fontsize
rcParams = PyDict(matplotlib["rcParams"])
rcParams["font.size"] = 16


# A helper to transform the coordinates into the unit rectangle
RtoP(x) = x./maximum(x)

# This function plots the nodes on the interrior of every block. The restriction to the interrior is needed as
# the boundary nodes of adjacent blocks, owing the weak formulation of the boundary/coupling conditions, are not uiquely 
# determined. As argument an array of size N x K x (N_i + N_o)  has to be given, where N and K are the dimensions of the
# Grid, while the last axis holds the interrior nodes.
function InnerHeatmapAr(p, ar, name, arrname)
        Gx = p.Grid[:, :, p.SBPOp.No+1:end, 1]                                  # Collect the interrior node coordinates in x
        Gy = p.Grid[:, :, p.SBPOp.No+1:end, 2]                                  # and y.

        u = ar[:, :, p.SBPOp.No+1:end]                                          # Collect the solution values
        clf()                                                                   # Clear the canvas
	tricontourf(RtoP(Gx[:]), RtoP(Gy[:]), u[:], 200, cmap="Reds")           # plot a triangulated contour
        cbar = colorbar()                                                       # add a colorbar and name
        cbar.set_label(arrname)
        xlabel("x")                                                             # together with axis labels
        ylabel("y")
	PyPlot.tight_layout()                                                   # trim the layout

        PyPlot.savefig(name)                                                    # and save the figure
end

# The previous function exists once more not for meshes built from rectangles, but from triangles. 
# the implementation is otherwise identical to InnerHeatmapAr
function InnerTriHeatmapAr(p, ar, name, arrname)
        Gx = p.Grid[:, :, :, p.SBPOp.No+1:end, 1]
        Gy = p.Grid[:, :, :, p.SBPOp.No+1:end, 2]
        u = ar[:, :, :, p.SBPOp.No+1:end]
        
        clf()
	tricontourf(RtoP(Gx[:]), RtoP(Gy[:]), u[:], 200, cmap="Reds")
        cbar = colorbar()
        cbar.set_label(arrname)
        xlabel("x")
        ylabel("y")
	PyPlot.tight_layout()

        PyPlot.savefig(name)

end

# Equivalent to InnerHeatmapAr, but plots a solution object sol instead of an array.
function InnerHeatmap(p, sol, t, name, ref = (x, y, t) -> 0.0, subref = false)
        Gx = p.Grid[:, :, p.SBPOp.No+1:end, 1]
        Gy = p.Grid[:, :, p.SBPOp.No+1:end, 2]
        u = sol(t)[:, :, p.SBPOp.No+1:end]
        if subref
                uref = kufunc.(Gx, Gy, t)
                u = u - uref
        end

        clf()
	tricontourf(RtoP(Gx[:]), RtoP(Gy[:]), u[:], 200, cmap="Reds")
        cbar = colorbar()
        xlabel("x")
        ylabel("y")
        PyPlot.tight_layout()

        PyPlot.savefig(name)
end

# Triangle version of InnerHeatmap
function InnerTriHeatmap(p, sol, t, name, ref = (x, y, t) -> 0.0, subref = false)
        Gx = p.Grid[:, :, :, p.SBPOp.No+1:end, 1]
        Gy = p.Grid[:, :, :, p.SBPOp.No+1:end, 2]
        u = sol(t)[:, :, :, p.SBPOp.No+1:end]
        if subref
                uref = kufunc.(Gx, Gy, t)
                u = u - uref
        end

        clf()
	tricontourf(RtoP(Gx[:]), RtoP(Gy[:]), u[:], 200, cmap="Reds")
        cbar = colorbar()
        cbar.formatter.set_powerlimits((0, 0))                                  # Enforce scientific notation for the colorbar
        xlabel("x")
        ylabel("y")
	PyPlot.tight_layout()
        PyPlot.savefig(name)

end


# This function produces diagrams of the node distributions of operator Op and saves them under "name".
# ms is the markesize for the markers indicating the nodes. 
function PlotPts(Op, name; ms=1)
        prange = range(0.0, 1.0, length=1000)                                                           # a parameter range to sample the boundary of the domain
        xkord = t->Op.Omega.param(t)[1]                                                                 # the corresponding coordinate functions
        ykord = t->Op.Omega.param(t)[2]
        clf()
        PyPlot.plot(xkord.(prange), ykord.(prange), c="k", label="Boundary")                            # plot the outer boundary
        xlabel("x")
        ylabel("y")
        PyPlot.scatter(Op.cpoints[:, 1], Op.cpoints[:, 2], c="k", label="Collocation Points", s=ms)     # and add the collocation points.
        PyPlot.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2)                               # place a legend
        ax = PyPlot.gca()
        ax.set_aspect("equal", "box")
        PyPlot.savefig(name)                                                                            # save the results.
end
