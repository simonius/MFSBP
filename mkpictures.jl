# This file serves as a batch processing script to produce test results


# Tests for the linear advection on triangular blocks.
include("triangtest.jl")

p, sol = TriangAdvNT(TriOp)
InnerTriHeatmap(p, sol, 10.0, "LinAdvTriang.jpg", kufunc)
InnerTriHeatmap(p, sol, 10.0, "LinAdvTriangDist.jpg", kufunc, true)

p, sol = TriangAdvNT(PolyTriOp)
InnerTriHeatmap(p, sol, 10.0, "PolyLinAdvTriang.jpg", kufunc)
InnerTriHeatmap(p, sol, 10.0, "PolyLinAdvTriangDist.jpg", kufunc, true)


# Tests for a steady state on triangular blocks.
include("triangSteadyTest.jl")
include("steadyTriOp.jl")
include("polyTriOp.jl")

p, sol = TriangStNT(TriOp)
InnerTriHeatmap(p, sol, 20.0, "Boltz.jpg", kufunc)
InnerTriHeatmap(p, sol, 20.0, "BoltzDist.jpg", kufunc, true)

p, sol = TriangStNT(PolyTriOp)
InnerTriHeatmap(p, sol, 20.0, "pBoltz.jpg", kufunc)
InnerTriHeatmap(p, sol, 20.0, "pBoltzDist.jpg", kufunc, true)

