# This file, when included, constructs a multidimensional SBP operator exact for polynomials of degree up to 3
# The operator is named simply Op and is constructed for a reference triangle with vertices (0.0, 0.0) (1.0, 0.0), (0.0, 1.0).

# inclusion of needed libraries
using QuadGK									# For numerical quadrature of the exact values for inner products

include("ptgen.jl")								# point insertion series.
include("FSBPcon2.jl")								# FSBP construction library
include("bfs.jl")								# basis functions for function spaces

xf(x) = x[1] + x[2] <= 1.0 && x[1] >= 0 && x[2] >= 0 ? 1.0 : 0.0		# Characteristic function of the domain

# parametrization of the outer parameter with domain [0.0, 1.0]
function triparam(t)								
        if t < 1.0/3.0
                lp = 3.0*t
                return [lp, 0.0]
        elseif t < 2.0/3.0
                lp = 3.0*t-1.0
                return [1.0 - lp, lp]
        else
                lp = 3.0*t-2.0
                return [0.0, 1.0-lp]
        end
end

# derivative of this parametrization
function dtparam(t)
          if t < 1.0/3.0
                return [3.0, 0.0]
        elseif t < 2.0/3.0
                return [-3.0, 3.0]
        else
                return [0.0, -3.0]
        end
end

# rotates a tangent into a normal
RotMat = [0 1;
        -1 0 ]

# parametrization of the normal, i.e. the rotated derivative of the parametrization of the circle
normalparam(t) = RotMat*dtparam(t)

# The used function space, defined by spanning set of functions. Included are constants (cons = 1), linear monomials (linx = x, liny = y),
# the quadratic monomials (q1 = x^2, q2 = xy, q3 = y^2)  and the cubic monomials, characterised by their multiindices
# (poly(0.0, 3.0) = y^3, poly(1.0, 2.0) = xy^2, poly(2.0, 1.0) = x^2y, poly(3.0, 0.0) = x^3)
fcts = [cons, linx, liny, q1, q2, q3, poly(0.0, 3.0), poly(1.0, 2.0), poly(2.0, 1.0), poly(3.0, 0.0)]
dxfcts = [x->0.0, dxlinx, dxliny, dxq1, dxq2, dxq3, dxpoly(0.0, 3.0), dxpoly(1.0, 2.0), dxpoly(2.0, 1.0), dxpoly(3.0, 0.0)]
dyfcts = [x->0.0, dylinx, dyliny, dyq1, dyq2, dyq3, dypoly(0.0, 3.0), dypoly(1.0, 2.0), dypoly(2.0, 1.0), dypoly(3.0, 0.0)]


No = 24							# Number of nodes on the boundary (from now termed outer nodes)
Ni = 21							# Additional nodes in the interrior (from now termed inner nodes)

# Save the definition of the domain into the struct Dom
Dom = Domain([1.0, 1.0], xf, triparam, dtparam, normalparam, [0.0, 1.0/3.0, 2.0/3.0, 1.0])
opoints = zeros(No, 2)					# Declare space for the outer nodes and their normals
normals = zeros(No, 2)

for k=1:No
        opoints[k, :] = triparam((k - 0.5) / No)		# Equal spaced nodes on the boundary, respecting rotation symmetry
        normals[k, :] = RotMat*dtparam((k-0.5)/No)		# with normal vectors
        normals[k, :] = normals[k, :] / norm(normals[k, :])	# that are normed.
end

ipoints = EnterHalton(2, 3,xf, [1.0, 1.0], Ni)			# Definition of the inner nodes

pts = vcat(opoints, ipoints)					# concat the inner and outer nodes into one list

# The following functions return values of the surface inner products in x and y direction using numerical quadrature.
function BxEl(k, l)
        integ(x) = fcts[k](triparam(x))*fcts[l](triparam(x))*dtparam(x)[2] 	# Define the integrand
        val, err = quadgk(integ, 0.0,1/3,2/3, 1.0, maxevals=10^6)		# integrate
        if err > 10E-12
                print("Warning inexact quadrature moment")
        end
        return val
end

function ByEl(k, l)
        integ(x) = -fcts[k](triparam(x))*fcts[l](triparam(x))*dtparam(x)[1]
        val, err = quadgk(integ, 0.0, 1/3, 2/3, 1.0, maxevals=10^6)
        if err > 10E-12
                print("Warning inexact quadrature moment")
        end
        return val
end

# These two functions return the values of the volume inner product of the x and y derivatives of the spanning functions.
# For efficiency reasons the integrals were converted to integrals over the the x and y axis using the fact, 
# that one integration can be carried out as the integral function of a derivative is (easily) known.
# Therefore, these are as BxEl and ByEl.
function Sx(k, l)
        integ(y) = fcts[k]([1.0-y, y])*fcts[l]([1.0-y, y]) - fcts[k]([0.0, y])*fcts[l]([0.0, y])
        val, err = quadgk(integ, 0.0, 1.0, maxevals = 10^6)
        if err > 10E-7
                print("Warning inexact quadrature moment")
        end
        return val
end

function Sy(k, l)
        integ(x) = fcts[k]([x, 1.0-x])*fcts[l]([x, 1.0-x]) - fcts[k]([x, 0.0])*fcts[l]([x, 0.0])
        val, err = quadgk(integ, 0.0, 1.0, maxevals = 10^6)
        if err > 10E-7
                print("Warning inexact quadrature moment")
        end
        return val
end

# Using the constructor of struct FSBP, the needed information for a MFSBP operator is save in Op
PolyTriOp = FSBP(Dom, No, Ni, pts, normals, fcts, dxfcts, dyfcts, Sx, Sy, BxEl, ByEl, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
# And the construction carried out.
M, Dx, Dy, Bx, By, res = MakeFSBPmimeticB(PolyTriOp, 30000)
