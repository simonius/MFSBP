# This file, when included, constructs a multidimensional SBP operator exact for polynomials of degree up to 3
# The operator is named simply Op and is constructed for a circular disk.


# inclusion of needed libraries
using QuadGK                                                                            # For numerical quadrature of the exact values for inner products
using LinearAlgebra                                                                     # Penrose-Inverse and SVD

include("FSBPcon2.jl")                                                                  # FSBP construction library
include("ptgen.jl")                                                                     # point insertion series.
include("bfs.jl")                                                                       # basis functions for function spaces
include("pfs.jl")                                                                       # plotting routines

xf(x) = norm(x - [0.5, 0.5]) <= 0.5 ? 1.0 : 0.0                                         # Characteristic function of the domain

circparam(t) = [0.5, 0.5] + 0.5 * [cos(2*pi*t), sin(2*pi*t)]                            # parametrization of the outer parameter with domain [0.0, 1.0]
dtcircparam(t) = 0.5 * [-2*pi*sin(2*pi*t), 2*pi*cos(2*pi*t)]                            # derivative of this parametrization
upper(x) = 0.5 + 0.5*sqrt(1-(2*(x-0.5))^2)                                              # parametrization of the upper half as x->(x, upper(x))
lower(x) = 0.5 - 0.5*sqrt(1-(2*(x-0.5))^2)                                              # and lower half as x->(x, lower(x)) for integration of the moments

RotMat = [0 1;                                                                          # rotates a tangent into a normal
        -1 0 ]

normalparam(t) = RotMat*dtcircparam(t)                                                  # parametrization of the normal, i.e. the rotated derivative of the parametrization of the circle

# The used function space, defined by spanning set of functions. Included are constants (cons = 1), linear monomials (linx = x, liny = y), 
# the quadratic monomials (q1 = x^2, q2 = xy, q3 = y^2)  and the cubic monomials, characterised by their multiindices 
# (poly(0.0, 3.0) = y^3, poly(1.0, 2.0) = xy^2, poly(2.0, 1.0) = x^2y, poly(3.0, 0.0) = x^3)
fcts = [cons, linx, liny, q1, q2, q3, poly(0.0, 3.0), poly(1.0, 2.0), poly(2.0, 1.0), poly(3.0, 0.0)]

# The derivatives in x and y direction of the spanning set are also needed to form the space of derivatives of the
# squared function space.
dxfcts = [x->0.0, dxlinx, dxliny, dxq1, dxq2, dxq3, dxpoly(0.0, 3.0), dxpoly(1.0, 2.0), dxpoly(2.0, 1.0), dxpoly(3.0, 0.0)]
dyfcts = [x->0.0, dylinx, dyliny, dyq1, dyq2, dyq3, dypoly(0.0, 3.0), dypoly(1.0, 2.0), dypoly(2.0, 1.0), dypoly(3.0, 0.0)]


No = 12                                                                                 # Number of nodes on the boundary (from now termed outer nodes)

Ni = 5*5                                                                                # Additional nodes in the interrior (from now termed inner nodes)

Dom = Domain([1.0, 1.0], xf, circparam, dtcircparam, normalparam, [0.0, 1.0])           # Save the definition of the domain into the struct Dom
opoints = zeros(No, 2)                                                                  # Declare space for the outer nodes and their normals
normals = zeros(No, 2)


for k=1:No
        opoints[k, :] = circparam((k - 0.5) / No)                                       # Equal spaced nodes on the boundary
        normals[k, :] = RotMat*dtcircparam((k-0.5)/No)                                  # with normal vectors
        normals[k, :] = normals[k, :] / norm(normals[k, :])                             # that are normed.
end

ipoints = EnterHalton(2, 3,xf, [1.0, 1.0], Ni)                                          # Definition of the inner nodes

pts = vcat(opoints, ipoints)                                                            # concat the inner and outer nodes into one list


# The following functions return values of the surface inner products in x and y direction using numerical quadrature.
function BxEl(k, l)
        integ(x) = fcts[k](circparam(x))*fcts[l](circparam(x))*dtcircparam(x)[2]	# Define the integrand
        val, err = quadgk(integ, 0.0, 1.0, maxevals=10^4)				# integrate
        if err > 1.0E-7
                print("Warning inexact quadrature moment")				# and warn if the estimated error is higher than 1.0E-7
        end
        return val
end

function ByEl(k, l)
        integ(x) = -fcts[k](circparam(x))*fcts[l](circparam(x))*dtcircparam(x)[1]
        val, err = quadgk(integ, 0.0, 1.0, maxevals=10^4)
        if err > 1.0E-7
                print("Warning inexact quadrature moment")
        end
        return val
end

# These two functions return the values of the volume inner product of the x and y derivatives of the spanning functions.
# For efficiency reasons the integrals were converted to integrals over the the x and y axis using the fact, 
# that one integration can be carried out as the integral function of a derivative is (easily) known.
# Therefore, these are as BxEl and ByEl.
function Sx(k, l)
        integ(y) = fcts[k]([upper(y), y])*fcts[l]([upper(y), y]) - fcts[k]([lower(y), y])*fcts[l]([lower(y), y])
        val, err = quadgk(integ, 0.0, 1.0, maxevals = 10^4)
        if err > 1.0E-7
                print("Warning inexact quadrature moment")
        end
        return val
end

function Sy(k, l)
        integ(x) = fcts[k]([x, upper(x)])*fcts[l]([x, upper(x)]) - fcts[k]([x, lower(x)])*fcts[l]([x, lower(x)])
        val, err = quadgk(integ, 0.0, 1.0, maxevals = 10^4)
        if err > 1.0E-7
                print("Warning inexact quadrature moment")
        end
        return val
end

# Using the constructor of struct FSBP, the needed information for a MFSBP operator is save in Op
Op = FSBP(Dom, No, Ni, pts, normals, fcts, dxfcts, dyfcts, Sx, Sy, BxEl, ByEl, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
# And the construction carried out.
M, Dx, Dy, Bx, By = MakeFSBPbB(Op, 3000)

