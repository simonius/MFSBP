# This file defines basis functions used in the construction of the operators
# All basis functions take one two-dimensional vector as argument.

# A constant basis function
cons(x) = 1.0

# The linear monomial in x
linx(x) = x[1]
# and its derivative in x and
dxlinx(x) = 1.0
# y direction
dylinx(x) = 0.0

# The linear monomial in y
liny(x) = x[2]
# derivative in x
dxliny(x) = 0.0
# and y direction
dyliny(x) = 1.0

# the quadratic monomial with multiindex (2, 0)
q1(x) = x[1]^2
# and its derivatives
dxq1(x) = 2.0*x[1]
dyq1(x) = 0.0

# As above, but with multiindex (1, 1)
q2(x) = x[2]^2
dxq2(x) = 0.0
dyq2(x) = 2.0.*x[2]

# And multiindex (0, 2)
q3(x) = x[1]*x[2]
dxq3(x) = x[2]
dyq3(x) = x[1]

# For higher degree polynomials higher class functions are used, this means that the functions below return a 
# function on their own and are therefore able to express arbitrary monomials when the corresponding multiindices (k, l) are given

# dpow(x, k) returns the derivative of the monomial x^k 
dpow(x, k) = k > 0 ? k*x^(k-1) : 0.0

# poly(k, l) returns the function function that takes the first and second entry of x to the k-th and l-th power - and NOT a value
poly(k, l) = x->(x[1])^k*(x[2])^l
# dxpoly and dypoly return a function, that when called returns the derivative of x[1]^k*x[2]^l.
dxpoly(k, l) = x->dpow(x[1], k)*x[2]^l
dypoly(k, l) =  x->(x[1])^k*dpow(x[2], l)

# The exponential function of the first entry of x and its derivatives in x and y-direction
expx(x) = exp(x[1])
dxexpx(x) = exp(x[1])
dyexpx(x) = 0.0

# The exponential function of the second entry of x and its derivatives in x and y-direction 
expy(x) = exp(x[2])
dxexpy(x) = 0.0
dyexpy(x) = exp(x[2])

# The mirrored version of expx. This is needed as - in contrast to monomials - the exponential is not symmetric or antisymmetric.
# Therefore, only by the inclusion of the exponential can a suitable amount of symmetry be guaranteet. As before, derivatives are
# also defined.
mexpx(x) = exp(-x[1])
dxmexpx(x) = -exp(-x[1])
dymexpx(x) = 0.0

# mirrored version of expy
mexpy(x) = exp(-x[2])
dxmexpy(x) = 0.0
dymexpy(x) = -exp(-x[2])

# The next basis functions allow us differentiate waves in the v direction exactly. TL is the typical length of the waves, 
# while the wavelength is TL/(pi||v||). 
TL = 10.0
v = [1.0, 1.0]

sinsol(x) =  sin(pi*dot(v, x)/TL)
dxsinsol(x) = cos(pi*dot(v, x)/TL)*pi/TL*v[1]
dysinsol(x) = cos(pi*dot(v, x)/TL)*pi/TL*v[2]

cossol(x) =  cos(pi*dot(v, x)/TL)
dxcossol(x) = -sin(pi*dot(v, x)/TL)*pi/TL*v[1]
dycossol(x) = -sin(pi*dot(v, x)/TL)*pi/TL*v[2]

# If one wants instead to be exact for a set of fourier modes the following set of function definitions is helpfull.
# first, as before, a one-dimensional fourier basis mode is defined, together with its derivative.
fb1(k, x) = mod(k, 2) == 1 ? sin(2*pi*floor(k/2)*x) : cos(2*pi*floor(k/2)*x)
dfb1(k, x) = mod(k, 2) == 1 ? cos(2*pi*floor(k/2)*x) * 2*pi*floor(k/2) : -sin(2*pi*floor(k/2)*x)*2*pi*floor(k/2)

# The two-dimensional fourier basis function is defined as product of the previous one-dimensional case. 
# This is again a higher class function, i.e. fb(3, 4) returns a function with one argument that when called
# returns the fourier mode with multiindex (3, 4)
fb2(k, l) = x->fb1(k, x[1])*fb1(l, x[2])
dxfb2(k, l) = x->dfb1(k, x[1])*fb1(l, x[2])
dyfb2(k, l) = x->fb1(k, x[1])*dfb1(l, x[2]) 


# The first fourier modes are defined here explicitly for a steady state test with non-unit wavelength.
p1(x) =  sin(2*pi*x[1]/TL)*sin(2*pi*x[2]/TL)
dxp1(x) = cos(2*pi*x[1]/TL)*sin(2*pi*x[2]/TL)*2*pi/TL
dyp1(x) = sin(2*pi*x[1]/TL)*cos(2*pi*x[2]/TL)*2*pi/TL

p2(x) =  cos(2*pi*x[1]/TL)*cos(2*pi*x[2]/TL)
dxp2(x) = -sin(2*pi*x[1]/TL)*cos(2*pi*x[2]/TL)*2*pi/TL
dyp2(x) = -cos(2*pi*x[1]/TL)*sin(2*pi*x[2]/TL)*2*pi/TL

p3(x) = sin(2*pi*x[1]/TL)*cos(2*pi*x[2]/TL)
dxp3(x) = cos(2*pi*x[1]/TL)*cos(2*pi*x[2]/TL)*2*pi/TL
dyp3(x) = -sin(2*pi*x[1]/TL)*sin(2*pi*x[2]/TL)*2*pi/TL

p4(x) = cos(2*pi*x[1]/TL)*sin(2*pi*x[2]/TL)
dxp4(x) = -sin(2*pi*x[1]/TL)*sin(2*pi*x[2]/TL)*2*pi/TL
dyp4(x) = cos(2*pi*x[1]/TL)*cos(2*pi*x[2]/TL)*2*pi/TL

# Including constant versions and double frequency versions allows to also esactly
# differentiate squares, as these will result in those linear combinations.
psq1(x) = sin(2*pi*x[1]/TL)
dxpsq1(x) = cos(2*pi*x[1]/TL)*2*pi/TL
dypsq1(x) = 0.0

psq2(x) = sin(2*pi*x[2]/TL)
dxpsq2(x) = 0.0
dypsq2(x) = cos(2*pi*x[2]/TL)*2*pi/TL


psq3(x) = cos(2*pi*x[1]/TL)
dxpsq3(x) = -sin(2*pi*x[1]/TL)*2*pi/TL
dypsq3(x) = 0.0


psq4(x) = cos(2*pi*x[2]/TL)
dxpsq4(x) = 0.0
dypsq4(x) = -sin(2*pi*x[2]/TL)*2*pi/TL




