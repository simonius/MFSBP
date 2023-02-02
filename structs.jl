# This file contains definitions of all structs used in 
# the implementation.

# FSBP operators are described using this struct.
# Only the last section of data is set in the construction procedure 
mutable struct FSBP
        Omega           # Domain definition, see struct below.
        No              # Number of nodes on the boundary
        Ni              # Number of nodes in the interrior
        cpoints         # (N_o + N_i) x 2 array of collocation point coordinates
        normals         #  N_0 x 2 array of normal directions for the surface nodes.

        funcs           # List of basis functions
        dxfuncs         # List of basis functions derived in x
        dyfuncs         # List of basis functions derived in y
    
        
        SxEl            # Moments of derivatives in x direction, Volume Integral
        SyEl            # Moments of derivatives in y direction, Volume Integral
        BxEl            # Moments on the surface, normal x component
        ByEl            # Moments on the surface, normal y component


        P               # Mass matrix
        Pi              # Inverted mass matrix
        Dx              # x derivative operator
        Dy              # y derivative operator 
        Bx              # x normal boundary operator
        By              # y normal boundary operator
end

# This struct describes the domain.
# It is assumed that the domain can be inscribed into a bounding box starting at (0.0).
# Please note that the last is only needed for efficient quadrature
mutable struct Domain
        maxv            # Top right end of bounding box
        CharFunc        # Characteristic function of the domain
        param           # Boundary parametrisation on (0.0, 1.0)
        dparam          # Derivative of boundary parametrisation
        nparam          # Parametrization of the outer normal
        discs           # Discontinuities in the derivatives of the parametrization
end

# Describes a linear advection with a single block.
mutable struct LinAdv
        vvec            # advection direction
        SBPOp           # Used SBP operator
        Bp              # Negative definite boundary operator
        Bn              # positive definite boundary operator
        sigma           # SAT constant
        ufunc           # u(x, t), initial condition and time evolution of the exact solution
end

# Describes linear advection with multiple rectangular blocks
# and boundary conditions. Uses struct LinAdv
mutable struct LinAdvMb
        N               #
        K               # Number of Elements in X and Y direction
        Adv             # Struct LinAdv that describes linear advection and base discretisation
        dt              # Timestep size, only used for classical LF flux.
        Grid            # Coordinates of all points in the Grid, N x K x (No+Ni) x 2
        Idp             # N x K x No x 3 Array of all neighboring points for all boundary nodes, i.e. 
                        # Idp[n_1, k_1, l_1, :] yields the node coordinates n_2, k_2, l_2 of the corresponding
                        # neighbor of node l_1 in block (n_1, k_1)
        SBPOp           # Convenience copy
end

# Describes burgers equation with multiple rectangular blocks
# and toroid periodic boundary conditions
mutable struct BurgersMB
        SBPOp           # Used MFSBP operator
        ufunc           # Unitial condition
        N               # Blocks in x
        K               # and y direction
        dt              # Timestep size
        Grid            # Coordinates of all points in the grid
        Idp             # # N x K x No x 3 Array of all neighboring points for all boundary nodes, i.e. 
                        # Idp[n_1, k_1, l_1, :] yields the node coordinates n_2, k_2, l_2 of the corresponding
                        # neighbor of node l_1 in block (n_1, k_1)
        split           # splitting constant of the split form volume term.
end


# Describes a multiblock discrectisation based on triangles
# In a regular grid of rectangles every rectangle is split
# from the top left to the lower right 
mutable struct TriangMB
        N               # Blocks in x
        K               # and y direction
        dt              # Timestep size
        ub              # boundary data
        Grid            # Coordinates of all points in the Grid, N x K x 2 x (No + Ni) x 2
                        # the third dimension enumerates if lower left or upper right traingle is used 
        SBPOp           # Original Operator, for lower left triangle
        MSBPOp          # Mirrored Operator, for upper right triangle
        vvec            # advection direction 
end
