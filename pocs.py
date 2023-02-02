import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})

import numpy as np

# Function returns an array of points from the pocs algorithm
# x0 is the starting point, p1 and p2 are the two projections
# and N the number of iterations.
def POCS(x0, p1, p2, N):
        xar = np.zeros([2, N])
        xar[:, 0] = x0
        for k in range(1, N):
                if k%2 == 1:
                        xar[:, k] = p1(xar[:, k-1])
                else:
                        xar[:, k] = p2(xar[:, k-1])
        
        return xar

# Projections onto the subspaces
# The two subspaces are x = 0 and y=x+1
def p1(x):
        return np.array([0.0, x[1]])

def p2(x):
        # Normal
        v = np.array([1.0, -1.0])
        # Support
        x0 = np.array([1.0, 2.0])
        delt = x - x0
        pdelt = delt - np.dot(v, delt)*v/np.dot(v, v)
        return pdelt + x0

x0 = np.array([-1.25, -0.25])
pw = POCS(x0, p1, p2, 100)

# Plotting the graphs
# First subspace
plt.plot([0.0, 0.0], [-0.5, 1.5], color="k")
plt.plot([-1.5, 0.5], [-0.5, 1.5], color="k")

plt.scatter(pw[0, :], pw[1, :], color="b")
for k in range(99):
    diff = np.linalg.norm(pw[:, k+1] - pw[:, k])
    plt.arrow(pw[0, k], pw[1, k], pw[0, k+1] - pw[0, k], pw[1, k+1] - pw[1, k], width = 0.03*diff, length_includes_head=True, fc="b", ec="b")

plt.scatter([0.0], [0.0], color="k")
plt.xlabel("x")
plt.ylabel("y")
ax = plt.gca()
ax.set_aspect('equal', 'box')
ax.text(-1.5, -0.3, "$C$")
ax.text(0.1, -0.4, "$D$")
ax.text(0.05, 0.0, "0")
plt.axis('off')
plt.savefig("pocsfig.pdf")

