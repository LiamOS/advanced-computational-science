#!/usr/bin/env python3
# Program to solve the wave equation with c=1 using central differences.
# Dirichlet boundaries are used, and imposed implicitly, by allocating an array
# of zeros and not writing to the edge points.
# Liam O'Sullivan

from matplotlib.pyplot import *
from numpy import *

# Define our space and time step sizes.
dx = 0.004
dt = dx

# Define our Temporal and Spatial lengths.
Time = [0,14]
Dist = [-7,7]

# Define variables with the number of space/time steps for ease later.
# They must be cast as ints for use with range(), but should be integers
# anyway, since I'll only be setting dx as a nice value, I swear...
nt = int((Time[1] - Time[0])/dt)
nx = int((Dist[1] - Dist[0])/dx)

# Define our initial arrays.
U = zeros((nt, nx), dtype=float64)
t = zeros(nt, dtype=float64)
x = zeros(nx, dtype=float64)
# Initialise the tand x arrays.
for n in range(0,nt):
    t[n] = Time[0] + n*dt
for j in range(0,nx):
    x[j] = Dist[0] + j*dx

# Initialise the t_0 part of the array.
for j in range(1,nx-1):
    U[0][j] = e**(-x[j]*x[j])

# c is not c from the wave equation:
c = (dt/dx)**2

# Since U_j^n-1 is unknown, I take is equivalent to U_j^n for the first run.
for j in range(1,nx-1):
    U[1][j] = c*(U[0][j-1] - 2*U[0][j] + U[0][j+1]) + U[0][j]

# Now perform our integration.
for n in range(2,nt):
    for j in range(1,nx-1):
        U[n][j] = c*(U[n-1][j-1] - 2*U[n-1][j] + U[n-1][j+1]) - U[n-2][j] + 2*U[n-1][j]
    if ((n-2)%600) is 0:
        plot(x,U[n],label="t = "+str(round(14*(n/nt),1)))
xlabel("$x$")
ylabel("$u(x,t)$")
xlim([-7,7])
legend()
savefig("dirichlet.pdf",format="pdf")
show()
