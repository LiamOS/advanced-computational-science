#!/usr/bin/env python3
# Program to solve the wave equation with c=1 using central differences.
# Dirichlet boundaries are used, and imposed implicitly, by allocating an array
# of zeros and not writing to the edge points. A variety of h values are used.
# Liam O'Sullivan
from matplotlib.pyplot import *
from numpy import *

def run(dx,T,D):
    # Define variables with the number of space/time steps for ease later.
    # They must be cast as ints for use with range(), but should be integers
    # anyway, since I'll only be setting dx as a nice value, I swear...
    dt = dx/sqrt(2)
    nt = int((T[1] - T[0])/dt)
    nx = int((D[1] - D[0])/dx)
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

    # c is not c from the wave equation, just a useful constant.
    c = (dt/dx)**2

    # Since U_j^n-1 is unknown, I take is equivalent to U_j^n for the first run.
    for j in range(1,nx-1):
        U[1][j] = c*(U[0][j-1] - 2*U[0][j] + U[0][j+1]) + U[0][j]

    # Now perform our norm computation.
    for n in range(2,nt):
        for j in range(1,nx-1):
            U[n][j] = c*(U[n-1][j-1] - 2*U[n-1][j] + U[n-1][j+1]) - U[n-2][j] + 2*U[n-1][j]

    Ltime = nt*10.5/(T[1] - T[0])
    norm = 0
    for k in U[Ltime]:
        norm += k**2
    norm = sqrt(dx*norm)
    return norm


# Define our Temporal and Spatial lengths.
Time = [0,14]
Dist = [-7,7]

# Define our space and time step sizes.
dxs = [(0.001 + o/2500) for o in range(0,2000)]
err = []

for i in range(0,len(dxs)):
    Lfull = run(dxs[i],Time,Dist)
    Lhalf = run(dxs[i]/2,Time,Dist)
    err.append(1 - (Lhalf/Lfull))
    print("h = "+str(round(dxs[i],4))+", err = "+str(err[-1]))

xlabel("$h$")
ylabel("$\epsilon(h)$")
semilogy(dxs,err,'.')
savefig("Error.pdf",format='PDF')
show()
