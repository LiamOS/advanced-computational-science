#!/usr/bin/env python3
# Program to solve the wave equation with c=1 and external potential using
# central differences. Neumann boundaries are used.
# Liam O'Sullivan
from matplotlib.pyplot import *
from numpy import *

def run(dx,Time,Dist):
    dt = dx/sqrt(2)
    # Define variables with the number of space/time steps for ease later.
    # They must be cast as ints for use with range().
    nt = int((Time[1] - Time[0])/dt)
    nx = int((Dist[1] - Dist[0])/dx)
    # Define our initial arrays.
    U = zeros((nt, nx), dtype=float64)
    V = zeros(nx, dtype=float64)
    t = zeros(nt, dtype=float64)
    x = zeros(nx, dtype=float64)
    # Initialise the tand x arrays.
    for n in range(0,nt):
        t[n] = Time[0] + n*dt
    for j in range(0,nx):
        x[j] = Dist[0] + j*dx
    # Initialise the t_0 part of the array; Gaussian with sigma=1 and mu=6.
    for j in range(1,nx-1):
        U[0][j] = e**(-((x[j]-6)**2)/2)
    # Initialise our potential array.
    for j in range(0,nx):
        V[j] = 1/(cosh(x[j])**2)

    # c is not c from the wave equation, but nu^2.
    c = (dt/dx)**2
    # I manually impose the initial derivative condition, by computing the second
    # time step as d/dx of the analytic u_0, using the Euler method(Forgive me).
    for j in range(1,nx-1):
        U[1][j] = U[0][j] + dt*(-(x[j]-6)*(e**(-((x[j]-6)**2)/2)))
        U[1][0] = U[1][1]
        U[1][nx-1] = U[1][nx-2]
    # Now perform our integration.
    for n in range(2,nt):
        for j in range(1,nx-1):
            U[n][j] = c*(U[n-1][j-1]-2*U[n-1][j]+U[n-1][j+1])-U[n-2][j]+2*U[n-1][j]-(dt**2)*V[j]*U[n-1][j]
        U[n][0] = U[n][1]
        U[n][nx-1] = U[n][nx-2]
    Ltime = 18.*nt/(Time[1] - Time[0]) # The time t_{max}.
    norm = 0
    for b in U[Ltime]:
        norm += b**2
    norm = sqrt(dx*norm)
    return norm

# Define our space and time step sizes.
dxs = linspace(0.001,0.751,700)
# Define our Temporal and Spatial lengths.
aTime = [0,19]
aDist = [-20,20]
hs = []
for i in range(0,len(dxs)):
    n1 = run(dxs[i],aTime,aDist)
    n2 = run(dxs[i]/2.,aTime,aDist)
    hs.append(1 - n2/n1)
    print("h = "+str(round(dxs[i],4))+", err = "+str(hs[-1]))
semilogy(dxs,hs,'.')
xlabel("$h$")
ylabel("$\epsilon(h)$")
savefig("epsilon.pdf",format="pdf")
show()
