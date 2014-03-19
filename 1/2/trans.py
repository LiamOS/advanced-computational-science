#!/usr/bin/env python3
# Program to solve the wave equation with c=1 and external potential using
# central differences. Neumann boundaries are used.
# This program computes the ratio of transmitted to reflected wave.
# Liam O'Sullivan
from matplotlib.pyplot import *
from numpy import *

def run(sig,dx,Time,Dist):
    # Define our space and time step sizes.
    dt = dx/sqrt(2)
    # Define variables with the number of space/time steps for ease later.
    # They must be cast as ints for use with range().
    nt = int((Time[1] - Time[0])/dt) + 1
    nx = int((Dist[1] - Dist[0])/dx) + 1
    # Define our initial arrays.
    U = zeros((nt, nx), dtype=float64)
    H = zeros((nt, nx), dtype=float64)
    V = zeros(nx, dtype=float64)
    t = zeros(nt, dtype=float64)
    x = zeros(nx, dtype=float64)

    # Initialise the tand x arrays.
    for n in range(0,nt):
        t[n] = Time[0] + n*dt
    for j in range(0,nx):
        x[j] = Dist[0] + j*dx
    # Initialise the t_0 part of the array; Gaussian with sigma=sig and mu=6.
    for j in range(1,nx-1):
        U[0][j] = e**(-((x[j]-6)**2)/(2*(sig**2)))
    # Initialise our potential array.
    for j in range(0,nx):
        V[j] = 1/(cosh(x[j])**2)
    # c is not c from the wave equation, but nu^2.
    c = (dt/dx)**2
    # I manually impose the initial derivative condition, by computing the second
    # time step as d/dx of the analytic u_0, using the Euler method(Forgive me).
    for j in range(1,nx-1):
        U[1][j] = U[0][j] + dt*(-((x[j]-6)/(sig**2))*(e**(-((x[j]-6)**2)/(2*(sig**2)))))
    # Now perform our integration.
    for n in range(2,nt):
        for j in range(1,nx-1):
            U[n][j] = c*(U[n-1][j-1]-2*U[n-1][j]+U[n-1][j+1])-U[n-2][j]+2*U[n-1][j]-(dt**2)*V[j]*U[n-1][j]
            # Impose the Neumann condition here:
        U[n][0] = U[n][1]
        U[n][nx-1] = U[n][nx-2]
    # No need to worry about t=0,tmax boundaries, start compuating H generally.
    for n in range(1,nt-1):
        # Sort out the j boundaries, ignore potential as it's really small.
        H[n][0] = 0.5*(((1/(2*dt))*(U[n+1][0] - U[n-1][0]))**2) + 0.5*(((1/dx)*(U[n][1] - U[n][0]))**2)
        H[n][nx-1] = 0.5*(((1/(2*dt))*(U[n+1][nx-1]-U[n-1][nx-1]))**2)+0.5*(((1/dx)*(U[n][nx-1]-U[n][nx-2]))**2)
        for j in range(1,nx-1):
            H[n][j] = 0.5*(((1/(2*dt))*(U[n+1][j]-U[n-1][j]))**2) + 0.5*(((1/(2*dx))*(U[n][j+1]-U[n][j-1]))**2) + 0.5*V[j]*(U[n][j]**2)
    Tmax = 18.*nt/(Time[1] - Time[0]) # Find the index for t_max.
    jmid = int(nx/2)
    T = dx*sum(H[Tmax][0:jmid])
    E = dx*sum(H[Tmax])
    return T/E

aTime = [0,19]
aDist = [-20,20]
sigmas = []
Ts = []
for s in linspace(4,0.4,200):
     sigmas.append(s)
     Ts.append(run(s,0.05,aTime,aDist))
     print("Sigma = "+str(round(s,5))+", T = "+str(round(Ts[-1],5)))

plot(sigmas,Ts,'-')
ylabel("$T(\sigma)$")
xlabel("$\sigma$")
savefig("sigma.pdf",format="pdf")
show()
