#!/usr/bin/env python3
# Program to solve the wave equation with c=1 using central differences.
# Dirichlet boundaries are used, and imposed implicitly, by allocating an array
# of zeros and not writing to the edge points.
# This program proceeds to caltulate the energy as a function of time and h.
# Liam O'Sullivan

from matplotlib.pyplot import *
from numpy import *

def run(dx,Time,Dist):
    # Define our space and time step sizes.
    dt = dx/sqrt(2)

    # Define variables with the number of space/time steps for ease later.
    # They must be cast as ints for use with range().
    nt = int((Time[1] - Time[0])/dt) + 1
    nx = int((Dist[1] - Dist[0])/dx) + 1

    # Define our initial arrays.
    U = zeros((nt, nx), dtype=float64)
    H = zeros((nt, nx), dtype=float64)
    E = zeros(nt, dtype=float64)
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

    # From here we compute H(x,t). The boundary terms are messy, but most
    # of this corresponds exactly to the analytical form.
    # Sort out the t = 0 boundary.
    H[0][0] = 0.5*(((1/dt)*(U[1][0] - U[0][0]))**2) + 0.5*(((1/dx)*(U[0][1] - U[0][0]))**2)
    H[0][nx-1] = 0.5*(((1/dt)*(U[1][nx-1] - U[0][nx-1]))**2) + 0.5*(((1/dx)*(U[0][nx-1] - U[0][nx-2]))**2)
    for j in range(1,nx-1):
        H[0][j] = 0.5*(((1/(dt))*(U[1][j]-U[0][j]))**2) + 0.5*(((1/(2*dx))*(U[0][j+1]-U[0][j-1]))**2)
    E[0] += dx*sum(H[0])
    # Now start compuating H more generally.
    for n in range(1,nt-1):
        # Sort out the j boundaries.
        H[n][0] = 0.5*(((1/(2*dt))*(U[n+1][0] - U[n-1][0]))**2) + 0.5*(((1/dx)*(U[n][1] - U[n][0]))**2)
        H[n][nx-1] = 0.5*(((1/(2*dt))*(U[n+1][nx-1]-U[n-1][nx-1]))**2)+0.5*(((1/dx)*(U[n][nx-1]-U[n][nx-2]))**2)
        for j in range(1,nx-1):
            H[n][j] = 0.5*(((1/(2*dt))*(U[n+1][j]-U[n-1][j]))**2) + 0.5*(((1/(2*dx))*(U[n][j+1]-U[n][j-1]))**2)
        E[n] += dx*sum(H[n])
    # Sort the t_max boundaries.
    H[nt-1][0] = 0.5*(((1/dt)*(U[nt-1][0] - U[nt-2][0]))**2) + 0.5*(((1/dx)*(U[nt-1][1] - U[nt-1][0]))**2)
    H[nt-1][nx-1]=0.5*(((1/dt)*(U[nt-1][nx-1]-U[nt-2][nx-1]))**2)+0.5*(((1/dx)*(U[nt-1][nx-1]-U[nt-1][nx-2]))**2)
    for j in range(1,nx-1):
        H[nt-1][j]=0.5*(((1/(dt))*(U[nt-1][j]-U[nt-1][j]))**2)+0.5*(((1/(2*dx))*(U[nt-1][j+1]-U[nt-1][j-1]))**2)
    E[nt-1] += dx*sum(H[nt-1])

    # F is the relative error of the energy.
    F = array([(e/E[0] - 1) for e in E])
    F[-1] = F[-2] # Due to the boundaries the final term is not well constrained.
    return array([t,F])

# Try four different values of h.
a1 = run(0.1,[0,14],[-7,7])
a2 = run(0.05,[0,14],[-7,7])
a3 = run(0.025,[0,14],[-7,7])
a4 = run(0.0125,[0,14],[-7,7])

# Now we plot everything.
plot(a1[0],a1[1],'-',label="$h$ = 0.1")
plot(a2[0],a2[1],'-',label="$h$ = 0.05")
plot(a3[0],a3[1],'-',label="$h$ = 0.025")
plot(a4[0],a4[1],'-',label="$h$ = 0.0125")
xlabel("$t$")
ylabel("$E(t)$")
legend()
title("Relative error of energy for Dirichlet boundaries")
savefig("Ed.pdf",format="pdf")
show()

