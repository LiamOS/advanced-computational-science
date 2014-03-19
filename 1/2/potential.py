#!/usr/bin/env python3
# Program to solve the wave equation with c=1 and external potential using
# central differences. Neumann boundaries are used.
# Liam O'Sullivan
from matplotlib.pyplot import *
from numpy import *
#sigma for gaussian
sig = 1
# Define our space and time step sizes.
dx = 0.05
dt = dx/sqrt(2)
# Define our Temporal and Spatial lengths.
Time = [0,20]
Dist = [-20,20]

# Define variables with the number of space/time steps for ease later.
# They must be cast as ints for use with range().
nt = int((Time[1] - Time[0])/dt) + 1
nx = int((Dist[1] - Dist[0])/dx) + 1
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
    U[0][j] = e**(-((x[j]-6)**2)/(2*(sig**2)))
# Initialise our potential array.
for j in range(0,nx):
    V[j] = 1/(cosh(x[j])**2)
# Plot the potential with a broken line for visual purposes.
plot(x,V,'--',color='grey')
# c is not c from the wave equation, but nu^2.
c = (dt/dx)**2
# I manually impose the initial derivative condition, by computing the second
# time step as d/dx of the analytic u_0, using the Euler method(Forgive me).
for j in range(1,nx-1):
    U[1][j] = U[0][j] + dt*(-((x[j]-6)/(sig**2))*(e**(-((x[j]-6)**2)/(2*(sig**2)))))

# Now perform our integration.
for n in range(2,nt):
    for j in range(1,nx-1):
        U[n][j] = c*(U[n-1][j-1] - 2*U[n-1][j] + U[n-1][j+1]) - U[n-2][j] + 2*U[n-1][j] - (dt**2)*V[j]*U[n-1][j]
        # Impose the Neumann condition here:
    U[n][0] = U[n][1]
    U[n][nx-1] = U[n][nx-2]
    if ((n-2)%80) is 0:
        plot(x,U[n],label="t = "+str(round(Time[1]*(n/nt),1)))
xlabel("$x$")
ylabel("$u(x,t)$")
xlim([-20,20])
ylim([-0.75,1])
legend()
savefig("cosh.pdf",format="pdf")
show()
