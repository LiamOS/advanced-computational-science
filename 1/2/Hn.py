#!/usr/bin/env python3
# Program to solve the wave equation with c=1 and an external potential.
# Neumann boundaries are used, and imposed explicitly.
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
      U[n][j] = c*(U[n-1][j-1] - 2*U[n-1][j] + U[n-1][j+1]) - U[n-2][j] +2*U[n-1][j]-(dt**2)*V[j]*U[n-1][j]
    U[n][0] = U[n][1]
    U[n][nx-1] = U[n][nx-2]

  # Sort out the t = 0 boundary.
  H[0][0] = 0.5*(((1/dt)*(U[1][0] - U[0][0]))**2) + 0.5*(((1/dx)*(U[0][1] - U[0][0]))**2)+0.5*(V[0]*(U[0][0]**2))
  H[0][nx-1] = 0.5*(((1/dt)*(U[1][nx-1] - U[0][nx-1]))**2) + 0.5*(((1/dx)*(U[0][nx-1] - U[0][nx-2]))**2)+0.5*(V[nx-1]*(U[0][nx-1]**2))
  for j in range(1,nx-1):
    H[0][j]=0.5*(((1/(dt))*(U[1][j]-U[0][j]))**2)+0.5*(((1/(2*dx))*(U[0][j+1]-U[0][j-1]))**2)+0.5*(V[j]*(U[0][j]**2))
  # Now start compuating H more generally.
  for n in range(1,nt-1):
    # Sort out the j boundaries.
    H[n][0] = 0.5*(((1/(2*dt))*(U[n+1][0] - U[n-1][0]))**2) + 0.5*(((1/dx)*(U[n][1] - U[n][0]))**2)+0.5*(V[0]*(U[n][0]**2))
    H[n][nx-1] = 0.5*(((1/(2*dt))*(U[n+1][nx-1]-U[n-1][nx-1]))**2)+0.5*(((1/dx)*(U[n][nx-1]-U[n][nx-2]))**2)+0.5*(V[nx-1]*(U[n][nx-1]**2))
    for j in range(1,nx-1):
      H[n][j]=0.5*(((1/(2*dt))*(U[n+1][j]-U[n-1][j]))**2)+0.5*(((1/(2*dx))*(U[n][j+1]-U[n][j-1]))**2)+0.5*(V[j]*(U[n][j]**2))
  for n in range(0,nt):
    E[n] = dx*sum(H[n])
  E[nt-1] = E[nt-2]

  F = array([(e/E[0] - 1) for e in E])
  return array([t,F])

a2 = run(0.05,[0,20],[-20,20])
a1 = run(0.025,[0,20],[-20,20])

plot(a2[0],a2[1],'-',label="$h$ = 0.05")
plot(a1[0],a1[1],'-',label="$h$ = 0.025")
xlabel("$t$")
ylabel("$E(t)$")
legend()
savefig("En.pdf",format="pdf")
show()

