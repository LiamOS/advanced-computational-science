#!/usr/bin/python3
# Program for exercise 2.1. Liam O'Sullivan --- 10309537

from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from numpy import *

# Some obvious globals.
dt = 1/10000.
dx = 1/200.
dy = dx
Time = 0.1

# X,Y are the coordinates for the grid.
X = arange(-1,1,dx)
Y = arange(-1,1,dy)
lX, lY = meshgrid(X,Y) # This is used for the 3d plots.
nx = len(X)
ny = len(Y)
nt = int(Time/dt)
# Define our U grid. 
U = zeros((nt,nx,ny))
Uhalf = zeros((nx,ny)) # Define one for Uhalf that gets reused.

# As it's asked, print here the memory used by U, in GB.
print("Memory footprint of U = "+str(round(8*U.size/(1024*1024*1024),5))+"GB.")

# Set the initial condition, leaving the bounaries at 0.
for i in range(1,nx-1):
    for j in range(1,ny-1):
        U[0][i][j] = (1 - X[i]**2)*(1 - Y[j]**4)

# Display the initial surface.
fig = figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(lX,lY,U[0],cmap=cm.coolwarm)
title("t = 0")
xlabel("y")
ylabel("x")
savefig("1/t0.00.pdf",format='PDF')
show()

# Equivalent to \alpha_x/y in the notes, with them divided by two for good measure.
cx = dt/(dx**2)
cy = dt/(dy**2)
cx2 = cx/2.
cy2 = cy/2.

# This looks a bit abusive, but it's not that bad. These are the upper, lower, and diagonal
# parts of the matrix for the Thomson algorithm. There's probably a clever way to do this
# with a function and save some memory, but function calls actually slow it down, and it's not
# like 1d arrays are our biggest concern.
Aux = zeros(nx)
for i in range(0,nx):   Aux[i] = -cx2
Adx = zeros(nx)
for i in range(0,nx):   Adx[i] = 1 + cx
Alx = zeros(nx)
for i in range(1,nx-1): Alx[i] = -cx2
Auy = zeros(ny)
for i in range(0,ny-1): Auy[i] = -cy2
Ady = zeros(ny)
for i in range(0,ny):   Ady[i] = 1 + cy
Aly = zeros(ny)
for i in range(1,ny):   Aly[i] = -cy2

# These are also for the Thomson algorithm, being the array with the explicit part.
Dx = zeros((nx,ny))
Dy = zeros((nx,ny))

def Tsolve(a, b, c, d):
    ''' This funtion solves the Thomson algorithm, basically. '''
    ntmp = len(d)
    a[0] = 0 # Make sure.
    c[ntmp-1] = 0 # Make sure here, too.
    e = zeros(ntmp) # Allocate arrays for e, f.
    f = e.copy()
    ret = e.copy() # The array to return.
    e[0] = c[0]/b[0]
    f[0] = d[0]/b[0]
    for i in range(1,ntmp): # Sweep up and set e, f.
        e[i] = c[i]/(b[i] - a[i]*e[i-1])
        f[i] = (d[i] - a[i]*f[i-1])/(b[i] - a[i]*e[i-1])
    ret[ntmp-1] = f[ntmp-1]
    for i in reversed(range(0,ntmp-1)): # Sweep back down.
        ret[i] = f[i] - e[i]*ret[i+1]
    ret[1] = 0; ret[ntmp-2] = 0; # This is necessary for some reason.
    return ret


# Now we start the big loop over time.
print("Total n = "+str(nt))
for n in range(1,nt):
    print("n = " + str(n)) # Printing the run number is handy, allowing me to play 2048 in the meantime.

    # Sweep 
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            Dy[i][j] = cy2*U[n-1][i][j-1] + (1-cy)*U[n-1][i][j] + cy2*U[n-1][i][j+1]
        TMP = Tsolve(Aly, Ady, Auy, Dy[i])
        for k in range(1,len(TMP-1)):
            Uhalf[i][k] = TMP[k]

    # This loop is basically the same as above, but with a load of indices flipped about.
    for j in range(1,ny-1):
        for i in range(1,nx-1):
            Dx[j][i] = cx2*Uhalf[i-1][j] + (1-cy)*Uhalf[i][j] + cx2*Uhalf[i+1][j]
        TMP = Tsolve(Alx, Adx, Aux, Dx[j])
        for k in range(1,len(TMP)-1):
           U[n][k][j] = TMP[k]

    # Print the max and max location of the function. This was unspeakably useful for debugging...
    # which may have taken a while. 2d arrays in Python are a little funny sometimes, so the syntax
    # in this is pretty mad.
    print(unravel_index(U[n].argmax(), U[n].shape), max(U[n].flat))

    if (n*dt % 0.005 == 0):
        fig = figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(lX,lY,U[n],cmap=cm.coolwarm)
        xlabel("y"); ylabel("x") # The graph showed me I got the axis labeling wrong. This is easier than fixing it `properly`.
        title("t = "+str(n*dt))
        savefig("1/t"+str(round(n*dt,5))+".pdf",format='PDF')
        show()


# vim: ts=4
# vim: expandtab
