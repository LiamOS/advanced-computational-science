#!/usr/bin/python3
# Program for exercise 2.4. Liam O'Sullivan --- 10309537

# Although the question suggests keeping 'snapshots', I did not find this to be necessary.
# This program uses only one matrix for U, time independant, and overwrites when necessary.
# After each integration, the L^2 norm is computed and put into an array for later plotting.

from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from numpy import *
from time import clock

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
    ret[1] = 0; ret[ntmp-2] = 0 # This is necessary for some reason.
    return ret



dt = 1/1000.
Time = 0.1
nt = int(Time/dt)
time = arange(0,Time,dt)

def run(h):
    inTime = clock()
    dx = float(h)
    dy = dx

    # X,Y are the coordinates for the grid.
    X = arange(-1,1,dx)
    Y = arange(-1,1,dy)
    lX, lY = meshgrid(X,Y) # This is used for the 3d plots.
    nx = len(X)
    ny = len(Y)
    # Define our U grid. 
    U = zeros((nx,ny))
    Uhalf = zeros((nx,ny)) # Define one for Uhalf that gets reused.
    L2norm = zeros(nt) # Store our L^2 norms

    # Set the initial condition, leaving the bounaries at 0.
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            U[i][j] = (1 - X[i]**2)*(1 - Y[j]**4)

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


    L2 = 0 # Compute the first L^2 norm.
    for i in range(0,nx):
        for j in range(0,ny):
            L2 += dx*dy*U[i][j]**2
    L2norm[0] = sqrt(L2)

    # Now we start the big loop over time.
    print("Total n = "+str(nt))
    for n in range(1,nt):
        print("n = " + str(n)) # Printing the run number is handy, allowing me to play 2048 in the meantime.

        # Sweep 
        for i in range(1,nx-1):
            for j in range(1,ny-1):
                Dy[i][j] = cy2*U[i][j-1] + (1-cy)*U[i][j] + cy2*U[i][j+1]
            TMP = Tsolve(Aly, Ady, Auy, Dy[i])
            for k in range(1,len(TMP-1)):
                Uhalf[i][k] = TMP[k]

        # This loop is basically the same as above, but with a load of indices flipped about.
        for j in range(1,ny-1):
            for i in range(1,nx-1):
                Dx[j][i] = cx2*Uhalf[i-1][j] + (1-cy)*Uhalf[i][j] + cx2*Uhalf[i+1][j]
            TMP = Tsolve(Alx, Adx, Aux, Dx[j])
            for k in range(1,len(TMP)-1):
               U[k][j] = TMP[k]

        # Compute the L^2 norms with naive integration.
        L2 = 0
        for i in range(0,nx):
            for j in range(0,ny):
                L2 += dx*dy*U[i][j]**2
        L2norm[n] = sqrt(L2)

    outTime = clock()
    runTime = outTime - inTime
    print("Run Time = "+str(runTime))
    return L2norm, runTime


hs = [1/200.,1/50.,1/100.,1/300.,1/400.,1/500.,1/600.,1/700.,1/800.]
Ls = []
Ts = []
L = []

o = 0
for i in hs:
    a, ti = run(i)
    Ts.append(ti)
    Ls.append(a)
    L.append(1 - Ls[-1]/Ls[0])
    for w in [1,2,4,8]:
        if (o == w):
            plot(time,L[-1],label="h_"+str(o)+" = 1/"+str(int(1/i)))
    o += 1

xlabel("Time")
ylabel(r"$e(t)$")
legend(loc='best')
savefig("4/plot.pdf",format='PDF')
show()
clf()
plot(hs,Ts,'o')
ylabel("Runtime [s]")
xlabel("Resolution")
savefig("4/resolution.pdf",format='PDF')
show()

# vim: ts=4
# vim: expandtab
