#!/usr/bin/python3

from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from numpy import *

dt = 1/800.
dx = 1/200.
dy = dx
Time = 0.1

X = arange(-1,1,dx)
Y = arange(-1,1,dy)
lX, lY = meshgrid(X,Y) # This is used for the 3d plots.
nx = len(X)
ny = len(Y)
nt = int(Time/dt)
U = zeros((nt,nx,ny))
Uhalf = zeros((nx,ny))

fig = figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(0,nx):
    for j in range(0,ny):
        U[0][i][j] = (1 - X[i]**2)*(1 - Y[j]**4)

ax.plot_surface(lX,lY,U[0],cmap=cm.coolwarm)
title("t = 0")
xlabel("x")
ylabel("y")
show()

# Equivalent to \alpha_x/y in the notes.
cx = dt/(dx**2)
cy = dt/(dy**2)
cx2 = cx/2.
cy2 = cy/2.

Aux = zeros(nx)
for i in range(0,nx-1):
    Aux[i] = -cx2
Adx = zeros(nx)
for i in range(0,nx):
    Adx[i] = 1 + 1*cx ; Adx[0] = 1.; Adx[nx-1] = 1.
Alx = zeros(nx)
for i in range(1,nx):
    Alx[i] = -cx2
Auy = zeros(ny)
for i in range(0,ny-1):
    Auy[i] = -cy2
Ady = zeros(ny)
for i in range(0,ny):
    Ady[i] = 1 + 1*cy ; Ady[0] = 1.; Ady[ny-1] = 1.
Aly = zeros(ny)
for i in range(1,ny):
    Aly[i] = -cy2

Dx = zeros((nx,ny))
Dy = zeros((nx,ny))

def Tsolve(a, b, c, d):
    ntmp = len(d)
    a[0] = 0 # Make sure.
    c[ntmp-1] = 0 # Make sure here, too.
    e = zeros(ntmp)
    f = e.copy()
    ret = e.copy()
    e[0] = c[0]/b[0]
    f[0] = d[0]/b[0]
    for i in range(1,ntmp):
        e[i] = c[i]/(b[i] - a[i]*e[i-1])
        f[i] = (d[i] - a[i]*f[i-1])/(b[i] - a[i]*e[i-1])
    ret[ntmp-1] = f[ntmp-1]
    for i in reversed(range(0,ntmp-1)):
        ret[i] = f[i] - e[i]*ret[i+1]
    return ret


for n in range(1,nt):
    print("n = " + str(n))
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            Dy[i][j] = U[n-1][i][j] + cy2*U[n-1][i][j-1] - cy*U[n-1][i][j] + cy2*U[n-1][i][j+1]
        TMP = Tsolve(Aly, Ady, Auy, Dy[i])
#        print(max(TMP))
        for k in range(1,len(TMP-1)):
            Uhalf[i][k] = TMP[k]

    for k in range(0,ny):
        Uhalf[nx-1][k] = 0; Uhalf[0][k] = 0
    for k in range(0,nx):
        Uhalf[k][ny-1] = 0; Uhalf[k][0] = 0

    for i in range(1,nx-1):
        for j in range(1,ny-1):
            Dx[i][j] = Uhalf[i][j] + cx2*Uhalf[i-1][j] - cy*Uhalf[i][j] + cx2*Uhalf[i+1][j]
        TMP = Tsolve(Alx, Adx, Aux, Dx[i])
#        print(max(TMP))
        for k in range(1,len(TMP-1)):
            U[n][i][k] = TMP[k]

    for k in range(0,ny):
        U[n][nx-1][k] = 0; U[n][0][k] = 0
    for k in range(0,nx):
        U[n][k][ny-1] = 0; U[n][k][0] = 0

    fig = figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(lX,lY,U[n],cmap=cm.coolwarm)
    xlabel("x")
    ylabel("y")
    show()


# vim: ts=4
# vim: expandtab
