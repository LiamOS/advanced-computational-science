#!/usr/bin/python

from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from numpy import *

LU = 0
w = pi/5
dx = 0.1
x = arange(-5,5,dx)
y = arange(-5,5,dx)
X, Y = meshgrid(x,y)
nx = len(X)
U = zeros((nx,nx))
O = zeros((nx,nx))

for i in range(0,nx):
    a = sin(w*x[i])
    U[0][i]    = a
    U[nx-1][i] = -a
    U[i][0]    = -a
    U[i][nx-1] = a

fig = figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X,Y,U,cmap=cm.coolwarm)
show(block=False)
while True: #for n in range(0,200):
    LO = LU
    LU = 0
    for i in range(1,nx-1):
        for j in range(1,nx-1):
            U[i][j] = 0.25*(U[i-1][j] + U[i+1][j] + U[i][j-1] + U[i][j+1])
            LU += (U[i][j]**2)*(dx**2)
    dL = (LU - LO)/LU
    print("Normalised L2-norm = "+str(round(dL,5)))
    ax.cla()
    surf = ax.plot_surface(X,Y,U,cmap=cm.coolwarm)
    draw()
    #if abs(dL) < 1e-5: break

