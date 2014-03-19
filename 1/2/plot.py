#!/usr/bin/env python3
from matplotlib.pyplot import *
from numpy import *

a = genfromtxt("ascii.dat").transpose()

semilogy(a[0],a[1],'.')
xlabel("$h$")
ylabel("$\epsilon(h)$")
savefig("epsilon.pdf",format="pdf")
show()
