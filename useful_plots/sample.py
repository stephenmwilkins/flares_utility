

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import cmasher as cmr

cm = cmr.torch

import flare.plt as fplt


fig, ax = fplt.simple()


ax.set_xlabel(r'$\rm x/units$')
ax.set_ylabel(r'$\rm y/units$')

N = 200

X = np.random.randn(N)
Y = np.random.randn(N)
s = 20*np.random.random(N)
c = np.sqrt(X**2+Y**2)

ax.scatter(X,Y,s=s,c=c,cmap=cm)

mx = 3

ax.set_xlim([-3,3])
ax.set_ylim([-3,3])

fig.savefig(f'figures/sample.pdf')

fig.clf()
