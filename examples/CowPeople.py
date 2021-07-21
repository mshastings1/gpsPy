#!/usr/bin/env python3.7

# add parent dir to python path
import sys
sys.path.append('../')

# import the gpsBruh class
from gpsBruh import gpsBruh, load_pickleBruh

# load the pickle file and plot the NEU comps
g = load_pickleBruh('./MANA.pkl')
fig = g.plotNEU()
fig.show()
