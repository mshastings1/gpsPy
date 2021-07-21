#!/usr/bin/env python3.7

# add parent dir to python path
import sys
sys.path.append('../')

# import the gpsBruh class
from gpsBruh import gpsBruh 

# make instance of gpsBruh class from UNR file
g = gpsBruh('./MANA.CA.tenv3')

# subset dataframe from 2015 to 2021
g.time_window(2015,2021)

# read in the earthquakes file 
g.read_eqDates('./MANA.eqs')

# convert YYMMMDD to decimal year
g.get_eqTimes()

# make a plot of the NEU componenets and display it
fig = g.plotNEU()
fig.show()

# weighted least squares model the north, east, and up componenets (see README for formulation)
g.model_gps()

# compute the rms for the models and store them 
g.rms() 

# make a plot for the data and the models
fig2 = g.plotModel()
fig2.show()

# print out the dictionary holding the model parameters and fit
print("*** MODEL PARAMETERS AND FIT ***")       # note that the parameters are in meters
print(g.model)

# print out all the attributes and methods that can investigated on the gpsBruh instance
print("*** METHODS AND ATTRIBUTES ***")
print("*** SEE gpsBruh CLASS DEF FOR MORE DETAILS ***") 
print(dir(g))

# save the gpsBruh object to a file
g.save_pickleBruh()

