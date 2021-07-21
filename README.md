# gpsPy

gpsPy is a small repository I made for working with data from the University of Nevada at Reno's MAGNET Program. 
Currently the core module is the gpsBruh.py file which contains a self-title class "gpsBruh". Don't ask why
I named it that, or why there are gibberish print statements throughout the methods. I was trying to make it 
entertaining for myself when I built it and did not intend to put it out to the public.

I developed this to make a better way of reading in the UNR data and storing it in a way that is easily readable
and loadable to the python workspace. Check out the gpsBruh class description for more information.

## Dependencies

There are not many dependecies as I like to create my own logic rather than relying on tons of other peoples modules.
The pandas dependency is to make use of the DataFrame object because I like how much more readable they are than 
2D arrays.

pandas v1.2.3 <br>
numpy v1.19.2 <br>
pickle-mixin v1.0.2 <br>
matplotlib v3.3.4 <br>

## The gpsBruh Class

This class was designed to hold the time series information as well as metadata associated with a gps station.
That includes the north, east, and up ocmponent time series as well as information like the station name, 
station locations, earthquake times that occur within the series, and more. Additionally I also added methods
to be able to model the time series as a linear model with annual and semiannual signals and heaviside functions
to represent earthquakes within the series.

<img src="https://render.githubusercontent.com/render/math?math=y = a %2B bt %2B csin(\omega t) %2B dcos(\omega t) %2B esin(2\omega t) %2B fcos(2\omega t) %2B H(t_{eq})">

where the <img src="https://render.githubusercontent.com/render/math?math=H(t_{eq})"> is a heaviside function for
each of the earthquakes in the time series. The earthquake lists can be copied from the station pages on the UNR site. 
The dates are in YYMMMDD format similar to how they are in the data and on the UNR site and the model adds step functions for 
every earthquake in the list.

I will add more later but for now see the examples if you are want practical implications.

## Supporting Files

getUNR_CA.sh - bash script for downloading data from the stn_CA.list file (note CA is for stable Caribbean plate) <br>
stn_CA.list - file with station names to download. <br>
exanples/BruhDotMP4.py - implementation of the gpsBruh class to read in a gps time series and make a plot of the data over a certain time frame <br>
examples/MANA.CA.tenv3 - data file downloaded from UNR MAGNET Program <br>
examples/MANA.eqs - earthquakes list from the MANA station on UNR page <br>
examples/MANA.pkl - pickle file of the gpsBruh object saved from BruhDotMP4.py <br>
examples/CowPeople.py - read in the MANA.pkl file to show its the same as the object from BruhDotMP4.py <br>

