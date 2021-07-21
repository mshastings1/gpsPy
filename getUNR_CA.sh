#!/bin/bash  
fstn='stn_CA.list' 
var=$(<${fstn})

for stni in ${var}
do
    wget -c http://geodesy.unr.edu/gps_timeseries/tenv3/plates/CA/${stni}.CA.tenv3
done