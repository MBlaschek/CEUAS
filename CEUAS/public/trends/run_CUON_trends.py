#!/usr/bin/env python
import os

input = "/mnt/users/scratch/leo/scratch/converted_v19/long/"
output = "./cuon_trends"
version = "2024-04-24"

for sd, ed in [("1994-01-01", "2023-12-31"), ("1973-01-01", "2002-12-31"), ("1958-01-01", "1987-12-31")]:
    for var in ['ta']:
        for adjustment in ['RAOBCORE_bias_estimate', 'RASE_bias_estimate', 'RICH_bias_estimate', 'RISE_bias_estimate']:
            pressure = "10000"
            cbl = "2"
        
            os.system("python CUON_trends.py -v "+ var +" -sd "+ sd +" -ed "+ ed +" -i "+ input +" -o "+ output +" -a "+ adjustment +" -p "+ pressure +" -cbl "+ cbl +" -vr "+ version)


for sd, ed in [("1940-01-01", "1959-12-31"), ("1994-01-01", "2023-12-31"), ("1973-01-01", "2002-12-31"), ("1958-01-01", "1987-12-31")]:
    for var in ['u', 'v', 'wd']:
        pressure = "70000"
        adjustment = "wind_bias_estimate"
        cbl = "2"
    
        os.system("python CUON_trends.py -v "+ var +" -sd "+ sd +" -ed "+ ed +" -i "+ input +" -o "+ output +" -a "+ adjustment +" -p "+ pressure +" -cbl "+ cbl +" -vr "+ version)

for sd, ed in [("1994-01-01", "2023-12-31"), ("1973-01-01", "2002-12-31"), ("1958-01-01", "1987-12-31")]:
    for var in ['rh', 'dp']:
        pressure = "10000"
        adjustment = "humidity_bias_estimate"
        if var == 'rh':
            cbl = "0.3"
        else:
            cbl = "10"
    
        os.system("python CUON_trends.py -v "+ var +" -sd "+ sd +" -ed "+ ed +" -i "+ input +" -o "+ output +" -a "+ adjustment +" -p "+ pressure +" -cbl "+ cbl +" -vr "+ version)
