# Fast analyze lat mismatch


import pandas as pd
import matplotlib.pyplot as plt

import numpy as np
import cartopy.crs as ccrs


df = pd.read_csv('inventories/1759_minusSignLatMismatch.txt',  sep = "_" , names = ["file", "dist", "lat","lon","lat_w","lon_w"] )

# primary_id      file_lat        wban_lat        file_lon        wban_lon        file
old = pd.read_csv("inventories/OLD_era5_1759_WBAN_latitude_mismatch.dat" , sep = "\t" , index_col = False )
old_f = open( "inventories/OLD_era5_1759_WBAN_latitude_mismatch.dat",'r').readlines()
old_files = [f.split("\t")[-1].replace("\n",'') for f in old_f if 'file' not in f ]

miss = []
for f in old_files:
    if f not in df['file'].values:
        miss.append(f)
        
        
MISS = { "lat_file" : [] , "lon_file" : [], "lat_inv" : [], "lon_inv": [] }

# primary_id      file_lat        wban_lat        file_lon        wban_lon        file
for l in old_f :
    if "primary" in l:
        continue
    ll = l.split("\t")
    
    f = ll[6].replace("\n","")
    
    if f not in miss:
        continue 
    MISS["lat_file"].append(float(ll[1]))
    MISS["lat_inv"].append(float(ll[2]))
    
    MISS["lon_file"].append(float(ll[3]))    
    MISS["lon_inv"].append(float(ll[4]))
         
         
    
    


ax = plt.axes(projection=ccrs.PlateCarree())
ax.stock_img()
plt.title("Latitude mismatch old 109 stations - ERA5 1759")
plt.scatter(old["wban_lon"], old["wban_lat"],
         color='gold', marker='o', s=3,
         transform=ccrs.PlateCarree(), label='WBAN Coord',
         )

plt.scatter(old["file_lon"], old["file_lat"],
         color='lime', marker='o', s=3,
         transform=ccrs.PlateCarree(), label='File Coord',
         )
plt.legend()
plt.show()
plt.close()
plt.savefig("All_Old.png")



ax = plt.axes(projection=ccrs.PlateCarree())
ax.stock_img()
plt.title("Latitude mismatch MISSED stations - ERA5 1759")


plt.scatter(MISS["lon_file"], MISS["lat_file"],
         color='orange', marker='o', s=3,  label='File Coord',
         transform=ccrs.PlateCarree(),
         )

plt.scatter( MISS["lon_inv"], MISS["lat_inv"], 
         color='red', marker='o', s=3,
         transform=ccrs.PlateCarree(), label='WBAN Coord',
         )
plt.legend()
plt.show()
plt.close()
plt.savefig("Missed.png")


"""
plt.scatter(df["lon_w"], df["lat_w"],
         color='red', marker='o', s=2,
         transform=ccrs.PlateCarree(),
         )

plt.scatter(df["lon"], df["lat"],
         color='magenta', marker='o', s=2,
         transform=ccrs.PlateCarree(),
         )
"""


print(0)