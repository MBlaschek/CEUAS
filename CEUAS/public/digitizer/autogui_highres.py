import pyautogui
import os
import glob
import time
from PIL import Image
'''
This script is built to automatically take screen shots only of the WView data displaying tool.
All the pixel values are optimized for a resolution of 2560 x 1440.
The script will only work, if WView is in maximized to fit the full screen and if a 
dataset as well as a variable are already selected.
'''

# Mapping the sonde manufacturer to the pixel values in the selection columns. 
name = {
    80:"Vaisala",
    120:"SRS",
    160:"Meisei",
    200:"Modem",
    250:"Sip",
    280:"Graw",
    320:"3therm",
    360:"MKII",
    420:"Vais_GPS",
    480:"Graw_GPS",
}

# Giving the user a moment to switch to the WView window, before the script starts.
time.sleep(5)

# Keeping track of the current screenshot and how far the curser needs to move.
offset = 19
scroll_tick = 0
scroll_need = 0
scroll_limit = 6

# Iterate through the rows:
for y in range(160,(150+68*offset),offset):
    # Iterate through the columns:
    for x in [80,120,160,200,250,280,320,360,420,480]:
        print(len(glob.glob("./new_screen_shots/sc_"+name[x]+"_"+str(int((y-130)/offset))+"_*_.png")))
        
        # Skip already done screen shots:
        if len(glob.glob("./new_screen_shots/sc_"+name[x]+"_"+str(int((y-130)/offset))+"_*_.png")) >= 1 :
            continue
        
        # Select ascent and load it into viewer:
        y_new = y - scroll_tick
        pyautogui.click(x, y_new)
        pyautogui.click(40, 100)
        
        # The script will zoom into the plot to get more information while keeping the same resolution. 
        for zoom in range(0,7):
            pyautogui.rightClick(1738, 627) # unzoom
            time.sleep(0.5)
            pyautogui.click(1780, 687) # unzoom
            if zoom == 0:
                pyautogui.click(800, 126) # change to zoom steps
            else:
                pyautogui.click(800, 126+180*zoom-40) # change to zoom steps
            time.sleep(1)
            pyautogui.dragTo(850, 126+180+180*zoom, 1, button="left") # make metadata larger
            time.sleep(1)
            pyautogui.screenshot("./new_screen_shots/sc_"+name[x]+"_"+str(int((y-130)/offset))+"_"+str(zoom)+"_.png")
        # Give the system time to compute and save the screenshot:
        time.sleep(0.5)

        # Clear screen for next graph:
        pyautogui.click(1028, 590) # if error message
        pyautogui.click(1285, 735) # if error message
        pyautogui.click(130, 100)
    
    # Keep track of the rows done, so the scrowling is done correctly.
    scroll_need += 1
    if (scroll_need == 2) and (scroll_limit > 0):
        scroll_need = 0
        scroll_limit -= 1
        scroll_tick += offset * 2
        pyautogui.click(600, 1385)
