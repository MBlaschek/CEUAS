import pyautogui
import time
from PIL import Image

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
time.sleep(5)
offset = 19
scroll_tick = 0
scroll_need = 0
scroll_limit = 6
for y in range(150,(150+10*offset),offset):
    for x in [80,120,160,200,250,280,320,360,420,480]:
        y_new = y - scroll_tick
        pyautogui.click(x, y_new)
        pyautogui.click(40, 100)
        #pyautogui.click(540, 70)
        for zoom in range(0,7):
            pyautogui.rightClick(1738, 627) # unzoom
            time.sleep(0.5)
            pyautogui.click(1780, 687) # unzoom
            if zoom == 0:
                pyautogui.click(800, 122) # change to zoom steps
            else:
                pyautogui.click(800, 122+180*zoom-40) # change to zoom steps
            time.sleep(1)
            pyautogui.dragTo(850, 122+180+180*zoom, 1, button="left") # make metadata larger
            time.sleep(1)
            pyautogui.screenshot("./sc_"+name[x]+"_"+str(int((y-130)/offset))+"_"+str(zoom)+"_.png")
            # unzoom
        # time.sleep(1)
        #pyautogui.click(540, 70)
        pyautogui.click(1028, 590) # if error message
        pyautogui.click(1285, 735) # if error message
        pyautogui.click(130, 100)
    scroll_need += 1
    if (scroll_need == 2) and (scroll_limit > 0):
        scroll_need = 0
        scroll_limit -= 1
        scroll_tick += offset * 2
        pyautogui.click(600, 1385)
