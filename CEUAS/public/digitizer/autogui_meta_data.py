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

offset = 19 # pixel between the buttons vertically
scroll_tick = 0
scroll_need = 0
scroll_limit = 6
for y in range((150),(150+69*offset),offset):
    pyautogui.click(130, 100) # clear
    for x in [80,120,160,200,250,280,320,360,420,480]: # select all the boxes
        y_new = y - scroll_tick
        pyautogui.click(x, y_new)
    pyautogui.click(40, 100) # load
    pyautogui.rightClick(857, 193) # right click screen
    pyautogui.click(950, 395) # click metadata
    pyautogui.click(1100, 395) # click as text
    pyautogui.click(1670, 900) # to metadata right edge
    time.sleep(1)
    pyautogui.dragTo(1900, 900, 1, button="left") # make metadata larger
    pyautogui.screenshot("./meta_"+str((y-150)/offset)+".png")
    time.sleep(1)
    pyautogui.click(1028, 590) # if error message
    pyautogui.click(1285, 735) # if error message
    pyautogui.click(1875, 485) # close metadata again
    scroll_need += 1
    if (scroll_need == 2) and (scroll_limit > 0):
        scroll_need = 0
        scroll_limit -= 1
        scroll_tick += offset * 2
        pyautogui.click(600, 1385)
