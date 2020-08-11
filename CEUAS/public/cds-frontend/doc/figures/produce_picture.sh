#!/bin/bash

python C3S_webpage_figure.py
# using graphics magick or imagemagick
gm convert -crop 600x1280+100 thetropicalu.jpg C3S_webpage_ballon.png
gm convert C3S_webpage_ballon.png C3S_webpage_logo.png +append ../overview.png

