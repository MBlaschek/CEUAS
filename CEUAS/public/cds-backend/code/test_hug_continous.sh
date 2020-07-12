#!/bin/bash

# This script is intended to run as cron job
# Run a default query to the hug server and check if it is running

curl -H "Content-Type: application/json" -X POST --digest --data '{"statid":"10393","date":["20000231"],"variable": "temperature"}' -o download.zip http://localhost:80
