#!/bin/bash

# This script is intended to run as cron job
# Run a default query to the hug server and check if it is running

report_to_master(){
  echo "Help Hug Server not running ?"
}
error=false
curl -H "Content-Type: application/json" -X POST --digest --data '{"statid":"10393","date":["20000231"],"variable": "temperature"}' -o download.zip http://localhost:80

if [ -e download.zip ]; then
  # might work, check inside
  # does not work yet, mail server rejects emails from us
  #
    unzip -l download.zip > report.txt
    if [ $? -eq 0 ]; then
      # unzip works
      echo "$(date) Hug Server Operational. Report back to you tomorrow" | sendmail "master@univie.ac.at"
    else
      error=true
    fi
fi

if $error; then
  # Error
  cat report.txt | sendmail "master@unvie.ac.at"
fi