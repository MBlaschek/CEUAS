import os,sys
import time
import datetime


ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print('Starting at: ', st)

os.system('/opt/anaconda3/bin/python3  build_311c_cdmfiles_ALL.py  -d all -o PROVA_TUTTI_logstationid ') 

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print ('Finished at: ', st )
