import datetime as dt
import os, glob, time, subprocess
import numpy as np
import pandas as pd
from io import StringIO

def ndays_in_month(datetime1):
  # return number of days in the month
  firstday_thismonth = dt.datetime(datetime1.year, datetime1.month, 1)
  oneday_nextmonth   = firstday_thismonth + dt.timedelta(days=32)
  firstday_nextmonth = dt.datetime(oneday_nextmonth.year, oneday_nextmonth.month, 1)
  ndays1             = (firstday_nextmonth - firstday_thismonth).days
  return ndays1

def open_clean_bufr_files(bufr_file_output_fmt0, yyyymmdd_dt0, nprocesses0, loop_hours_to_add=[0,6,12,18,24]):

  # create file array
  if type(nprocesses0) is int:
    filearray1 = {iproc:{} for iproc in range(1,nprocesses0+1)}
    procstars  = ['*']
  else:
    filearray1 = {iproc:{} for iproc in nprocesses0}
    procstars  = ['{0}'.format(iproc) for iproc in nprocesses0]
  # make sure the initial date is at 00 UTC
  yyyymmdd_dt_start = dt.datetime(yyyymmdd_dt0.year, yyyymmdd_dt0.month, yyyymmdd_dt0.day, 0, 0, 0)
  # consider the date of the data, plus the following day 00 UTC
  yyyymmdd_dts = [yyyymmdd_dt_start + dt.timedelta(hours=h1) for h1 in loop_hours_to_add]
  for yyyymmdd_dt1 in yyyymmdd_dts:
    for procstar1 in procstars:
      bufr_file_hh_ipstar = bufr_file_output_fmt0.format(yyyymmdd_dt1.strftime('%Y'), yyyymmdd_dt1.strftime('%m'), yyyymmdd_dt1.strftime('%d'), yyyymmdd_dt1.strftime('%H'), yyyymmdd_dt_start.strftime('%Y%m%d'), procstar1) # each file is named according to the data source date (filename) and the destination date (pathname), plus a processor number
      bufr_file_path = os.path.dirname(bufr_file_hh_ipstar)
      if not os.path.exists(bufr_file_path):
        # create directory
        os.makedirs(bufr_file_path)
      else:
        pre_existing_files = glob.glob(bufr_file_hh_ipstar)
        for fic1 in pre_existing_files:
          # remove existing files
          os.remove(fic1)
    # create empty files
    for iproc in filearray1.keys():
      bufr_file_hh_ip1 = bufr_file_output_fmt0.format(yyyymmdd_dt1.strftime('%Y'), yyyymmdd_dt1.strftime('%m'), yyyymmdd_dt1.strftime('%d'), yyyymmdd_dt1.strftime('%H'), yyyymmdd_dt_start.strftime('%Y%m%d'), iproc)
      filearray1[iproc][yyyymmdd_dt1.strftime('%Y%m%d%H')] = open( bufr_file_hh_ip1, 'wb' )

  # return array of files opened for binary writing
  return filearray1

def close_clean_bufr_files(filearray1, keep_first_proc_even_if_empty=True):

  for proc1 in filearray1.keys():
    for dt1 in filearray1[proc1].keys():
      mydatafile1 = filearray1[proc1][dt1]
      # close file
      mydatafile1_name = mydatafile1.name
      mydatafile1.close()
      # verify file contents: empty or not?
      dwait = 0
      while (not os.path.exists(mydatafile1_name)) and (dwait<10):
        # added a wait max. 10 seconds in case the filesystem isn't responsive immediately
        os.sleep(1)
        dwait += 1
      # we waited long enough... can't hold the whole machine for this!
      if os.path.exists(mydatafile1_name) and (os.path.getsize(mydatafile1_name)==0) and ((not keep_first_proc_even_if_empty) or (keep_first_proc_even_if_empty and (proc1>1))): # remove empty files, leaving only the first one to keep track processing was done, if needed
        os.remove(mydatafile1_name)

def get_synoptic_date_and_time(obs_dt1):
  obs_dt1_hhmmss = int(obs_dt1.strftime('%H%M%S'))
  if obs_dt1_hhmmss<=30000:
    synobs_dt1 = dt.datetime(obs_dt1.year, obs_dt1.month, obs_dt1.day,  0) # 00 UTC
  elif obs_dt1_hhmmss<=90000:
    synobs_dt1 = dt.datetime(obs_dt1.year, obs_dt1.month, obs_dt1.day,  6) # 06 UTC
  elif obs_dt1_hhmmss<=150000:
    synobs_dt1 = dt.datetime(obs_dt1.year, obs_dt1.month, obs_dt1.day, 12) # 12 UTC
  elif obs_dt1_hhmmss<=210000:
    synobs_dt1 = dt.datetime(obs_dt1.year, obs_dt1.month, obs_dt1.day, 18) # 18 UTC
  else:
    synobs_dt1 = dt.datetime(obs_dt1.year, obs_dt1.month, obs_dt1.day,  0) + dt.timedelta(days=1) # 00 UTC next day
  return synobs_dt1

def bufr_count(bufr_files):
  nmsgs = 0
  for bufr_file1 in bufr_files:
    nmsgs += int(subprocess.run(['bufr_count',bufr_file1], capture_output=True).stdout.strip().decode('utf-8'))
  return nmsgs

def bufr_ls(bufr_files,bufr_elements):

  HD = os.environ['SCRATCHDIR']
  big_df = None
  nbufr_elements=len(bufr_elements)
  empty_csv_line='thislineisempty'
  for bufr_file1 in bufr_files:
    fic_ls1 = HD +'/'+os.path.basename(bufr_file1)+'.list'
    if os.path.exists(fic_ls1):
      os.remove(fic_ls1)
    nmsgs = bufr_count([bufr_file1])
    if nmsgs==0:
      continue # skip this file
    os.system("bufr_ls -p {0} {1} > {2}".format(','.join(bufr_elements),bufr_file1,fic_ls1))

    readlines=open(fic_ls1,'r').readlines()

    def process_one_line(line1):
      if 'message' not in line1 and bufr_file1 not in line1 and 'localLatitude' not in line1:
        elms = line1.replace('\n','').split(' ')
        while '' in elms:
          elms.remove('')
        if len(elms)==nbufr_elements:
          return ','.join(elms)
      return empty_csv_line
    bigarray = np.vectorize(process_one_line)(np.array(readlines))

    wok=np.where(bigarray!=empty_csv_line)[0]
    csv_text="{0}\n".format(','.join(bufr_elements))
    csv_text += "\n".join(bigarray[wok])
    df1 = pd.read_csv(StringIO(csv_text),sep=",")
    if len(df1)!=nmsgs:
      print("ERROR... we did not read all the messages...",nmsgs,len(df1))
      stop
    if len(df1)>0:
      if big_df is None:
        big_df = df1.copy()
      else:
        big_df = pd.concat([ big_df, df1], axis=0)

  return big_df

def bufr_get(bufr_files,bufr_elements):
  import eccodes

  bufr_types_mapping={'localLatitude':float,'localLongitude':float,'ident':str,'typicalDate':int,'typicalTime':int}
  packed_elms=['pressureReducedToMeanSeaLevel','nonCoordinatePressure', 'oceanographicWaterTemperature', 'airTemperature', 'windDirection', 'windSpeed', 'significantWaveHeight']

  lunpack=False
  for x in packed_elms:
    if x in bufr_elements:
      bufr_types_mapping[x]=float
      lunpack=True
  bufr_types={x:bufr_types_mapping[x] for x in bufr_elements}
 
  if lunpack: 
    def msg_unpack(ibufr):
      eccodes.codes_set(ibufr,"unpack",1)
  else:
    def msg_unpack(ibufr):
      return

  def check_message(ibufr):
    msg_unpack(ibufr)
    nsubsets = eccodes.codes_get(ibufr,'numberOfSubsets')
    if nsubsets!=1:
      print("ERROR... numberOfSubsets is not 1...",nsubsets)
      stop
    values = []
    for bufr_element1 in bufr_elements:
      type1 = bufr_types[bufr_element1]
      try:
        values.append(eccodes.codes_get(ibufr,bufr_element1,ktype=type1))
      except:
        values.append({float:np.nan,int:None,str:''}[type1])
      if (type1 is float) and (~np.isnan(values[-1])) and (np.abs(values[-1])>1e11):
        values[-1]=np.nan
    return values

  nbufr_elements=len(bufr_elements)
  bufr_values_big = []
  for bufr_file1 in bufr_files:
    bufr_file_in1 = open(bufr_file1, 'rb')
    while 1==1:
      ibufr_in = eccodes.codes_bufr_new_from_file(bufr_file_in1)
      if ibufr_in is None: # last message, we are done
        break
      bufr_values = check_message(ibufr_in)
      del ibufr_in
      nbufr_values= len(bufr_values)
      if nbufr_values>0:
        if nbufr_values!=nbufr_elements:
          print("ERROR... number of values different from number of elements",nbufr_values,nbufr_elements)
          stop
      bufr_values_big.append(bufr_values)
  mydict = {}
  for ib1,bufr_element1 in enumerate(bufr_elements):
    myvals = [bufr_val1[ib1] for bufr_val1 in bufr_values_big]
    mydict[bufr_element1] = myvals
  big_df = pd.DataFrame(mydict)
  return big_df
 

