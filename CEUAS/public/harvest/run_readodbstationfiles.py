import os,sys
import configparser

config = configparser.ConfigParser()
config.read('input/readodbstationfiles_input_test.ini')
githome = config['PATHS']['githome'] # path of the base git dir

""" Read the data to extract from the parameter file """
datasets         = config['DATA']  ['datasets']
station           = config['DATA']  ['station']
variables        = config['DATA']  ['variables']
outdir             = config['OUTPUT']['outdir']
databasepath = config['PATHS'] ['databasepath']
gribdir            = config['PATHS'] ['gribdir']


from try_read import *


""" Running the converter """
if __name__ == "__main__":

    for e in datasets.split(','):
        exp = e
        for v in variables.split(','):
            print('v is', v )
            run_converter(dataset=e, single_stat= station, pool=False, varno= int(v), debug = True )
            print('Finished with the database', e , ' **** for the variable: ', str(v))

    exit()

print("*** Finished converting to netCDF files !")
