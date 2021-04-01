#!/usr/bin/env python
# By Michael B
# Purpose: Process all CEUAS files in parallel
import sys
import time
import os
import logging

# def job(ifile):
#     #square = my_number * my_number
#     time.sleep(1)
#     return ifile.split('/')[-1]

def job(ifile):
    """Run adjustment procedure and save into backend file + raobcore format to ftp
    """
    import cds_eua3 as eua
    import raso_adj_cdm_v1 as adj
    import xarray as xr

    # Activate logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s | %(funcName)s - %(levelname)s - %(message)s')
    # logger.handlers = [] # reset remove all handlers
    ch = logging.FileHandler('logs/{}.log'.format(ifile.split('/')[-1]))
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info("Logging %s", ifile)
    # needs to be before 
    status = ""
    try:
        if False:
            # DEBUG
            iofile = eua.CDMDataset(ifile)
            iofile.reopen(write_to_filename='monkey.nc', mode='r+', strict=True)
        else:
            iofile = eua.CDMDataset(ifile)
            newfilename = '/raid60/scratch/mblaschek/CEUAS_db/' + ifile.split('/')[-1]
            # copy only parts
            iofile.copy(newfilename, force=True, 
                        variables=['/observations_table/observation_value',
                                   '/observations_table/date_time',
                                   '/observations_table/z_coordinate', 
                                   '/observations_table/observed_variable',
                                   '/era5fb/fg_depar@body'], 
                        groups=['advanced_homogenisation', 'recordindices', 'header_table'])
            #iofile.reopen(write_to_filename='/raid60/scratch/mblaschek/CEUAS_db/' + ifile.split('/')[-1], mode='r+', strict=True)  # is that ok? or copy first?
            iofile.close()
            iofile = eua.CDMDataset(newfilename, mode='r+')
        status += "R+"
        variable = 'relative_humidity'
        # Code 34, dew point departure
        # variable = 'dew_point_departure'
        data = iofile.read_data_to_cube(variable,
                                        feedback='fg_depar@body',
                                        feedback_group='era5fb',
                                        )
        status += "D"
        # should contain variables
        # e.g. temperature, temperature_an_depar
        variable, depar = list(data.keys())
        data = xr.Dataset(data)
        # run adjustment procedure
        data = adj.adjustment_procedure(data,
                                    obs_name=variable,
                                    dep_name=depar,
                                    metadata=False,
                                    times=[0, 12],
                                    dim='time',
                                    plev='plev',
                                    return_dataset=False,
                                    )
        status += "A"
        # TODO Convert adjustments to other variables?
        #
        #
        # Write back adjusted (interpolation, extrapolation)
        #
        iofile.write_observed_data('humidity_bias_estimate',
                                   varnum=eua.cdm_codes[variable],
                                   cube=data['adjustments'],
                                   group='advanced_homogenisation',
                                   attributes={'version':'1.0'}
                                   )
        status += "WB"
        #
        # Write as well to RAOBCORE FORMAT
        #
        ftppath = '/raid/home/srvx7/ftpsecure/'
        try:
            os.remove(ftppath + 'pub/RAOBCOREv1.8/humidity/' + iofile.name)
        except:
            pass
        iofile.convert_to_raobcore(variable, 
                                   ftppath + 'pub/RAOBCOREv1.8/humidity/' + iofile.name, 
                                   feedback=['humidity_bias_estimate', 'fg_depar@body'],
                                   feedback_group=['era5fb', 'advanced_homogenisation'],
                                   source='RAOBCORE v1.8, relative humidity', 
                                   title='Station daily relative humidity series with ERA5 background departure statistics and bias estimates')
        status += "CR"
    except Exception as e:
        logger.error(e)
    return status
    
if __name__ == "__main__":
    import argparse
    from tqdm.contrib.concurrent import process_map  # or thread_map
    
    parser = argparse.ArgumentParser(description="Run standardized radiosonde homogenisation software")
    parser.add_argument("-n", "--nprocs", type=int, help="Number of processes")
    parser.add_argument("files", nargs='*', help="Input files")
    args, unknown = parser.parse_known_args()
    if args.files is None:
        print("Missing files")
        parser.print_help()
        sys.exit(1)
    
    n = len(args.files)
    print(n, args.nprocs)
    
    if True:
        r = process_map(job, args.files, max_workers=args.nprocs, chunksize=1)
    else:
        print(job(args.files[0]))
    
    with open('adjustments.results.list', 'w') as f:
        for i,j in zip(r, args.files):
            f.write("{} : {}\n".format(i,j))
    #with Pool(args.nprocs) as p:
    #    r = list(tqdm.tqdm(p.imap(job, range(n)), total=n))    
    