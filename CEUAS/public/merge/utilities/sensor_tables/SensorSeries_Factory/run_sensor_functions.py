import sys
sys.path.append('modules')

from sensor_functions import *
from plot_functions_sensor import *



def run_wrapper(save_fig, station):
    """ Wrapper to full run on a station file """
  
    # getting the data
    data_clean_all, all_sensor_station, data_all_wmo = get_data(station, force_create=False)
 
    plot = Plot(station.split('_')[-1], save=save_fig)
 
    # making the plots
    series = plot.time_series( data_clean_all, label='')
    table = plot.sensor_table( all_sensor_station)
    # wmo_table = plot.wmo_bar_plot(data_clean_all)

    print("*** COMPLETED ***" , station )
    return series, table



""" Running """


stations = [s.split('_')[0] for s in os.listdir(merged) if 'before' not in s]

stations = [s for s in stations if '10393' in s or '06610' in s or '11035' in s or '82930' in s ]

# single or multiprocessing
POOL = True
save_fig = True



print(stations)

if __name__ == "__main__":

    if POOL:
        p = Pool(30)
        func = partial(run_wrapper, save_fig,)
        out = p.map(func, stations)   

    else:
        for stat in stations:
            s,t = run_wrapper(save_fig, stat)  ### to debug

            try:
                s,t = run_wrapper(save_fig, stat)
            except:
                print('Failed +++ ' , stat )
                pass

