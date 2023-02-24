import numpy
import numpy as np

from pyproj import Geod
g = Geod(ellps="WGS84")

def calc_height(t, p):
    '''
    t
    
    isotherm height formula
    z = -R*t0/g * ln(p/p0)
    z = -287.053*t0/9.80665 * ln(p/p0)
    
    polytrop height forumula
    z = t0/L * ((p/p0)**(-L*R/g) -1)
    L = −0.0065 K/m
    R = 287.053 J/(kg K)
    g = 9.80665 m/s2
    z = t0/−0.0065 * ((p/p0)**(0.0065*287.053/9.80665) -1)
    
    international height formula
    z = 288.15/0.0065 * (1- (p/1013.25)**(1/5.255))
    
    '''
    # from: https://www.cesm.ucar.edu/models/cesm1.1/cesm/cesmBbrowser/html_code/cam/tropopause.F90.html
    # https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2003GL018240
    
    SHR_CONST_AVOGAD  = 6.02214e26
    SHR_CONST_BOLTZ   = 1.38065e-23
    SHR_CONST_MWDAIR  = 28.966
    SHR_CONST_RGAS    = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ
    SHR_CONST_RDAIR   = SHR_CONST_RGAS/SHR_CONST_MWDAIR
    rair = SHR_CONST_RDAIR

    SHR_CONST_G       = 9.80616
    gravit = SHR_CONST_G

    SHR_CONST_CPDAIR  = 1.00464e3
    cappa        = (SHR_CONST_RGAS/SHR_CONST_MWDAIR)/SHR_CONST_CPDAIR
    cnst_kap = cappa

    cnst_faktor = -gravit/rair
    cnst_ka1    = cnst_kap - 1.
        
    z = []
    for i in range(len(t)):
        
        if i == 0:
            L = -0.0065
            height = t[i]/L * ((p[i]/101325)**(-L*287.053/9.80665) -1)
            z.append(height)
        else:
                   
            # dt/dz
            pmk= .5 * (p[i-1]**cnst_kap+p[i]**cnst_kap)
            pm = pmk**(1/cnst_kap)               
            a = (t[i-1]-t[i])/(p[i-1]**cnst_kap-p[i]**cnst_kap)
            b = t[i]-(a*p[i]**cnst_kap)
            tm = a * pmk + b               
            dtdp = a * cnst_kap * (pm**cnst_ka1)
            L = cnst_faktor*dtdp*pm/tm # dtdz
            if L == 0:
                L = -0.001

            height = t[i-1]/L * ((p[i]/p[i-1])**(-L*287.053/9.80665) -1)
            if np.isnan(height):
                z.append(z[-1])
            else:
                z.append(z[-1] + height)
    return z


def haversine(lon1, lat1, lon2, lat2):
    """
    lat1               starting latitude [°] [float]
    lon1               starting longitude [°] [float]
    lat2               end latitude [°] [float]
    lon2               end longitude [°] [float]
    
    returns:
        distance
        
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    
    # convert decimal degrees to radians 
    lon1 = numpy.radians(lon1)
    lat1 = numpy.radians(lat1)
    lon2 = numpy.radians(lon2)
    lat2 = numpy.radians(lat2)

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return c * r

def inverse_haversine(lat, lon, distance, direction):
    '''
    lat                actual latitude [°] [float]
    lon                actual longitude [°] [float]
    distance           distance to move [km] [float]
    direction          direction to move ['NORTH', 'EAST']
    
    returns:
        new_latitude, new_longitude
        
    inverse haversine calculation - point and distance to new point
    '''
    
    lat = numpy.radians(lat)
    lon = numpy.radians(lon)
    d = numpy.array(distance)
    r = 6371 #[km]
    if direction == "NORTH":
        brng = numpy.radians(0)
    elif direction == "EAST":
        brng = numpy.radians(90)
    else:
        return "error - not a valid direction"
    return_lat = numpy.arcsin(numpy.sin(lat) * numpy.cos(d / r) + numpy.cos(lat) * numpy.sin(d / r) * numpy.cos(brng))
    return_lon = lon + numpy.arctan2(numpy.sin(brng) * numpy.sin(d / r) * numpy.cos(lat), numpy.cos(d / r) - numpy.sin(lat) * numpy.sin(return_lat))

    return numpy.degrees(return_lat), numpy.degrees(return_lon)


def transport(lat, lon, u_dist, v_dist, transport_type):
    if transport_type == 'sphere':
        new_lat, new_lon = transport_sphere(lat, lon, u_dist, v_dist)
    else:
        new_lat, new_lon = transport_geod(lat, lon, u_dist, v_dist)
    return new_lat, new_lon


def transport_sphere(lat, lon, u_dist, v_dist):
    '''
    lat                actual latitude [°] [float]
    lon                actual longitude [°] [float]
    u_dist             longitudinal distance added to position [km] [float]
    v_dist             meridional distance added to position [km] [float]
    
    returns:
        new_latitude, new_longitude
    '''
    new_lat, new_lon = inverse_haversine(lat, lon, u_dist, "EAST")
    new_lat, new_lon = inverse_haversine(new_lat, new_lon, v_dist, "NORTH")
    return new_lat, new_lon


def transport_geod(lat, lon, u_dist, v_dist):
    '''
    lat                actual latitude [°] [float]
    lon                actual longitude [°] [float]
    u_dist             longitudinal distance added to position [km] [float]
    v_dist             meridional distance added to position [km] [float]
    
    returns:
        new_latitude, new_longitude
    '''
    new_lon, new_lat, backward_azim = g.fwd(lons=lon, lats=lat, az=90., dist=u_dist*1000.0)
    new_lon, new_lat, backward_azim = g.fwd(lons=new_lon, lats=new_lat, az=0., dist=v_dist*1000.0)
    return new_lat, new_lon


def trajectory(lat, lon, u, v, pressure, temperature, w_rs = 5.0, wind = 'mean', output='degree', transport_type='ellipsoid'):
    '''
    main function to calculate trajectories
    
    lat                station latitude [°] [int]
    lon                station longitude [°] [int]
    u                  eastward wind speed [m/s] [array - float]
    v                  northward wind speed [m/s] [array - float]
    pressure           pressure for given levels, ascending, [Pa] [array - float]
    temperature        temperature [K] [array - float]
    w_rs               radio sonde rising speed [m/s] [float]
    wind               wind calculation option ['mean', 'upper', 'lower']
    output             displacement output unit ['degree', 'km']
    transport_type     distance calculation ['sphere', ellipsoid'] 
    
    returns: 
        latitude displacement [output - array], longitude displacement [output - array],
        u wind shear [m/s - array], v wind shear [m/s - array],
        seconds since launch [s - array]
    '''
    
    # check if sorted correctly
    if pressure[0] < pressure[-1]:
        print("Please resort the input data - ascending order is necessary!")
        return 0,0,0,0
        
    z = calc_height(temperature, pressure) # m from K and Pa
    
    new_lat = lat
    new_lon = lon
    
    lat_displacement = [0.]
    lon_displacement = [0.]
    
    u_shear=[0.]
    v_shear=[0.]
    
    rts = [0]
    
    for i in range(len(z)):
        if i != 0:
            rising_time = (z[i]-z[i-1]) / w_rs
            rts.append(rts[-1] + rising_time)
            u_shear.append(u[i]-u[i-1])
            v_shear.append(v[i]-v[i-1])

            if wind == 'mean':
                if output == 'degree':
                    new_lat, new_lon = transport(new_lat, new_lon, (np.mean([u[i],u[i-1]]) * rising_time)/1000. , (np.mean([v[i],v[i-1]]) * rising_time)/1000., transport_type)
                elif output == 'km':
                    new_lon = (np.mean([u[i],u[i-1]]) * rising_time)/1000.
                    new_lat = (np.mean([v[i],v[i-1]]) * rising_time)/1000.
                                        
            elif wind == 'upper':
                new_lat, new_lon = transport(new_lat, new_lon, (u[i] * rising_time)/1000., (v[i] * rising_time)/1000., transport_type)
            elif wind == 'lower':
                new_lat, new_lon = transport(new_lat, new_lon, (u[i-1] * rising_time)/1000., (v[i-1] * rising_time)/1000., transport_type) 
            else:
                print('error: not a valid wind request')


            if output == 'degree':
                lat_displacement.append(new_lat - lat)
                lon_displacement.append(new_lon - lon)
            elif output == 'km':
                lat_displacement.append(new_lat)
                lon_displacement.append(new_lon)

    return lat_displacement, lon_displacement, np.array(u_shear), np.array(v_shear), rts


def trajectory_height(lat, lon, u, v, height, w_rs = 5.0, wind = 'mean', output='degree', transport_type='ellipsoid'):
    '''
    main function to calculate trajectories
    
    lat                station latitude [°] [int]
    lon                station longitude [°] [int]
    u                  eastward wind speed [m/s] [array - float]
    v                  northward wind speed [m/s] [array - float]
    height             height for given levels, ascending, [m] [array - float]
    temperature        temperature [K] [array - float]
    w_rs               radio sonde rising speed [m/s] [float]
    wind               wind calculation option ['mean', 'upper', 'lower']
    output             displacement output unit ['degree', 'km']
    transport_type     distance calculation ['sphere', ellipsoid'] 
    
    returns: 
        latitude displacement [output - array], longitude displacement [output - array],
        u wind shear [m/s - array], v wind shear [m/s - array],
        seconds since launch [s - array]
    '''
    
    # check if sorted correctly
    if height[0] > height[-1]:
        print("Please resort the input data - ascending order is necessary!")
        return 0,0,0,0
    
    z = height
    
    new_lat = lat
    new_lon = lon
    
    lat_displacement = [0.]
    lon_displacement = [0.]
    
    u_shear=[0.]
    v_shear=[0.]
    
    rts = [0]
    
    for i in range(len(z)):
        if i != 0:
            rising_time = (z[i]-z[i-1]) / w_rs
            rts.append(rts[-1] + rising_time)
            u_shear.append(u[i]-u[i-1])
            v_shear.append(v[i]-v[i-1])

            if wind == 'mean':
                if output == 'degree':
                    new_lat, new_lon = transport(new_lat, new_lon, (np.mean([u[i],u[i-1]]) * rising_time)/1000., (np.mean([v[i],v[i-1]]) * rising_time)/1000., transport_type)
                elif output == 'km':
                    new_lon = (np.mean([u[i],u[i-1]]) * rising_time)/1000.
                    new_lat = (np.mean([v[i],v[i-1]]) * rising_time)/1000.
                                        
            elif wind == 'upper':
                new_lat, new_lon = transport(new_lat, new_lon, (u[i] * rising_time)/1000., (v[i] * rising_time)/1000., transport_type)
            elif wind == 'lower':
                new_lat, new_lon = transport(new_lat, new_lon, (u[i-1] * rising_time)/1000., (v[i-1] * rising_time)/1000., transport_type) 
            else:
                print('error: not a valid wind request')


            if output == 'degree':
                lat_displacement.append(lat - new_lat)
                lon_displacement.append(lon - new_lon)
            elif output == 'km':
                lat_displacement.append(new_lat)
                lon_displacement.append(new_lon)

    return lat_displacement, lon_displacement, np.array(u_shear), np.array(v_shear), rts


