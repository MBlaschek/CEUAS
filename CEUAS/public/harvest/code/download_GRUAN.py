import cdsapi

import cdsapi

c = cdsapi.Client()

c.retrieve(
    'insitu-observations-gruan-reference-network',
    {
        'format': 'csv-lev.zip',
        'year': '2007',
        'month': '01',
        'day': [
            '01', '03', '04',
            '06', '07', '13',
            '15', '16', '17',
            '18', '19', '20',
            '21', '22', '23',
            '24', '25', '26',
            '27', '28', '29',
        ],
        'variable': [
            'air_temperature', 'air_temperature_post_processing_radiation_correction', 'air_temperature_random_uncertainty',
            'air_temperature_systematic_uncertainty', 'air_temperature_total_uncertainty', 'altitude',
            'altitude_total_uncertainty', 'eastward_wind_component', 'frost_point_temperature',
            'geopotential_height', 'northward_wind_component', 'relative_humidity',
            'relative_humidity_effective_vertical_resolution', 'relative_humidity_post_processing_radiation_correction', 'relative_humidity_random_uncertainty',
            'relative_humidity_systematic_uncertainty', 'relative_humidity_total_uncertainty', 'shortwave_radiation',
            'shortwave_radiation_total_uncertainty', 'vertical_speed_of_radiosonde', 'water_vapor_volume_mixing_ratio',
            'wind_from_direction', 'wind_from_direction_total_uncertainty', 'wind_speed',
            'wind_speed_total_uncertainty',
        ],
    },
    'download.csv-lev.zip')
