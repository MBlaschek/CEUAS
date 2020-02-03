import calendar
import re
from collections import defaultdict

import cdstoolbox as ct
import numpy as np

REANALYSIS_PERIOD = (1979, 2018)
REFERENCE_PERIOD = (1981, 2010)

REANALYSIS_PERIOD_STR = '-'.join(str(year) for year in REANALYSIS_PERIOD)
REFERENCE_PERIOD_STR = '-'.join(str(year) for year in REFERENCE_PERIOD)

DESCRIPTION = (
    f'## Click on a point or search for a city in the map below to discover a '
    f'range of local climate statistics for the period '
    f'{REANALYSIS_PERIOD_STR}.\n'
    f'### **Note:** The performance of this application is optimised for the '
    f'most populous cities in Europe. When clicking on the map or searching '
    f'for lower population and/or non-European cities, you might need to wait '
    f'a while for the data to be retrieved and processed.\n&nbsp;\n'
    f'This application is driven by '
    f'[ERA5](https://cds.climate.copernicus.eu/cdsapp#!/dataset/'
    f'reanalysis-era5-single-levels-monthly-means?tab=overview), the fifth '
    f'generation ECMWF atmospheric reanalysis of the global climate.'
)


class Koppen:
    """
    Class for applying the Köppen climate classification to a location based
    on its monthly average temperature and precipitation totals.

    """

    #: The five possible Koppen climate classifications
    CLIMATE_TYPES = ('dry', 'tropical', 'temperate', 'continental', 'polar')

    def __init__(self, tas_clim, precip_clim, latitude):
        self.tas_monthly = ct.basic.get_values(tas_clim)['values']

        self.precip_monthly = ct.basic.get_values(precip_clim)['values']

        self.lat = latitude

    @property
    def summer_precip(self):
        return self._get_summer(self.precip_monthly)

    @property
    def winter_precip(self):
        return self._get_winter(self.precip_monthly)

    @property
    def summer_tas(self):
        return self._get_summer(self.tas_monthly)

    @property
    def winter_tas(self):
        return self._get_winter(self.tas_monthly)

    def _get_summer(self, monthly_data):
        summer = [3, 4, 5, 6, 7, 8] if self.lat > 0 else [9, 10, 11, 0, 1, 2]
        return [monthly_data[m] for m in summer]

    def _get_winter(self, monthly_data):
        winter = [3, 4, 5, 6, 7, 8] if self.lat < 0 else [9, 10, 11, 0, 1, 2]
        return [monthly_data[m] for m in winter]

    @property
    def mean_tas(self):
        return sum(self.tas_monthly) / 12

    @property
    def tot_precip(self):
        return sum(self.precip_monthly)

    def classify(self):
        """
        Given monthly average temperatures and precipitation totals, classify
        the location into one of the five Köppen climate groups.

        """
        for climate_type in Koppen.CLIMATE_TYPES:
            climate = getattr(self, climate_type)()
            if climate:
                break
        else:
            climate = ''

        return climate

    def tropical(self):
        """
        Climate type which is hot year-round, and often very wet.

        """
        if min(self.tas_monthly) > 18:
            if min(self.precip_monthly) > 60:
                climate = 'wet equatorial'
            elif min(self.precip_monthly) > (100 - sum(self.precip_monthly) / 25):
                climate = 'tropical monsoon'
            else:
                climate = 'tropical savanna'
        else:
            climate = False

        return climate

    def dry(self):
        """
        Climate type which sees little precipitation.

        """
        threshold = 20 * self.mean_tas
        summer_precip_fraction = sum(self.summer_precip) / self.tot_precip
        if summer_precip_fraction >= 0.7:
            threshold += 280
        elif summer_precip_fraction >= 0.3:
            threshold += 140

        if self.tot_precip < threshold:
            if self.tot_precip < (20 * self.mean_tas + 280) / 2:
                climate = 'desert'
            else:
                climate = 'semi-arid'
            climate = f'{"hot" if self.mean_tas > 18 else "cold"} {climate}'
        else:
            climate = False

        return climate

    def temperate(self):
        """
        Climate type which sees mild temperatures and moderate precipitation.

        """
        if max(self.tas_monthly) >= 10 and -3 <= min(self.tas_monthly) <= 18:

            min_summer_precip = min(self.summer_precip)
            min_winter_precip = min(self.winter_precip)
            max_summer_precip = max(self.summer_precip)
            max_winter_precip = max(self.winter_precip)

            if max(self.tas_monthly) > 22:
                summer = 'hot'
            elif all([tas > 10 for tas in sorted(self.tas_monthly)[-4:]]):
                summer = 'warm'
            else:
                summer = 'cold'

            if all(min_summer_precip < t for t in (30, max_winter_precip / 3)):
                climate = f'{summer}-summer Mediterranean'
            elif min_winter_precip < max_summer_precip / 10:
                climate = {
                    'hot': 'humid subtropical',
                    'warm': 'subtropical highland',
                    'cold': 'cold subtropical highland',
                }[summer]
            else:
                climate = {
                    'hot': 'humid subtropical',
                    'warm': 'temperate oceanic',
                    'cold': 'subpolar oceanic',
                }[summer]

        else:
            climate = False

        return climate

    def continental(self):
        """
        Climate type which sees mild temperatures with cold winters.

        """
        if max(self.tas_monthly) >= 10 and min(self.tas_monthly) < -3:
            if min(self.tas_monthly) <= -38:
                climate = 'extremely cold subarctic climate'
            elif max(self.tas_monthly) > 22:
                climate = 'hot-summer continental'
            elif all([t > 10 for t in sorted(self.tas_monthly)[-4:]]):
                climate = 'warm-summer continental'
            else:
                climate = 'subarctic continental'
        else:
            climate = False

        return climate

    def polar(self):
        """
        Climate type which is very cold year-round.

        """
        if max(self.tas_monthly) < 10:
            if max(self.tas_monthly) > 0:
                climate = 'tundra'
            else:
                climate = 'ice cap'
        else:
            climate = False

        return climate


@ct.child(position='right')
@ct.output.markdown()  # Title
@ct.output.markdown()  # Koppen climate text
@ct.output.markdown()  # Temperature climatology title
@ct.output.markdown()  # Temperature climatology text
@ct.output.livefigure()  # Temperature climatology plot
@ct.output.markdown()  # Warming stripes title
@ct.output.livefigure()  # Warming stripes plot
@ct.output.markdown()  # Warming stripes text
@ct.output.markdown()  # Temperature indices title
@ct.output.livefigure()  # Temperature indices plot
@ct.output.markdown()  # Temperature indices text
@ct.output.markdown()  # Temperature indices (yearly) title
@ct.output.livefigure()  # Temperature indices (yearly) plot
@ct.output.markdown()  # Precipitation climatology title
@ct.output.livefigure()  # Precipitation climatology plot
@ct.output.markdown()  # Precipitation climatology text
@ct.output.markdown()  # Precipitation anomaly title
@ct.output.livefigure()  # Precipitation anomaly plot
@ct.output.markdown()  # Wind speed, direction & gust climatology title
@ct.output.livefigure()  # Wind speed, direction & gust climatology plot
@ct.output.markdown()  # Wind speed, direction & gust climatology text
def site_breakdown(**kwargs):
    """
    Generate a breakdown of climatological statistics for the given location.

    Application main steps:

    - Extract location-specific metadata (lat/lon, location name)
    - Extract a single grid point of data at the location specified
    - Generate temperature, precipitation and wind climate statistics and plots

    """
    location = kwargs.get('location')
    data = kwargs.get('data')
    print(location)
    lat = location['lat']
    lon = location['lon']

    try:
        # Extract the city name from the location dict
        loc = re.sub("[\(\[].*?[\)\]]", "", location['label']).rstrip(' ')
    except KeyError:
        # If there is no city name in the location dict, use lat and lon
        loc = (f'{abs(lat):.2f}°{"N" if lat > 0 else "S"}, '
               f'{abs(lon):.2f}°{"E" if lon > 0 else "W"}')

    site_data = extract_site_data(data, lat, lon)

    # Generate a title
    title = f'# Climate reanalysis for {loc}'

    # Calculate the Koppen climate classification of the given location
    koppen = Koppen(site_data['tas']['ref']['mean'],
                    site_data['precip']['ref']['mean'], lat).classify()
    if koppen:
        word = 'an' if koppen[0].lower() in 'aeiou' else 'a'
        koppen = (f'## {loc} had {word} **{koppen}** climate over the '
                  f'{REFERENCE_PERIOD_STR} period.')

    # Generate a breakdown of temperature statistics
    tas_output = tas_breakdown(site_data['tas'], loc)

    # Generate a breakdown of precipitation statistics
    precip_output = precip_breakdown(site_data['precip'], loc)

    # Generate a breakdown of wind statistics
    wind_output = wind_breakdown(site_data['wind'], loc)

    return (title, koppen, *tas_output, *precip_output, *wind_output)


def extract_site_data(data, lat, lon):
    """
    Extract data at the closest grid point to the given lat/lon values, from
    the pre-calculated climate statistics.

    """
    extract = lambda x: ct.geo.extract_point(x, lat=lat, lon=lon)

    site_data = dict()
    site_data['tas'] = dict()
    site_data['tas']['ref'] = defaultdict(dict)
    site_data['tas']['ref']['mean'] = extract(data['tas']['ref']['mean'])
    site_data['tas']['ref']['daily_min']['mean'] = extract(
        data['tas']['ref']['daily_min']['mean'])
    site_data['tas']['ref']['daily_min']['max'] = extract(
        data['tas']['ref']['daily_min']['max'])
    site_data['tas']['ref']['daily_min']['min'] = extract(
        data['tas']['ref']['daily_min']['min'])
    site_data['tas']['ref']['daily_max']['mean'] = extract(
        data['tas']['ref']['daily_max']['mean'])
    site_data['tas']['ref']['daily_max']['max'] = extract(
        data['tas']['ref']['daily_max']['max'])
    site_data['tas']['ref']['daily_max']['min'] = extract(
        data['tas']['ref']['daily_max']['min'])
    site_data['tas']['ref']['tropical_nights'] = extract(
        data['tas']['ref']['tropical_nights'])
    site_data['tas']['ref']['summer_days'] = extract(
        data['tas']['ref']['summer_days'])
    site_data['tas']['ref']['frost_days'] = extract(
        data['tas']['ref']['frost_days'])
    site_data['tas']['anomaly'] = extract(data['tas']['anomaly'])
    site_data['tas']['annual'] = dict()
    site_data['tas']['annual']['tropical_nights'] = extract(
        data['tas']['annual']['tropical_nights'])
    site_data['tas']['annual']['summer_days'] = extract(
        data['tas']['annual']['summer_days'])
    site_data['tas']['annual']['frost_days'] = extract(
        data['tas']['annual']['frost_days'])
    site_data['tas']['annual']['5y_tropical_nights'] = extract(
        data['tas']['annual']['5y_tropical_nights'])
    site_data['tas']['annual']['5y_summer_days'] = extract(
        data['tas']['annual']['5y_summer_days'])
    site_data['tas']['annual']['5y_frost_days'] = extract(
        data['tas']['annual']['5y_frost_days'])

    site_data['precip'] = dict()
    site_data['precip']['ref'] = dict()
    site_data['precip']['ref']['mean'] = extract(data['precip']['ref']['mean'])
    site_data['precip']['anomaly'] = extract(data['precip']['anomaly'])

    site_data['wind'] = dict()
    site_data['wind']['ref'] = defaultdict(dict)
    site_data['wind']['ref']['max_gust_mean'] = extract(
        data['wind']['ref']['max_gust_mean'])
    site_data['wind']['ref']['speed_mean'] = extract(
        data['wind']['ref']['speed_mean'])
    site_data['wind']['ref']['dir_mean'] = extract(
        data['wind']['ref']['dir_mean'])

    return site_data


def tas_breakdown(tas, loc):
    """
    Generate all of the temperature statistics and plots.

    """
    # Generate text describing the climatology
    climatology_title = (f'## **Monthly minimum, mean and maximum '
                         f'temperatures for {loc} ({REFERENCE_PERIOD_STR})**')
    climatology_text = tas_climatology_text(tas['ref']['mean'], loc)
    # Generate a temperature climatology plot
    climatology_plot = tas_climatology_plot(tas['ref'], loc)

    # Generate temperate anomaly stripes
    anomaly_title = (f'## **Warming stripes for {loc} '
                     f'({REANALYSIS_PERIOD_STR})**')
    anomaly_stripes = tas_anomaly_stripes(tas['anomaly'], loc)
    anomaly_text = tas_anomaly_text(loc)

    # Generate a plot of temperature indices
    indices_title = (f'## **Monthly temperature indices for {loc} '
                     f'({REFERENCE_PERIOD_STR})**')
    anomaly_stripes = tas_anomaly_stripes(tas['anomaly'], loc)
    indices_plot = tas_indices_plot(tas['ref'], loc)

    # Generate some text describing the temperature indices plot
    indices_text = tas_indices_text(loc)

    # Generate a plot of temperature index trends
    trends_title = (f'## **Annual temperature indices for {loc} '
                    f'({REANALYSIS_PERIOD_STR})**')
    trends_plot = tas_indices_trends_plot(tas['annual'], loc)

    return (climatology_title, climatology_text, climatology_plot,
            anomaly_title, anomaly_stripes, anomaly_text,
            indices_title, indices_plot, indices_text, trends_title,
            trends_plot)


def tas_climatology_text(mean, loc):
    """
    Generate a couple of sentences of text which describe the typical monthly
    and yearly temperatures at the given location.

    """
    monthly_means = ct.basic.get_values(mean)['values']

    mean_tas = np.mean(monthly_means)

    min_tas = min(monthly_means)
    min_month = calendar.month_name[monthly_means.index(min_tas) + 1]

    max_tas = max(monthly_means)
    max_month = calendar.month_name[monthly_means.index(max_tas) + 1]

    text = (f'## For the {REFERENCE_PERIOD_STR} reference period, the annual '
            f'average temperature in {loc} was **{mean_tas:.1f}**°C.\n'
            f'## Monthly average temperatures ranged from **{min_tas:.1f}**°C '
            f'({min_month}) to **{max_tas:.1f}**°C ({max_month}).')

    return text


def tas_climatology_plot(tas, loc):
    """
    Generate and plot some statistics which describe the typical monthly
    temperatures at the given location.

    """
    print(tas)
    # Set up the plot style
    layout_kwargs = {
        'xaxis': {
            'tickmode': 'array',
            'tickvals': list(range(1, 13)),
            'tickangle': 45,
            'ticktext': [calendar.month_name[i] for i in range(1, 13)],
            'title': {
                'text': 'Month of year',
            },
            'automargin': True,
        },
        'yaxis': {
            'title': {
                'text': 'Temperature (°C)',
            },
            'automargin': True,
        },
        'legend': {
            'xanchor': 'left',
            'yanchor': 'bottom',
            'y': 1.01,
            'x': 0,
        },
        'margin': {
            'r': 40,
        }
    }
    print(tas['daily_max']['min'])
    # Plot the statistics
    fig = ct.chart.line(tas['daily_max']['min'],
                        layout_kwargs=layout_kwargs,
                        scatter_kwargs=dict(
                            showlegend=False, hoverinfo='skip',
                            line={'shape': 'spline', 'color': '#f66',
                                  'width': 0},
                            mode='lines'))

    ct.chart.line(tas['daily_max']['max'], fig=fig,
                  name='Monthly range of daily maximum',
                  scatter_kwargs=dict(
                      fill='tonexty', fillcolor='rgba(255, 102, 102, 0.2)',
                      customdata=tas['daily_max']['min'],
                      hovertemplate=('%{customdata:.1f} - %{y:.1f}°C'
                                     '<extra></extra>'),
                      line={'shape': 'spline', 'color': '#ffe0e0', 'width': 0},
                      mode='lines'))

    ct.chart.line(tas['daily_min']['max'], fig=fig,
                  scatter_kwargs=dict(
                      showlegend=False, hoverinfo='skip',
                      line={'shape': 'spline', 'color': '#6cf', 'width': 0},
                      mode='lines'))

    ct.chart.line(tas['daily_min']['min'], fig=fig,
                  name='Monthly range of daily minimum',
                  scatter_kwargs=dict(
                      fill='tonexty', fillcolor='rgba(102, 204, 255, 0.2)',
                      customdata=tas['daily_min']['max'],
                      hovertemplate=('%{y:.1f} - %{customdata:.1f}°C'
                                     '<extra></extra>'),
                      line={'shape': 'spline', 'color': '#e0f5ff', 'width': 0},
                      mode='lines'))

    ct.chart.line(tas['daily_min']['mean'], fig=fig,
                  name='Monthly average of daily minimum',
                  scatter_kwargs=dict(
                      hovertemplate='%{y:.1f}°C<extra></extra>',
                      line={'shape': 'spline', 'color': '#6cf'},
                      mode='lines'))

    ct.chart.line(tas['mean'], fig=fig, name='Monthly average of daily average',
                  scatter_kwargs=dict(
                      hovertemplate='%{y:.1f}°C<extra></extra>',
                      line={'shape': 'spline', 'color': '#cc9804'},
                      mode='lines'))

    ct.chart.line(tas['daily_max']['mean'], fig=fig,
                  name='Monthly average of daily maximum',
                  scatter_kwargs=dict(
                      hovertemplate='%{y:.1f}°C<extra></extra>',
                      line={'shape': 'spline', 'color': '#f66'},
                      mode='lines'))

    return fig


def tas_anomaly_stripes(anomaly, loc):
    """
    Generate and plot yearly temperature anomalies as blue and red "warming"
    stripes.

    """
    customdata = [_inclusive_range(*REANALYSIS_PERIOD)]
    fig = ct.chart.warming_stripes(anomaly,
                                   height=100,
                                   customdata=customdata,
                                   hovertemplate=('%{customdata} (%{z:+.1f}°C)'
                                                  '<extra></extra>'))

    return fig


def tas_anomaly_text(loc):
    """
    Generate a couple of sentences of text which describe (and credit) warming
    stripes.

    """
    text = (f'Warming stripes inspired by '
            f'[#ShowYourStripes](https://showyourstripes.info/).\n&nbsp;\n'
            f'## Warming stripes provide an at-a-glance view of yearly '
            f'temperature trends in {loc} for the period '
            f'{REANALYSIS_PERIOD_STR}.\n'
            f'## The colour of each stripe represents the **temperature '
            f'anomaly** for a given year, or how much warmer (red) or '
            f'colder (blue) that year was relative to the **long-term '
            f'reference** period of {REFERENCE_PERIOD_STR}.\n&nbsp;\n')

    return text


def tas_indices_plot(tas, loc):
    """
    Generate and plot statistics which describe the percentage of days in each
    monthly which typically fall into temperature-related categories.

    """
    # Set up the plot style
    layout_kwargs = {
        'xaxis': {
            'tickmode': 'array',
            'tickvals': list(range(1, 13)),
            'tickangle': 45,
            'ticktext': [calendar.month_name[i] for i in range(1, 13)],
            'title': {
                'text': 'Month of year',
            },
            'automargin': True,
        },
        'yaxis': {
            'tickmode': 'array',
            'tickformat': '%{n:.0%}',
            'title': {
                'text': 'Percentage of days',
            },
        },
        'legend': {
            'xanchor': 'left',
            'yanchor': 'bottom',
            'y': 1.01,
            'x': 0,
        },
        'margin': {
            'r': 40,
        }
    }
    print(tas['frost_days'])
    print(tas['tropical_nights'])
    print(tas['summer_days'])
    # Generate the plot
    fd = tas['frost_days']
    fig = ct.chart.line(fd, layout_kwargs=layout_kwargs,
                        name=('Frost days'),
                        scatter_kwargs=dict(
                            hovertemplate='%{y:.1%}<extra>Frost days</extra>',
                            line={'shape': 'spline', 'color': '#6cf'},
                            mode='lines'))

    ct.chart.line(tas['tropical_nights'], fig=fig,
                  name='Tropical nights',
                  scatter_kwargs=dict(
                      hovertemplate='%{y:.1%}<extra>Tropical nights</extra>',
                      line={'shape': 'spline', 'color': '#d66'},
                      mode='lines'))

    ct.chart.line(tas['summer_days'], fig=fig,
                  name='Summer days',
                  scatter_kwargs=dict(
                      hovertemplate='%{y:.1%}<extra>Summer days</extra>',
                      line={'shape': 'spline', 'color': '#fb3'},
                      mode='lines'))

    return fig


def tas_indices_text(loc):
    """
    Generate a few lines of text which describe the aboove temperature indices
    plot.

    """
    text = (f'## The above graph shows the typical monthly percentage of '
            f'days in {loc} which are classified as:\n'
            f'## **Frost days**: a day in which the minimum temperature is '
            f'below 0°C;\n'
            f'## **Summer days**: a day in which the maximum temperature is '
            f'above 25°C;\n'
            f'## **Tropical nights**: a day in which the minimum temperature '
            f'is above 20°C.\n'
            f'## The graph below shows the **yearly** percentage of each of '
            f'the day types listed above, along with a **five year rolling '
            f'average**, for the {REANALYSIS_PERIOD_STR} period.')

    return text


def tas_indices_trends_plot(tas, loc):
    """
    Generate and plot statistics which describe the percentage of days in each
    year (including a five-year running average) which typically fall into
    temperature-related categories.

    """
    # Set up the plot style
    layout_kwargs = {
        'xaxis': {
            'title': {
                'text': 'Year',
            },
            'tickmode': 'array',
            'tickvals': _inclusive_range(*REANALYSIS_PERIOD),
            'ticktext': _inclusive_range(*REANALYSIS_PERIOD),
            'hoverformat': '%Y',
            'automargin': True,
        },
        'yaxis': {
            'tickmode': 'array',
            'tickformat': '%{n:.0%}',
            'title': {
                'text': 'Percentage of days',
            },
        },
        'legend': {
            'xanchor': 'left',
            'yanchor': 'bottom',
            'y': 1.01,
            'x': 0,
        },
        'margin': {
            'r': 30,
        }
    }

    five_year_means = list(zip(
        _inclusive_range(REANALYSIS_PERIOD[0], REANALYSIS_PERIOD[-1] - 5),
        _inclusive_range(REANALYSIS_PERIOD[0] + 5, REANALYSIS_PERIOD[-1])
    ))
    five_year_means = ['-'.join(str(i) for i in pd) for pd in five_year_means]

    # Generate the plot
    fig = ct.chart.line(tas['frost_days'], layout_kwargs=layout_kwargs,
                        name=('Percentage of days with minimum temperature < '
                              '0°C'),
                        scatter_kwargs=dict(
                            showlegend=False,
                            hoverinfo='skip',
                            line={'shape': 'spline',
                                  'color': 'rgba(102, 204, 255, 0.3)'},
                            mode='lines'))

    ct.chart.line(tas['tropical_nights'], fig=fig,
                  name='Percentage of days with minimum temperature > 20°C',
                  scatter_kwargs=dict(
                      showlegend=False,
                      hoverinfo='skip',
                      line={'shape': 'spline',
                            'color': 'rgba(221, 102, 102, 0.3)'},
                      mode='lines'))

    ct.chart.line(tas['summer_days'], fig=fig,
                  name='Percentage of days with maximum temperature > 25°C',
                  scatter_kwargs=dict(
                      showlegend=False,
                      hoverinfo='skip',
                      line={'shape': 'spline',
                            'color': 'rgba(255, 187, 51, 0.3)'},
                      mode='lines'))

    ct.chart.line(tas['5y_summer_days'], fig=fig,
                  name='Summer days (5 year average)',
                  scatter_kwargs=dict(
                      customdata=five_year_means,
                      hovertemplate='%{y:.1%}<extra>%{customdata}</extra>',
                      line={'shape': 'spline', 'color': '#fb3'},
                      mode='lines'))

    ct.chart.line(tas['5y_tropical_nights'], fig=fig,
                  name='Tropical nights (5 year average)',
                  scatter_kwargs=dict(
                      customdata=five_year_means,
                      hovertemplate='%{y:.1%}<extra>%{customdata}</extra>',
                      line={'shape': 'spline', 'color': '#d66'},
                      mode='lines'))

    ct.chart.line(tas['5y_frost_days'], fig=fig,
                  name='Frost days (5 year average)',
                  scatter_kwargs=dict(
                      customdata=five_year_means,
                      hovertemplate='%{y:.1%}<extra>%{customdata}</extra>',
                      line={'shape': 'spline', 'color': '#6cf'},
                      mode='lines'))

    return fig


def precip_breakdown(precip, loc):
    """
    Generate all of the precipitation statistics and plots.

    """
    # Generate a precipitation climatology plot
    climatology_title = (f'## **Monthly average total precipitation for {loc} '
                         f'({REFERENCE_PERIOD_STR})**')
    climatology_plot = precip_climatology_plot(precip['ref']['mean'], loc)
    # Generate text describing the climatology
    climatology_text = precip_climatology_text(precip['ref']['mean'], loc)

    # Generate a precipitation anomaly plot
    anomaly_title = (f'## **Annual precipitation anomaly for {loc} '
                     f'({REANALYSIS_PERIOD_STR})**')
    anomaly_plot = precip_anomaly_plot(precip['anomaly'],
                                       precip['ref']['mean'], loc)

    return (climatology_title, climatology_plot, climatology_text,
            anomaly_title, anomaly_plot)


def precip_climatology_plot(mean, loc):
    """
    Generate and plot statistics which describe the typical monthly
    precipitation totals at the given location.

    """
    # Set up the style of the plot
    layout_kwargs = {
        'xaxis': {
            'tickmode': 'array',
            'tickvals': list(range(1, 13)),
            'ticktext': [calendar.month_name[i] for i in range(1, 13)],
            'title': {
                'text': 'Month of year',
            },
            'automargin': True,
        },
        'yaxis': {
            'title': {
                'text': 'Precipitation (mm)',
            },
        },
    }

    fig = ct.chart.bar(mean, layout_kwargs=layout_kwargs,
                       bar_kwargs={'name': '',
                                   'hovertemplate': '%{y:.1f}mm'})

    return fig


def precip_climatology_text(mean, loc):
    """
    Generate a couple of sentences of text which describe the typical monthly
    and yearly precipitation totals at the given location.

    """
    # Convert the precip date from m/s rate to monthly totals
    mean = ct.basic.get_values(mean)['values']

    total_precip = sum(mean)

    min_precip = min(mean)
    min_month = calendar.month_name[mean.index(min_precip) + 1]

    max_precip = max(mean)
    max_month = calendar.month_name[mean.index(max_precip) + 1]

    text = (f'## For the {REFERENCE_PERIOD_STR} reference period, the average '
            f'annual total precipitation in {loc} was **{total_precip:.0f}** '
            f'mm.\n'
            f'## Monthly average precipitation totals ranged from '
            f'**{min_precip:.0f}** mm ({min_month}) to **{max_precip:.0f}** '
            f'mm ({max_month}).\n'
            f'## The plot below shows the **precipitation anomaly** for each '
            f'year in the {REANALYSIS_PERIOD_STR} period, or how much more '
            f'(blue) or less (red) precipitation fell each year as a '
            f'percentage relative to the long-term reference period of '
            f'{REFERENCE_PERIOD_STR}.\n')

    return text


def precip_anomaly_plot(anomaly, mean, loc):
    """
    Generate and plot yearly precipitation anomalies.

    """
    mean = ct.basic.get_values(mean)['values']

    total_precip = sum(mean)

    anomaly = anomaly / total_precip

    # Set up the plot style
    layout_kwargs = {
        'xaxis': {
            'tickmode': 'array',
            'tickvals': _inclusive_range(*REANALYSIS_PERIOD),
            'ticktext': _inclusive_range(*REANALYSIS_PERIOD),
            'hoverformat': '%Y',
            'title': {
                'text': 'Year',
            },
            'automargin': True,
        },
        'yaxis': {
            'tickmode': 'array',
            'tickformat': '%{n:+.0%}',
            'title': {
                'text': 'Anomaly (%)',
            },
        },
        'margin': {
            'r': 10
        }
    }

    bar_color = ct.cube.where(anomaly > 0, 'rgb(0.19,0.49,0.72)',
                              'rgb(0.77,0.23,0.23)')

    # Generate the plot
    fig = ct.chart.bar(anomaly, layout_kwargs=layout_kwargs,
                       bar_kwargs={'name': '',
                                   'marker': {'color': bar_color},
                                   'hovertemplate': '%{y:+.0%}'})

    return fig


def wind_breakdown(wind, loc):
    """
    Generate all of the wind speed/direction/gust statistics and plots.

    """
    # Generate a wind speed, direction and gust climatology plot
    climatology_title = (f'## **Monthly average maximum wind gust, wind speed '
                         f'and wind direction for {loc} '
                         f'({REFERENCE_PERIOD_STR})**')
    climatology_plot = wind_climatology_plot(wind['ref']['max_gust_mean'],
                                             wind['ref']['speed_mean'],
                                             wind['ref']['dir_mean'], loc)

    climatology_text = wind_climatology_text(wind['ref']['max_gust_mean'],
                                             wind['ref']['speed_mean'],
                                             wind['ref']['dir_mean'], loc)

    return climatology_title, climatology_plot, climatology_text


def wind_climatology_plot(max_gust_clim, speed_clim, dir_clim, loc):
    """
    Generate and plot statistics which describe the typical monthly
    wind speed, direction and maximum gust speed at the given location.

    """
    # Generate wind direction arrows as chart annotations
    arrow_dict = {
        'xref': 'x',
        'yref': 'y',
        'axref': 'x',
        'ayref': 'y',
        'text': '?',
        'align': 'center',
        'valign': 'middle',
        'font': {'size': 30},
    }

    dirs = ct.basic.get_values(dir_clim)['values']

    # The y coordinate of each arrow should be 25% higher than the max gust
    y = ct.basic.get_values(ct.math.max(max_gust_clim))['values'][0] * 5 / 4

    arrows = []
    for x in range(1, 13):
        angle_dict = {'x': x, 'y': y, 'ax': x, 'ay': y,
                      'textangle': dirs[x - 1] + 90}
        arrows.append({**arrow_dict, **angle_dict})

    # Set up the style of the plot
    layout_kwargs = {
        'xaxis': {
            'tickmode': 'array',
            'tickvals': list(range(1, 13)),
            'ticktext': [calendar.month_name[i] for i in range(1, 13)],
            'title': {
                'text': 'Month of year',
            },
            'automargin': True,
        },
        'yaxis': {
            'title': {
                'text': 'Wind speed (ms<sup>-1</sup>)',
            },
            'rangemode': 'tozero',
        },
        'legend': {
            'xanchor': 'left',
            'yanchor': 'bottom',
            'y': 1.01,
            'x': 0,
        },
        'annotations': arrows,
    }

    # Generate the plot
    hovertemplate = '%{y:.1f}ms<sup>-1</sup><extra></extra>'

    fig = ct.chart.line(max_gust_clim, layout_kwargs=layout_kwargs,
                        scatter_kwargs={'name': 'Average maximum wind gust speed',
                                        'hovertemplate': hovertemplate,
                                        'line': {'shape': 'spline',
                                                 'color': '#c95eed'},
                                        'mode': 'lines'})

    ct.chart.line(speed_clim, fig=fig,
                  scatter_kwargs={'name': 'Average wind speed',
                                  'hovertemplate': hovertemplate,
                                  'line': {'shape': 'spline',
                                           'color': '#845095'},
                                  'mode': 'lines'})

    return fig


def wind_climatology_text(max_gust_clim, speed_clim, dir_clim, loc):
    """
    Generate a couple of sentences of text which describe the typical monthly
    and yearly wind speed/diredtion/gust at the given location.

    """
    gusts = ct.basic.get_values(max_gust_clim)['values']
    speeds = ct.basic.get_values(speed_clim)['values']
    dirs = ct.basic.get_values(dir_clim)['values']

    rounded_dirs = [45 * round(x / 45) % 360 for x in dirs]
    rounded_dirs, counts = np.unique(rounded_dirs, return_counts=True)
    prevailing = rounded_dirs[np.argmax(counts)]
    prevailing = {
        0.0: 'northerly',
        45.0: 'northeasterly',
        90.0: 'easterly',
        135.0: 'southeasterly',
        180.0: 'southerly',
        225.0: 'southwesterly',
        270.0: 'westerly',
        315.0: 'northwesterly',
    }[prevailing]

    min_speed = min(speeds)
    min_speed_month = calendar.month_name[speeds.index(min_speed) + 1]

    max_speed = max(speeds)
    max_speed_month = calendar.month_name[speeds.index(max_speed) + 1]

    min_gust = min(gusts)
    min_gust_month = calendar.month_name[gusts.index(min_gust) + 1]

    max_gust = max(gusts)
    max_gust_month = calendar.month_name[gusts.index(max_gust) + 1]

    text = (f'## For the {REFERENCE_PERIOD_STR} period, the modal wind '
            f'direction in {loc} was **{prevailing}**.\n'
            f'## Monthly average wind speeds ranged from **{min_speed:.1f}** '
            f'm/s ({min_speed_month}) to **{max_speed:.1f}** m/s '
            f'({max_speed_month}), and monthly average maximum wind gust '
            f'speeds ranged from **{min_gust:.1f}** m/s ({min_gust_month}) to '
            f'**{max_gust:.1f}** m/s ({max_gust_month}).')

    return text


layout = {
    'output_ncols': 1,
    'output_align': 'left',
    'output_width': 0.9,
}


@ct.application(layout=layout)
@ct.input.constant('hidden_divider', '', label=' ')  # Hacky layout fix
@ct.output.markdown()
@ct.output.livemap(geocode_city=site_breakdown, click_on_map=site_breakdown,
                   height=70)
def climate_explorer(hidden_divider='', cache_child=None):
    """
    An application for exploring historical climate statistics for many cities
    around the world through an interactive map interface.

    The application consists of:

    - Some introductory text
    - A livemap (interactive map) which plots global average temperatures,
      precipitation and wind speeds for the 1981-2010 preiod
    - A panel describing climate information for locations selected in the map,
      generated by the `site_breakdown` child application

    """
    # Download the input ERA5 data from the CDS catalogue
    data = [ct.catalogue.retrieve(
        'reanalysis-era5-single-levels',
        {
            'variable': '2m_temperature',
            'grid': ['3', '3'],
            'product_type': 'reanalysis',
            'year': _inclusive_range(*REANALYSIS_PERIOD),
            'month': [f'{month:02d}' for month in _inclusive_range(1, 12)],
            'day': [f'{day:02d}' for day in _inclusive_range(1, 31)],
            'time': [f'{hour:02d}:00' for hour in range(0, 24, 3)],
        }
    )]
    data += ct.catalogue.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'variable': ['2m_temperature',
                         '10m_u_component_of_wind', '10m_v_component_of_wind',
                         'instantaneous_10m_wind_gust',
                         'mean_total_precipitation_rate'],
            'grid': ['3', '3'],
            'product_type': 'monthly_averaged_reanalysis',
            'year': _inclusive_range(*REANALYSIS_PERIOD),
            'month': [f'{month:02d}' for month in _inclusive_range(1, 12)],
            'time': '00:00',
        }
    )

    # Generate global averages to plot on the livemap
    map_data = generate_map_data(data)

    # Generate stats to cache and use in the child app
    app_data = generate_app_data(data)

    fig = ct.livemap.plot(map_data, crs='EPSG3857',
                          click_kwargs={'data': app_data})

    if cache_child is not None:
        site_breakdown(location=cache_child, data=app_data)

    return DESCRIPTION, fig


def _inclusive_range(*args):
    """
    Utility function for handling inclusive date ranges.

    """
    args = list(args)
    if len(args) == 0:
        args[0] += 1
    elif len(args) >= 2:
        args[1] += 1
    return list(range(*args))


def generate_map_data(data):
    """
    Calculate the mean data fields to plot on the interactive map.

    """
    map_data = []

    # Temperature mean
    map_data.append(ct.cube.average(data[1], dim='time'))

    # Wind speed mean
    u, v = data[2], data[3]
    speed = ct.math.sqrt(ct.operator.add(ct.math.square(u), ct.math.square(v)))

    map_data.append(ct.cube.average(speed, dim='time'))

    # Precipitation mean
    map_data.append(ct.cube.average(data[5], dim='time') * 3.154e7)

    # Modify the metadata so that the right colour schemes are used
    map_data[-1] = ct.cdm.rename(map_data[-1], 'precipitation_amount')
    map_data[-1] = ct.cdm.update_attributes(
        map_data[-1], {'standard_name': 'foobar'})
    map_data[-1] = ct.cdm.update_attributes(
        map_data[-1], {'cds_hidden_name': 'era5_explorer_yearly_precipitation'}
    )
    map_data[-1] = ct.cdm.update_attributes(
        map_data[-1], {'long_name': 'Average yearly precipitation'})

    map_data[1] = ct.cdm.rename(map_data[1], '10m_wind_speed')
    map_data[1] = ct.cdm.update_attributes(
        map_data[1], {'standard_name': '10m_wind_speed'})
    map_data[1] = ct.cdm.update_attributes(
        map_data[1], {'long_name': 'Average wind speed'})

    map_data[0] = ct.cdm.update_attributes(
        map_data[0], {'cds_hidden_name': 'era5_explorer_air_temperature'})
    map_data[0] = ct.cdm.update_attributes(
        map_data[0], {'standard_name': 'foobar'})
    map_data[0] = ct.cdm.update_attributes(
        map_data[0], {'long_name': 'Average temperature'})
    map_data[0] = ct.units.convert_units(map_data[0], 'celsius')

    return map_data


def generate_app_data(data):
    """
    Calculate all of the climate statistics to be extracted and plotted in the
    site_breakdown.

    """
    hrly_tas, tas, u, v, gust, precip = data

    app_data = defaultdict(dict)

    # Calculate temperature statistics
    hrly_tas = ct.units.convert_units(hrly_tas, 'celsius', source_units='K')
    tas = ct.units.convert_units(tas, 'celsius', source_units='K')
    reference_hrly_tas = ct.cube.select(hrly_tas,
                                        start_time=f'{REFERENCE_PERIOD[0]}-01-01',
                                        stop_time=f'{REFERENCE_PERIOD[1]}-12-31')
    reference_tas = ct.cube.select(tas,
                                   start_time=f'{REFERENCE_PERIOD[0]}-01-01',
                                   stop_time=f'{REFERENCE_PERIOD[1]}-12-31')

    tas_data = defaultdict(dict)
    tas_data['ref'] = defaultdict(dict)

    daily_min = ct.climate.daily_min(reference_hrly_tas)
    daily_max = ct.climate.daily_max(reference_hrly_tas)

    tas_data['ref']['mean'] = ct.climate.climatology_mean(reference_tas)

    tas_data['ref']['daily_min']['mean'] = (
        ct.climate.climatology_mean(daily_min))
    tas_data['ref']['daily_min']['max'] = (
        ct.climate.climatology_perc(daily_min, percentiles=[100])[0])
    tas_data['ref']['daily_min']['min'] = (
        ct.climate.climatology_perc(daily_min, percentiles=[0])[0])

    tas_data['ref']['daily_max']['mean'] = (
        ct.climate.climatology_mean(daily_max))
    tas_data['ref']['daily_max']['max'] = (
        ct.climate.climatology_perc(daily_max, percentiles=[100])[0])
    tas_data['ref']['daily_max']['min'] = (
        ct.climate.climatology_perc(daily_max, percentiles=[0])[0])

    tas_data['ref']['tropical_nights'] = (
        ct.climate.climatology_mean(daily_min > 20))
    tas_data['ref']['summer_days'] = (
        ct.climate.climatology_mean(daily_max > 25))
    tas_data['ref']['frost_days'] = ct.climate.climatology_mean(daily_min < 0)

    tas_data['anomaly'] = ct.climate.anomaly(
        ct.cube.resample(tas, freq='year'),
        interval=[str(year) for year in REFERENCE_PERIOD])

    daily_min = ct.climate.daily_min(hrly_tas)
    daily_max = ct.climate.daily_max(hrly_tas)

    tas_data['annual']['tropical_nights'] = ct.cube.resample(
        daily_min > 20, freq='year')
    tas_data['annual']['summer_days'] = ct.cube.resample(
        daily_max > 25, freq='year')
    tas_data['annual']['frost_days'] = ct.cube.resample(
        daily_min < 0, freq='year')

    tas_data['annual']['5y_tropical_nights'] = ct.cube.moving_window(
        ct.cube.resample(daily_min > 20, freq='year'), time=5, center=True)
    tas_data['annual']['5y_summer_days'] = ct.cube.moving_window(
        ct.cube.resample(daily_max > 25, freq='year'), time=5, center=True)
    tas_data['annual']['5y_frost_days'] = ct.cube.moving_window(
        ct.cube.resample(daily_min < 0, freq='year'), time=5, center=True)

    app_data['tas'] = tas_data

    # Calculate precipitation statistics
    ref_precip = ct.cube.select(precip,
                                start_time=f'{REFERENCE_PERIOD[0]}-01-01',
                                stop_time=f'{REFERENCE_PERIOD[1]}-12-31')

    precip_data = defaultdict(dict)

    precip_data['ref']['mean'] = ct.climate.climatology_mean(
        ref_precip) * (3.154e7) / 12

    precip_data['anomaly'] = ct.climate.anomaly(
        ct.cube.resample(precip * 3.154e7, freq='year'),
        interval=[str(year) for year in REFERENCE_PERIOD])

    app_data['precip'] = precip_data

    # Calculate wind statistics
    reference_u = ct.cube.select(u,
                                 start_time=f'{REFERENCE_PERIOD[0]}-01-01',
                                 stop_time=f'{REFERENCE_PERIOD[1]}-12-31')
    reference_v = ct.cube.select(v,
                                 start_time=f'{REFERENCE_PERIOD[0]}-01-01',
                                 stop_time=f'{REFERENCE_PERIOD[1]}-12-31')
    reference_gust = ct.cube.select(gust,
                                    start_time=f'{REFERENCE_PERIOD[0]}-01-01',
                                    stop_time=f'{REFERENCE_PERIOD[1]}-12-31')

    wind_data = defaultdict(dict)

    wind_data['ref']['max_gust_mean'] = ct.climate.climatology_mean(
        ct.climate.daily_max(reference_gust))

    wind_data['ref']['speed_mean'] = ct.climate.climatology_mean(
        ct.math.sqrt(ct.operator.add(ct.math.square(reference_u),
                                     ct.math.square(reference_v))))

    u_clim = ct.climate.climatology_mean(reference_u)
    v_clim = ct.climate.climatology_mean(reference_v)
    wind_data['ref']['dir_mean'] = (
        ct.math.rad2deg(ct.math.arctan2(-u_clim, -v_clim)))

    app_data['wind'] = wind_data

    return app_data
