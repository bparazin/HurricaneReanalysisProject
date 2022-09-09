from HurdatReader import *
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import cartopy.crs as ccrs
from geopy import distance
import scipy.stats as st
import datetime as dt
SECONDS_PER_HOUR = 3600


def positionBeforeArrival(storm, tMinusTime, lon_min, lon_max, lat_min, lat_max, category_statistics=False):

    if type(tMinusTime) == int or type(tMinusTime) == float:
        tMinusTime = dt.timedelta(hours=-tMinusTime)

    keyTime = None
    for i, (lon, lat) in enumerate(zip(storm['data']['lon'], storm['data']['lat'])):
        if lon_min <= lon <= lon_max and lat_min <= lat <= lat_max:
            keyTime = storm['data']['datetime'][i] + tMinusTime
            break

    if keyTime is None:  # this means it doesn't enter our bounding box and can't be plotted
        return None, None

    bounds = getSynopticBounds(storm, keyTime)

    if bounds is None:  # this means we are too early/late in data and it can't be plotted
        return None, None

    if type(bounds) == int:

        if bounds != len(storm['data']['lat']) - 1:

            velocity = (distance.distance((storm['data']['lat'][bounds], storm['data']['lon'][bounds]),
                                          (storm['data']['lat'][bounds + 1], storm['data']['lon'][bounds + 1])).km /
                        (storm['data']['datetime'][bounds + 1] - storm['data']['datetime'][
                            bounds]).total_seconds()) * SECONDS_PER_HOUR

            bearing = twoPointBearing(storm['data']['lon'][bounds], storm['data']['lat'][bounds],
                                      storm['data']['lon'][bounds + 1], storm['data']['lat'][bounds + 1])
        else:
            velocity = (distance.distance((storm['data']['lat'][bounds - 1], storm['data']['lon'][bounds - 1]),
                                          (storm['data']['lat'][bounds], storm['data']['lon'][bounds])).km /
                        (storm['data']['datetime'][bounds] - storm['data']['datetime'][
                            bounds - 1]).total_seconds()) * SECONDS_PER_HOUR

            bearing = twoPointBearing(storm['data']['lon'][bounds - 1], storm['data']['lat'][bounds - 1],
                                      storm['data']['lon'][bounds], storm['data']['lat'][bounds])

        if not category_statistics:
            return ((storm['data']['lon'][bounds], storm['data']['lat'][bounds]), storm), (velocity, bearing)
        else:
            location_row = storm['data'].iloc[bounds].copy()
            location_row['system status'] = statusNearLocation(storm, lon_min, lon_max, lat_min, lat_max, returnTime=False)
            
            return location_row
    else:
        # god i hate spherical coordinates
        velocity = (distance.distance((storm['data']['lat'][bounds[0]], storm['data']['lon'][bounds[0]]),
                                      (storm['data']['lat'][bounds[1]], storm['data']['lon'][bounds[1]])).km /
                    (storm['data']['datetime'][bounds[1]] - storm['data']['datetime'][
                        bounds[0]]).total_seconds()) * SECONDS_PER_HOUR
        time = (keyTime - storm['data']['datetime'][bounds[0]]).total_seconds() / SECONDS_PER_HOUR
        stormDistance = velocity * time
        bearing = twoPointBearing(storm['data']['lon'][bounds[0]], storm['data']['lat'][bounds[0]],
                                  storm['data']['lon'][bounds[1]], storm['data']['lat'][bounds[1]])
        dest = distance.distance(kilometers=stormDistance).destination(
            (storm['data']['lat'][bounds[0]], storm['data']['lon'][bounds[0]]),
            bearing=bearing)
        if not category_statistics:
            return ((dest.longitude, dest.latitude), storm), (velocity, bearing)
        else:
            previous_entry = storm['data'].iloc[bounds[0]]
            next_entry = storm['data'].iloc[bounds[1]]
            intermediate_interp = intermediate_series(previous_entry, next_entry, keyTime - storm['data']['datetime'][bounds[0]])
            intermediate_interp['system status'] = statusNearLocation(storm, lon_min, lon_max, lat_min, lat_max, returnTime=False)
            return intermediate_interp
        
def intermediate_series(previous_entry, next_entry, time_gap):
    intermediate_step = pd.Series(dtype = object)
    intermediate_step['datetime'] = previous_entry['datetime'] + time_gap
    intermediate_step['record identifier'] = None
    velocity = (distance.distance((previous_entry['lat'], previous_entry['lon']),
                                  (next_entry['lat'], next_entry['lon'])).km /
                (next_entry['datetime'] - previous_entry['datetime']).total_seconds()) * SECONDS_PER_HOUR
    time = (time_gap).total_seconds() / SECONDS_PER_HOUR
    stormDistance = velocity * time
    bearing = twoPointBearing(previous_entry['lon'], previous_entry['lat'],
                              next_entry['lon'], next_entry['lat'])
    dest = distance.distance(kilometers=stormDistance).destination(
        (previous_entry['lat'], previous_entry['lon']),
        bearing=bearing)
    
    time_between = (next_entry['datetime'] - previous_entry['datetime']).total_seconds() * SECONDS_PER_HOUR
    
    intermediate_step['lat'] = dest.latitude
    intermediate_step['lon'] = dest.longitude
    
    intermediate_step['Max sustained wind (knots)'] = intermediate_observation(previous_entry['Max sustained wind (knots)'], 
                                                                               next_entry['Max sustained wind (knots)'],
                                                                               time_between,
                                                                               time_gap.total_seconds() * SECONDS_PER_HOUR)
    intermediate_step['34 knot wind radii max extent (NM)'] = intermediate_observation(previous_entry['34 knot wind radii max extent (NM)'], 
                                                                                       next_entry['34 knot wind radii max extent (NM)'],
                                                                                       time_between,
                                                                                       time_gap.total_seconds() * SECONDS_PER_HOUR)
    intermediate_step['50 knot wind radii max extent (NM)'] = intermediate_observation(previous_entry['50 knot wind radii max extent (NM)'], 
                                                                                       next_entry['50 knot wind radii max extent (NM)'],
                                                                                       time_between,
                                                                                       time_gap.total_seconds() * SECONDS_PER_HOUR)
    intermediate_step['64 knot wind radii max extent (NM)'] = intermediate_observation(previous_entry['64 knot wind radii max extent (NM)'], 
                                                                                       next_entry['64 knot wind radii max extent (NM)'],
                                                                                       time_between,
                                                                                       time_gap.total_seconds() * SECONDS_PER_HOUR)
    intermediate_step['max wind radius (NM)'] = intermediate_observation(previous_entry['max wind radius (NM)'], 
                                                                         next_entry['max wind radius (NM)'],
                                                                         time_between,
                                                                         time_gap.total_seconds() * SECONDS_PER_HOUR)
    return intermediate_step
                                                                               
    
def intermediate_observation(previous_obs, next_obs, time_between, current_time):
    slope = (next_obs - previous_obs) / time_between
    return previous_obs + current_time * slope
    


# This gets the index or two indicies of the 2 hurricane observation steps bounding a given datetime
def getSynopticBounds(storm, datetime):
    for i in range(len(storm['data']['datetime']) - 1):
        if datetime == storm['data']['datetime'][i]:
            return i
        if datetime == storm['data']['datetime'][i + 1]:
            return i + 1  # here if the datapoint we want is the last one overall
        if storm['data']['datetime'][i] < datetime < storm['data']['datetime'][i + 1]:
            return i, i + 1


# This looks at the status of the storm at the first point it enters a lonlat bounding box, returns None if the path
# never crosses through that box
def statusNearLocation(storm, lon_min=-180, lon_max=180, lat_min=-90, lat_max=90, returnTime=False):
    for i, (lon, lat) in enumerate(zip(storm['data']['lon'], storm['data']['lat'])):
        if lon_min <= lon <= lon_max and lat_min <= lat <= lat_max:
            if not returnTime:
                return storm['data']['system status'][i]
            else:
                return storm['data']['system status'][i], storm['data']['datetime'][i]


def lonLatToXYZ(lon, lat):
    return (np.cos(np.pi / 180 * lon) * np.cos(np.pi / 180 * lat),
            np.sin(np.pi / 180 * lon) * np.cos(np.pi / 180 * lat),
            np.sin(np.pi / 180 * lat))


def XYZtoLonLat(x, y, z):
    return (180 / np.pi * np.arctan2(y, x),
            90 - 180 / np.pi * np.arctan2(np.sqrt(x ** 2 + y ** 2), z))


def twoPointBearing(lon1, lat1, lon2, lat2):
    x_bearing = (np.cos(np.pi / 180 * lat2) *
                 np.sin(np.pi / 180 * (lon2 - lon1)))
    y_bearing = (np.cos(np.pi / 180 * lat1) * np.sin(np.pi / 180 * lat2) -
                 (np.sin(np.pi / 180 * lat1) * np.cos(np.pi / 180 * lat2) *
                  np.cos(np.pi / 180 * (lon2 - lon1))))
    bearing = 180 / np.pi * np.arctan2(x_bearing, y_bearing)
    return bearing


# from https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
# and https://stackoverflow.com/questions/12301071/multidimensional-confidence-intervals

def cov_ellipse(cov, q=None, nsig=None, **kwargs):
    """
    Parameters
    ----------
    cov : (2, 2) array
        Covariance matrix.
    q : float, optional
        Confidence level, should be in (0, 1)
    nsig : int, optional
        Confidence level in unit of standard deviations. 
        E.g. 1 stands for 68.3% and 2 stands for 95.4%.

    Returns
    -------
    width, height, rotation :
         The lengths of two axises and the rotation angle in degree
    for the ellipse.
    """

    if q is not None:
        q = np.asarray(q)
    elif nsig is not None:
        q = 2 * st.norm.cdf(nsig) - 1
    else:
        raise ValueError('One of `q` and `nsig` should be specified.')
    r2 = st.chi2.ppf(q, 2)

    val, vec = np.linalg.eigh(cov)
    width, height = 2 * np.sqrt(val[:, None] * r2)
    rotation = np.degrees(np.arctan2(*vec[::-1, 0]))

    return width, height, rotation


# potential todo: use a kent distribution rather than what I'm doing now

def plotTMinusPosition(hurdat, tMinus, stormType, lon_min=-180, lon_max=180, lat_min=-90, lat_max=90, storm_alpha=0.5,
                       bounds=None, saveFig=None):
    KM_TO_MILES = 0.621371

    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    ax.coastlines()
    if bounds is None:
        bounds = [-100, -30, 10, 90]

    ax.set_extent(bounds, crs=ccrs.Geodetic())
    # ax.set_global()

    lon_list = []
    lat_list = []
    storm_list = []
    status_list = []
    velocity_list = []
    bearing_list = []

    for key in hurdat:
        position, velocity = positionBeforeArrival(hurdat[key], tMinus, lon_min, lon_max, lat_min, lat_max)

        if position is not None:
            lon_list.append(position[0][0])
            lat_list.append(position[0][1])
            storm_list.append(position[1])
        if velocity is not None:
            velocity_list.append(velocity[0])
            bearing_list.append(velocity[1])

    lon_plot = []
    lat_plot = []
    v_plot = []

    for i, (sLon, sLat, storm) in enumerate(zip(lon_list, lat_list, storm_list)):
        if statusNearLocation(storm, lon_min, lon_max, lat_min, lat_max) in stormType:
            ax.plot(sLon, sLat, 'o', transform=ccrs.Geodetic(), alpha=storm_alpha)
            lon_plot.append(sLon)
            lat_plot.append(sLat)
            v_plot.append(velocity_list[i])

    lonLatList = np.stack([lon_plot, lat_plot], axis=1)

    if len(lonLatList) > 2:
        mean = np.mean(lonLatList, axis=0)
        cov = np.cov(lonLatList, rowvar=0)
        mvNorm = st.multivariate_normal(mean, cov=cov)

        S3width, S3height, S3rotation = cov_ellipse(cov, nsig=3)
        ax.add_patch(Ellipse(mean, S3width, S3height, transform=ccrs.Geodetic(),
                             angle=S3rotation, facecolor='White', edgecolor='red', label='3 sigma interval'))

        S2width, S2height, S2rotation = cov_ellipse(cov, nsig=2)
        ax.add_patch(Ellipse(mean, S2width, S2height, transform=ccrs.Geodetic(),
                             angle=S2rotation, facecolor='White', edgecolor='blue', label='2 sigma interval'))

        S1width, S1height, S1rotation = cov_ellipse(cov, nsig=1)
        ax.add_patch(Ellipse(mean, S1width, S1height, transform=ccrs.Geodetic(),
                             angle=S1rotation, facecolor='White', edgecolor='green', label='1 sigma interval'))

        ax.plot([lon_min, lon_max], [lat_min, lat_min], transform=ccrs.PlateCarree(), color='pink')
        ax.plot([lon_min, lon_max], [lat_max, lat_max], transform=ccrs.PlateCarree(), color='pink')
        ax.plot([lon_max, lon_max], [lat_min, lat_max], transform=ccrs.PlateCarree(), color='pink')
        ax.plot([lon_min, lon_min], [lat_min, lat_max], transform=ccrs.PlateCarree(), color='pink',
                label='Bounding Box')

        ax.legend()

        ax.plot(*mean, '+', transform=ccrs.PlateCarree(), color='black')

    ax.set_title(f'Position of systems {tMinus} hours before entering the bounding box')

    plt.figtext(0.5, 0.08,
                f'Mean Velocity of all storms plotted is {np.round(np.mean(v_plot) * KM_TO_MILES, 2)} miles/hr\nwith a standard deviation of {np.round(np.std(v_plot) * KM_TO_MILES, 2)} miles/hr',
                wrap=True, horizontalalignment='center', fontsize=12)

    print(f'{len(lonLatList)} storms plotted')

    if saveFig is not None:
        plt.savefig(saveFig)

    plt.show()
