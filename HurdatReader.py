import pandas as pd
import datetime as dt
import numpy as np


# most of these are a result of the hurdat format, linked in the readme
def cleanInputLine(line):
    return [entry.strip() for entry in line.split(',')]


# Some longitudes are reported as degrees West, even if they are east of prime merdiain which is an odd choice,
# this fixes this up to be (-180, 180]
def handleLon(lon):
    val = float(lon[:-1]) * -1 if lon[-1] == 'W' else float(lon[:-1])
    if val < -180:
        val += 360
    return val


def readStorm(idcode, name, lines):
    entry_list = pd.DataFrame([cleanInputLine(line) for line in lines])
    result = {}
    result['name'] = name
    result['idcode'] = idcode
    resultData = pd.DataFrame()
    resultData['datetime'] = [
        dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:]), int(time[0:2]), int(time[2:4]), tzinfo=dt.timezone.utc)
        for day, time in zip(entry_list[0], entry_list[1])]
    resultData['record identifier'] = entry_list[2]
    resultData['system status'] = entry_list[3]
    resultData['lat'] = [float(coord[:-1]) * -1 if coord[-1] == 'S' else float(coord[:-1]) for coord in entry_list[4]]
    resultData['lon'] = [handleLon(coord) for coord in entry_list[5]]
    resultData['Max sustained wind (knots)'] = [float(val) if val != '-999' else np.nan for val in entry_list[6]]
    statusWithSpeed = []
    for i, maxSpeed in enumerate(resultData['Max sustained wind (knots)']):
        # classify hurricanes by windspeed
        if maxSpeed >= 137 and resultData['system status'][i] == 'HU':
            statusWithSpeed.append('HU5')
        elif maxSpeed >= 113 and resultData['system status'][i] == 'HU':
            statusWithSpeed.append('HU4')
        elif maxSpeed >= 96 and resultData['system status'][i] == 'HU':
            statusWithSpeed.append('HU3')
        elif maxSpeed >= 83 and resultData['system status'][i] == 'HU':
            statusWithSpeed.append('HU2')
        elif maxSpeed >= 64 and resultData['system status'][i] == 'HU':
            statusWithSpeed.append('HU1')
        elif resultData['system status'][i] == 'HU':
            statusWithSpeed.append('HUU')
        else:
            statusWithSpeed.append(resultData['system status'][i])
    resultData['system status'] = statusWithSpeed
    resultData['Minimum Pressure (millibars)'] = [float(val) if val != '-999' else np.nan for val in entry_list[7]]
    # all max wind extents are given as (NE, SE, SW, NW)
    resultData['34 knot wind radii max extent (NM)'] = [np.asarray((float(val1) if val1 != '-999' else np.nan,
                                                         float(val2) if val2 != '-999' else np.nan,
                                                         float(val3) if val3 != '-999' else np.nan,
                                                         float(val4) if val4 != '-999' else np.nan))
                                                        for val1, val2, val3, val4 in
                                                        zip(entry_list[8], entry_list[9], entry_list[10],
                                                            entry_list[11])]
    resultData['50 knot wind radii max extent (NM)'] = [np.asarray((float(val1) if val1 != '-999' else np.nan,
                                                         float(val2) if val2 != '-999' else np.nan,
                                                         float(val3) if val3 != '-999' else np.nan,
                                                         float(val4) if val4 != '-999' else np.nan))
                                                        for val1, val2, val3, val4 in
                                                        zip(entry_list[12], entry_list[13], entry_list[14],
                                                            entry_list[15])]
    resultData['64 knot wind radii max extent (NM)'] = [np.asarray((float(val1) if val1 != '-999' else np.nan,
                                                         float(val2) if val2 != '-999' else np.nan,
                                                         float(val3) if val3 != '-999' else np.nan,
                                                         float(val4) if val4 != '-999' else np.nan))
                                                        for val1, val2, val3, val4 in
                                                        zip(entry_list[16], entry_list[17], entry_list[18],
                                                            entry_list[19])]
    resultData['max wind radius (NM)'] = [float(val) if val != '-999' else np.nan for val in entry_list[20]]
    result['data'] = resultData
    return result


def readHurdat(path, line=0):
    with open(path) as f:
        hurdat_raw = f.readlines()

    hurdat_final = {}
    while line < len(hurdat_raw):
        header = cleanInputLine(hurdat_raw[line])
        hurdat_final[header[0]] = readStorm(header[0], header[1], hurdat_raw[line + 1:line + 1 + int(header[2])])
        line += 1 + int(header[2])
    return hurdat_final


# All bounds are [min, max), and if the storm passes through the bounding box/timebox in at least 1 datapoint,
# it is included Min and max bounds for latlon come from globe geometery, min time is 1700, since Hurdat won't
# go back that far as it is before the United States, and max time is now since you can't see the future. Time bounds
# should be datetime objects, lat and lon should be floats in the range [-180, 180] for lon and [-90, 90] for lat
def trimHurdat(hurdat,
               lonMin=-180, lonMax=180,
               latMin=-90, latMax=91,
               timeMin=dt.datetime(1700, 1, 1, tzinfo=dt.timezone.utc), timeMax=dt.datetime.now(dt.timezone.utc)):
    trimmedHurdat = {}
    for key in hurdat:
        for lon, lat, time in zip(hurdat[key]['data']['lon'], hurdat[key]['data']['lat'],
                                  hurdat[key]['data']['datetime']):
            if lonMax >= lon >= lonMin and latMin <= lat <= latMax and timeMax >= time >= timeMin:
                trimmedHurdat[key] = hurdat[key]
                break
    return trimmedHurdat
