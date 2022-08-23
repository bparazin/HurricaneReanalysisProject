import pandas as pd
import datetime as dt
    
def cleanInputLine(line):
    return [entry.strip() for entry in line.split(',')]

def readStorm(idcode, name, lines):
    entry_list = pd.DataFrame([cleanInputLine(line) for line in lines])
    result = {}
    result['name'] = name
    result['idcode'] = idcode
    resultData = pd.DataFrame()
    resultData['datetime'] = [dt.datetime(int(day[0:4]),int(day[4:6]),int(day[6:]), int(time[0:2]), int(time[2:4]), tzinfo = dt.timezone.utc) 
           for day, time in zip(entry_list[0], entry_list[1])]
    resultData['record identifier'] = entry_list[2]
    resultData['system status'] = entry_list[3]
    resultData['lat'] = [float(coord[:-1]) * -1 if coord[-1]=='S' else float(coord[:-1]) for coord in entry_list[4]]
    resultData['lon'] = [float(coord[:-1]) * -1 if coord[-1]=='W' else float(coord[:-1]) for coord in entry_list[5]]
    resultData['Max sustained wind (knots)'] = [float(val) if val!='-999' else None for val in entry_list[6]]
    resultData['Minimum Pressure (millibars)'] = [float(val) if val!='-999' else None for val in entry_list[7]]
    #all max wind extents are given as (NE, SE, SW, NW)
    resultData['34 knot wind radii max extent (NM)'] = [(float(val1) if val1!='-999' else None, 
                                                                  float(val2) if val2!='-999' else None, 
                                                                  float(val3) if val3!='-999' else None, 
                                                                  float(val4) if val4!='-999' else None) 
                                                                  for val1, val2, val3, val4 in zip(entry_list[8], entry_list[9], entry_list[10], entry_list[11])]
    resultData['50 knot wind radii max extent (NM)'] = [(float(val1) if val1!='-999' else None, 
                                                                    float(val2) if val2!='-999' else None, 
                                                                    float(val3) if val3!='-999' else None, 
                                                                    float(val4) if val4!='-999' else None) 
                                                                    for val1, val2, val3, val4 in zip(entry_list[12], entry_list[13], entry_list[14], entry_list[15])]
    resultData['64 knot wind radii max extent (NM)'] = [(float(val1) if val1!='-999' else None, 
                                                                    float(val2) if val2!='-999' else None, 
                                                                    float(val3) if val3!='-999' else None, 
                                                                    float(val4) if val4!='-999' else None) 
                                                                    for val1, val2, val3, val4 in zip(entry_list[16], entry_list[17], entry_list[18], entry_list[19])]
    resultData['max wind radius (NM)'] = [val if val !='-999' else None for val in entry_list[20]]
    result['data'] = resultData
    return result

def readHurdat(path, line = 0):
    
    with open(path) as f:
        hurdat_raw = f.readlines()
    
    hurdat_final = {}
    while line < len(hurdat_raw):
        header = cleanInputLine(hurdat_raw[line])
        hurdat_final[header[0]] = readStorm(header[0], header[1], hurdat_raw[line+1:line+1+int(header[2])])['data']
        line += 1 + int(header[2])
    return hurdat_final