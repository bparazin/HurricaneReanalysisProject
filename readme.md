# Hurdat2 visualization tools

## Required Dependencies
* [Matplotlib](https://matplotlib.org/)
* [Cartopy](https://scitools.org.uk/cartopy/docs/latest/)
* [Geopy](https://geopy.readthedocs.io/en/stable/)
* [Scipy](https://scipy.org/)
* [Numpy](https://numpy.org/)
* [Pandas](https://pandas.pydata.org/)
* [Jupyter Notebooks](https://jupyter.org/)

## How to Use
Clone this github repo onto your computer and download one of the two
HURDATS from the [NOAA database](https://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html),
saving it as a .txt file. The easiest way to do this is open the HURDAT in your
browser, select all with CTRL-A, copy it with CTRL-C and paste it into a .txt file of
your choice, making sure to delete any blank lines added at the start and end of the
file. An example file, accessed August 24th, 2022 can be found in data/HURDAT.txt

### Reading the HURDAT

To read the HURDAT file, use the methods readHurdat() and trimHurdat:

`readHurdat(path, line=0)`

> path: str, the path to the location of the hurdat.txt file on your machine

> line: int, optional, the first line of the hurdat to start reading at

Outputs the hurdat as a large dict, with each key corresponding to the basin, ATCF cyclone number and 
year of each cyclone, forming a unique identifier string. Each storm is itself a dict with
the keys name, idcode and data. Name is the storm's name if it was given one, idcode is the previously
mentioned identifier string, and data is the information about the storm contained in the hurdat.
Details about the information formatting can be found at the HURDAT2 format description 
[here](https://www.aoml.noaa.gov/hrd/hurdat/hurdat2-format.pdf). 

The only major change is that hurricanes are also
classified by maximum sustained wind speed at each timestamp, changing the system status HU to HU1 through HU5
corresponding to the hurricane status, as well as HUU, meaning "hurricane unclassified" for older storms without
reliable maximum sustained wind speed measurements. 

Wind radii max extent is reported in Nautical Miles as a tuple 
formatted as (Northeast quadrant, southeastern quadrant, southwestern quadrant, northwestern quadrant)

`trimHurdat(hurdat, 
               lonMin = -180, lonMax = 180, 
               latMin = -90, latMax = 91, 
               timeMin = dt.datetime(1700, 1, 1, tzinfo = dt.timezone.utc), timeMax = dt.datetime.now(dt.timezone.utc))`

> hurdat: Dict, hurdat as returned by readHurdat

> lonMin, lonMax: Float in range [-180, 180], the min and max longitudes the storm passes through

> latMin, latMax: Float in range [-90, 90], the min and max latitudes the storm passes through

> timeMin, timeMax: Datetime objects with a timezone, the min and max times the storm can pass through the bounding box

Returns a hurdat restricted only to storms that pass through the bounding box defined during the timezone given

## Plotting Storms

There are 2 very programs for plotting storm positions, plot paths.ipynb and PlotPaths.py. Both contain the same code
and only differ in how you use them (importing plotPaths.py or directly running the jupyter notebook).
In both, the key function is plotTMinusPosition, which plots all the storms that pass through a bounding box a given
amount of time before they first enter that bounding box.

`plotTMinusPosition(hurdat, tMinus, stormType, lon_min=-180, lon_max=180, lat_min=-90, lat_max=90, 
storm_alpha=0.5, bounds = None)`

>hurdat: dict, hurdat as returned by readHurdat

> tMinus: float, how many hours before the storms enter the given bounding box to plot them. A negative value will plot
> storms |tMinus| hours *after* entering the bounding box

> stormType: array-like, this is a collection of strings corresponding to the storm system status when first entering
> the bounding box you want to restrict your plotting to

> lon_min, lon_max, lat_min, lat_max: floats defining the bounding box that all storms plotted pass through

> storm_alpha: float in [0,1], the opacity of all storms being plotted

> bounds: [lon_min, lon_max, lat_min, lat_max] floats defining the geographical bounds of the image

> saveFig: None or string, None to not save the picture created, or a string of the path and filename to save the image
> produced under

This does not return anything, but shows a cartopy plot of all the storms that pass through the bounding box, tMinus
hours before entering it, along with a black cross at the mean position and 3 ellipses representing the 1, 2 and 3 sigma
confidence areas of where storms that enter that bounding box are found. This is calculated using a multivariate normal
distribution, which assumes the points are distributed on a flat surface, rather than a spherical surface, meaning
there is some error there, but we determined that the development cost of using a kent distribution was higher than the
gain caused by the increased precision. The velocity data given at the bottom of the image is found by using the average
velocity between the current storm position, and the position at the next observation unless the current storm position
is the last one recorded, in which case it is the average velocity the storm sustained between the previous and current
positions.

### Bar Chart

This is a way of visualizing the temporal distribution of the storms which pass through a certain bounding box, and is 
found in BarChart.ipynb

`monthlyBarChartByGroup(hurdat, lon_min=-180, lon_max=180, lat_min=-90, lat_max=90, saveFig = None)`

> hurdat: dict, hurdat as read by readHurdat

> lon_min, lon_max, lat_min, lat_max: floats, defining the bound box that restricts the storm list to ones that pass 
> through it

> safeFig: None or string. None to not save the figure, and a string of the path/filename under which to save the figure