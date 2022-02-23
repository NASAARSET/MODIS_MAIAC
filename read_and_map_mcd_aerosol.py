'''
Module: read_and_map_mod_aerosol.py
==========================================================================================
Disclaimer: The code is for demonstration purposes only. Users are responsible to check for accuracy and revise to fit their objective.

Originally Developed by:    Justin Roberts-Pierel & Pawan Gupta, 2015 
Organization:               NASA ARSET

Modified for Cartopy & MCD19A2 by: Amanda Markert, June 2019
Organization: University of Alabama in Huntsville

Tested on Python Version: 3.7

Purpose: To read from a MODIS MAIAC HDFEOS dataset and generate a plot of AOD data

See the README associated with this module for more information.
==========================================================================================
'''

import pyproj
import numpy as np
import numpy.ma as ma
from collections import OrderedDict
import sys
from pyhdf import SD
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


# =============================================================================
# Inputs
#This uses the file "fileList.txt", containing the list of files, in order to read the files

try:
    fileList=open('fileList.txt','r')
except:
    print('Did not find a text file containing file names (perhaps name does not match)')
    sys.exit()
    
# =============================================================================
# Funtstions that pull metadata from HDFEOS file and extract lat/lon coordinates information
# to plot data.

def parse_hdfeos_metadata(string):
  out = OrderedDict()
  lines = [i.replace('\t','') for i in string.split('\n')]
  i = -1
  while i<(len(lines))-1:
      i+=1
      line = lines[i]
      if "=" in line:
          key,value = line.split('=')
          if key in ['GROUP','OBJECT']:
              endIdx = lines.index('END_{}={}'.format(key,value))
              out[value] = parse_hdfeos_metadata("\n".join(lines[i+1:endIdx]))
              i = endIdx
          else:
              if ('END_GROUP' not in key) and ('END_OBJECT' not in key):
                   try:
                       out[key] = eval(value)
                   except NameError:
                       out[key] = str(value)
  return out

def construct_coords(ds,grid='GRID_1'):
    attrs = ds.attributes()
    metadata = parse_hdfeos_metadata(attrs['StructMetadata.0'])

    gridInfo = metadata['GridStructure'][grid]

#    gridName = gridInfo['GridName']

    x1,y1 = gridInfo['UpperLeftPointMtrs']
    x2,y2 = gridInfo['LowerRightMtrs']
    yRes = (y1-y2)/gridInfo['YDim']
    xRes = (x1-x2)/gridInfo['XDim']

    #setting up coordinate grids along x and y axis
    x = np.arange(x2,x1,xRes)
    y = np.arange(y2,y1,yRes)[::-1]
    #set up 2D grid for plotting
    xx,yy = np.meshgrid(x,y)
    #get projection information
    if 'soid' in gridInfo['Projection'].lower():
        pp = 'sinu'
    else:
        pp = gridInfo['Projection'].lower()

    #formating projection name from metadata to pyproj for sinusoidal projection
    projStr = "+proj={} +lon_0=0 +x_0=0 +y_0=0 +a={} +units=m +no_defs".format(
      pp,gridInfo['ProjParams'][0])
    proj = pyproj.Proj(projStr)
    gcs = proj.to_latlong()
    #Convert between sinusoidal project to lat/lon coord projection
    lon,lat = pyproj.transform(proj,gcs,xx,yy)

    return lon,lat

# =============================================================================
# Extract and plot the data using the spatial extent and projection information from the functions

for FILE_NAME in fileList:

    user_input=input('Would you like to process\n' + FILE_NAME + '\n\n(Y/N)')

    if(user_input == 'N' or user_input == 'n'):
        continue
    else:
        # open the hdf file for reading
        hdf=SD.SD(FILE_NAME)
        attrs = hdf.attributes()
        
        longitude,latitude = construct_coords(hdf)
        min_lat=latitude.min()
        max_lat=latitude.max()
        min_lon=longitude.min()
        max_lon=longitude.max()

        #Selecte values from file for plot
        SDS_NAME = 'Optical_Depth_055' #Enter desired SDS name here
        sds = hdf.select('Optical_Depth_055') 
        data=sds.get()
        datamsk=ma.masked_where(data < 0, data) #Need to remove fill values <0 because they will affect taking mean of orbits
        sdsgrid = np.mean(datamsk,axis=0) #Flattens array into 1-dim by taking mean of 4 orbit SDS values

        sdsAttrs=sds.attributes()
        SDS_NAME = sdsAttrs['long_name'] #Extract SDS longname for plot title
        scale_factor = sdsAttrs['scale_factor']
        sdsgrid = sdsgrid * scale_factor #Scale factor to correct SDS values
    
        #Plot the data
        m = plt.axes(projection=ccrs.PlateCarree())
        extent = (min_lon, max_lon, min_lat, max_lat)
        m.set_extent(extent)
        plt.pcolormesh(longitude, latitude, sdsgrid, cmap=plt.cm.jet, transform=ccrs.PlateCarree())
        
        grd = m.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
        grd.xlabels_top = None
        grd.ylabels_right = None
        grd.xformatter = LONGITUDE_FORMATTER
        grd.yformatter = LATITUDE_FORMATTER
        plotTitle=FILE_NAME[:-4]
        plt.title('{0}\n {1}'.format(plotTitle, SDS_NAME))
        m.coastlines()
        plt.autoscale()
        cb = plt.colorbar() #create colorbar
        cb.set_label('AOD')

        plt.show() #Show the plot
        
print('\nAll valid files have been processed.')