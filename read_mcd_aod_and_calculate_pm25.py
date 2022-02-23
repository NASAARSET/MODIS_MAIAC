#!/usr/bin/python
'''
Module: pm25_modis_mcd.py
==========================================================================================
Disclaimer: The code is for demonstration purposes only. Users are responsible to check for accuracy and revise to fit their objective.

Originally Developed by:    Justin Roberts-Pierel & Pawan Gupta, 2015 
Organization:               NASA ARSET


Modified for Cartopy by: Amanda Markert, June 2019
Organization: University of Alabama in Huntsville

Tested on Python Version: 3.7

Purpose: To extract AOD data from a MODIS (MCD19A2) HDFEOS file (or series of files), calculate PM 2.5 from the data, and create a map of the results

Updated: Amanda Markert June 2019
Python Version: 3.7

See the README associated with this module for more information.
==========================================================================================
'''

#import necessary modules
from pyhdf import SD
import numpy as np
import pyproj
import sys
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.colors import LinearSegmentedColormap
from collections import OrderedDict

    
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
    projStr = "+proj={} +lon_0=0 +x_0=0 +y_0=0 +a={} +units=m +no_defs".format(pp,gridInfo['ProjParams'][0])
    proj = pyproj.Proj(projStr)
    gcs = proj.to_latlong()
    #Convert between sinusoidal project to lat/lon coord projection
    lon,lat = pyproj.transform(proj,gcs,xx,yy)

    return lon,lat

# =============================================================================

#loops through all files listed in the text file
for FILE_NAME in fileList:
    FILE_NAME=FILE_NAME.strip()
    user_input=input('\nWould you like to process\n' + FILE_NAME + '\n\n(Y/N)')
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
        
        SDS_NAME = 'Optical_Depth_055' #Enter desired SDS name here
        sds=hdf.select(SDS_NAME)

        #get scale factor for AOD SDS
        attributes=sds.attributes()
        scale_factor=attributes['scale_factor']
        #get valid range for AOD SDS
        range=sds.getrange()
        min_range=min(range)
        max_range=max(range)
        
        #get SDS data
        data=sds.get()
        #get data within valid range
        valid_data=data.ravel()
        valid_data=[x for x in valid_data if x>=min_range]
        valid_data=[x for x in valid_data if x<=max_range]
        valid_data=np.asarray(valid_data)
        #scale the valid data
        valid_data=valid_data*scale_factor
        #find the average
        average=sum(valid_data)/len(valid_data)
        #find the standard deviation
        stdev=np.std(valid_data)
        #print information
        print('\nThe valid range of values is: ',round(min_range*scale_factor,3), ' to ',round(max_range*scale_factor,3),'\nThe average is: ',round(average,3),'\nThe standard deviation is: ',round(stdev,3))
        print('The range of latitude in this file is: ',min_lat,' to ',max_lat, 'degrees \nThe range of longitude in this file is: ',min_lon, ' to ',max_lon,' degrees')
        
        #asks user if they want to set PM2.5 calculation parameters
        user_input=input('\nWould you like to enter a slope and intercept for PM 2.5 calculation?')
        if user_input == 'Y' or user_input == 'y': 
            slope=input('Please enter a slope: ')
            intercept=input('Please enter an intercept: ')
        else:
            #if not, choose the following:
            slope=29.4
            intercept=8.8
        valid_data=data*scale_factor
        pm25=float(slope)*valid_data+float(intercept)
        
        #Asks user if they would like to see a map
        is_map=input('\nWould you like to create a map of this data? Please enter Y or N \n')
        if is_map == 'Y' or is_map == 'y':
            
            #turn fillvalues to NaN
            data=pm25.astype(float)
            data[np.logical_and(data>=0,data <= 12)]=0
            data[np.logical_and(data>12,data <= 35.4)]=1
            data[np.logical_and(data>35.4,data <= 55.4)]=2
            data[np.logical_and(data>55.4,data <= 150.4)]=3
            data[np.logical_and(data>150.4,data <= 250.4)]=4
            data[data>250.4]=5
            data[data < 0] = np.nan
            
            #create the map
            data = np.ma.masked_array(data, np.isnan(data))
            data = np.mean(data,axis=0) #Flattens array into 1-dim by taking mean of 4 orbit SDS values
            
            extent = (min_lon, max_lon, min_lat, max_lat)
            m = plt.axes(projection=ccrs.PlateCarree())
            m.set_extent(extent)
            
            my_cmap=LinearSegmentedColormap.from_list('mycmap', ['green','yellow','orange','red','purple','brown'],6)
            plt.pcolormesh(longitude, latitude, data,cmap=my_cmap, transform=ccrs.PlateCarree())
            plt.clim(0,6)
             
            #create colorbar
            cb = plt.colorbar()
            cb.set_label('AQI Category')
            cb.set_ticks([.5, 1.5,2.5,3.5,4.5,5.5])  # force there to be only 7 ticks
            cb.set_ticklabels(['Good', 'Moderate', 'Unhealthy for \nSensitive Groups','Unhealthy','Very Unhealthy','Hazardous'])  # put text labels on them
            
            m.coastlines()
            grd = m.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
            grd.xlabels_top = None
            grd.ylabels_right = None
            grd.xformatter = LONGITUDE_FORMATTER
            grd.yformatter = LATITUDE_FORMATTER
            plt.autoscale()
            
            #title the plot
            plotTitle=FILE_NAME[:-4]
            plt.title('{0}\n {1}'.format(plotTitle, 'PM 2.5'))
            fig = plt.gcf()
            # Show the plot window.
            plt.show()
            
            
#            #once you close the map it asks if you'd like to save it
#            is_save=str(input('\nWould you like to save this map? Please enter Y or N \n'))
#            if is_save == 'Y' or is_save == 'y':
#                #saves as a png if the user would like5
#                pngfile = '{0}.png'.format(plotTitle)
#                fig.savefig(pngfile)

print('\nAll valid files have been processed.')
