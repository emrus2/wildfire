# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 14:13:40 2022

@author: emrus2

This script is used to create a subplot of 2-meter dewpoint temperature
for 3 timesteps of each day for three days, colormap

bounds are [-20,30]

UPDATED ON 08/15/2022
"""

#IMPORT MODULES
import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
from palettable.colorbrewer.diverging import BrBG_10


#DEFINE DATES OF ANALYSIS 
#define hour, day, month, year of analysis
year = '2020'
month = '09'
#define integers of date
year_i,month_i = (int(year[2:4]),int(month))

#create subplot for mapping multiple timesteps
fig = plt.figure(figsize=(7,7.7))
fig.suptitle('DP2M',fontsize=13,fontweight="bold",y=0.992)

#define subplot labels
sublabel = [False,'a)','e)','i)','b)','f)','j)','c)','g)','k)','d)','h)','l)']

#LOOP THROUGH DESIRED DAYS
n = 1
for day_i in np.arange(6,9,1):
    i = n
    #define string version of days
    day = str(day_i)
    if len(day) < 2:
        day = '0' + day
    #define MERRA2 location
    filename = (f'MERRA2_401.tavg1_2d_slv_Nx.{year}{month}{day}.SUB.nc')
    base_dir = 'I:\\MERRA2\\Daily_and_Subdaily\\Dewpoint_Temperature_at_2_meters_hourly'
    filepath = os.path.join(base_dir,filename)
    
    
    #COLLECT VARIABLE DATA FROM MERRA2 FILE
    #open the netcdf file in read mode
    gridfile = nc.Dataset(filepath,mode='r')
    gridlat = gridfile.variables['lat'][:]
    gridlon = gridfile.variables['lon'][:]
    time = gridfile.variables['time'][:]
    dewpoint = gridfile.variables['T2MDEW'][:] #temperature in K
    gridfile.close()
    
    #convert dewpoint temperature to celsius
    dewpoint_C = dewpoint - 273.15
    
    
    #REDUCE VARIABLES TO DESIRED AREA
    #convert height to a 3D array
    dewpointsm = np.squeeze(dewpoint_C)
    #print(np.shape(heightsm)) #output shows 8x361x576, time,lat,lon
    #define lat lon restrictions
    latmin = 29.5
    latmax = 60.5
    lonmin = -142.5
    lonmax = -101.5
    #reduce lat
    latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
    latind = np.where(latlims)[0]
    gridlatreduced = gridlat[latind]
    #reduce lon
    lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
    lonind = np.where(lonlims)[0]
    gridlonreduced = gridlon[lonind]
    #reduce dewpoint
    dewpointreduced = dewpointsm[:,latind,:]
    dewpointreduced = dewpointreduced[:,:,lonind]
    #print(np.shape(dewpointreduced)) #output shows 24x62x64, time,lat,lon
    
    
    #LOOP THROUGH EACH TIMESTEP
    for timestep in np.arange(0,24,6): #plotting 330,1230,2130Z timesteps
        zulu = str(timestep)
        if len(zulu) < 2:
            zulu = '0' + zulu
    
        #REDUCE VARIABLES TO DESIRED TIMESTEP
        timenew = time[timestep]
        dewpointnew = dewpointreduced[timestep,:,:]
    
        #MAP DESIRED VARIABLE
        #convert lat and lon into a 2D array
        lon, lat = np.meshgrid(gridlonreduced,gridlatreduced) 
        #define area threshold for basemap
        area_thresh = 1E4
        #create equidistant cylindrical projection basemap
        map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
                  urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
        xi, yi = map(lon,lat)
        #add figure to subplot
        print(i)
        ax = fig.add_subplot(4,3,i)
        #add zulu labels on top row
        if day_i == 6:
            ax.set_ylabel(f'{zulu} UTC',fontsize=10,fontweight="bold",labelpad=0.3)
            ax.yaxis.set_label_coords(-0.01,0.5)
        #add date labels on left column
        if timestep == 0:
            ax.set_title(f'{day} SEP',fontsize=11,fontweight="bold",pad=2)
        #add sublabels to single plots
        sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
        ax.text(0.0, 1.0, sublabel[i], transform=ax.transAxes + sublabel_loc,
            fontsize=10, verticalalignment='top', fontweight='bold',
            bbox=dict(facecolor='1', edgecolor='none',alpha=0.85, pad=1.5),zorder=10)
        #define z-axis bounds for colormap
        lowlim = -20
        highlim = 30
        #create colormap of MERRA2 data with contours
        colorm = map.pcolor(xi,yi,dewpointnew,shading='auto',cmap=BrBG_10.mpl_colormap,vmin=lowlim,vmax=highlim)
        #define border color and thickness
        border_c = '0.3'
        border_w = 0.4
        map.drawcoastlines(color=border_c,linewidth=border_w)
        map.drawstates(color=border_c, linewidth=border_w)
        map.drawcountries(color=border_c, linewidth=border_w)
        gridlinefont = 9
        parallels = np.arange(35.,60.,10.)
        meridians = np.arange(-134.,-109.,12.)
        if i == 3 or i == 6 or i == 9:
            map.drawparallels(parallels, labels=[0,1,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
            map.drawmeridians(meridians,color=border_c,linewidth=border_w)
        elif i == 10 or i == 11:
            map.drawparallels(parallels, color=border_c,linewidth=border_w)
            map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
        elif i == 12:
            map.drawparallels(parallels, labels=[0,1,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
            map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
        else:
            map.drawparallels(parallels, color=border_c,linewidth=border_w)
            map.drawmeridians(meridians,color=border_c,linewidth=border_w)

        i += 3
    n += 1

#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.028,right=0.948,bottom=0.092, top=0.945,hspace=0.04, wspace=0.04) #bottom colorbar
#fig.add_axis([left,bottom, width,height])
cbar_ax = fig.add_axes([0.05,0.045,0.876,0.025]) #bottom colorbar
cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(lowlim,highlim+1,10),orientation='horizontal')
cbar.ax.tick_params(labelsize=9)
cbar.set_label(u'\N{DEGREE SIGN}C',fontsize=9.5,labelpad=-0.9,fontweight='bold')

#show map
save_dir='I:\\Emma\\LaborDayWildfires\\Figures\\1FinalFigures'
os.chdir(save_dir)
#print(os.getcwd())
plt.savefig('F5_DP2M_rev.png',dpi=300)
plt.show()