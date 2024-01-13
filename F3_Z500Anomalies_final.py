# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 14:13:40 2022

@author: emrus2

This script is used to create a subplot of 500mb geopotential height
for 3 timesteps of each day for three days, colormap
with daily climatology anomalies in contours

bounds are [5300,6050]

UPDATED ON 08/15/2022

- CHANGE HATCH WIDTH AND COLORS
"""

#IMPORT MODULES
import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
    
#GENERATE CUSTOM COLORMAP
def center_colormap(lowlim, highlim, center=0):
    dv = max(-lowlim, highlim) * 2
    N = int(256 * dv / (highlim-lowlim))
    bwr = cm.get_cmap('seismic', N)
    newcolors = bwr(np.linspace(0, 1, N))
    beg = int((dv / 2 + lowlim)*N / dv)
    end = N - int((dv / 2 - highlim)*N / dv)
    newmap = ListedColormap(newcolors[beg:end])
    return newmap

#DEFINE DATES OF ANALYSIS 
#define hour, day, month, year of analysis
year = '2020'
month = '09'
#define integers of date
year_i,month_i = (int(year[2:4]),int(month))

#create subplot for mapping multiple timesteps
fig = plt.figure(figsize=(7,7.7))
fig.suptitle('Z500 Anomalies',fontsize=13,fontweight="bold",y=0.992)

#define subplot labels
sublabel = [False,'a)','e)','i)','b)','f)','j)','c)','g)','k)','d)','h)','l)']


#LOOP THROUGH DESIRED DAYS
n = 1
for day_i in np.arange(6,9,1): #(6,9,1)
    #define string version of days
    i = n
    day = str(day_i)
    if len(day) < 2:
        day = '0' + day
        
    #define MERRA2 location
    filename = (f'MERRA2.inst3_3d_asm_Np.{year}{month}{day}.SUB.nc')
    base_dir = 'I:\\MERRA2\\Daily_and_Subdaily\\500_hPa_Geopotential_Height_3hourly'
    filepath = os.path.join(base_dir,filename)
    
    
    #COLLECT VARIABLE DATA FROM MERRA2 FILE
    #open the netcdf file in read mode
    gridfile = nc.Dataset(filepath,mode='r')
    gridlat = gridfile.variables['lat'][:]
    gridlon = gridfile.variables['lon'][:]
    time = gridfile.variables['time'][:]
    height = np.squeeze(gridfile.variables['H'][:]) #time x height x lat x lon
    gridfile.close()
    
    #DEFINE CLIMATOLOGY LOCATIONS
    clim_dir = 'I:\\Emma\\LaborDayWildfires\\Data\\Z500_Climatology'
    clim_filename = (f'Z500mean_{month}{day}.nc')
    clim_filepath = os.path.join(clim_dir,clim_filename)
    clim_gridfile = nc.Dataset(clim_filepath,mode='r')
    meanheight = np.squeeze(clim_gridfile.variables['H'][:])
    clim_gridfile.close()
    #define std dev. location
    std_filename = (f'Z500std_{month}{day}.nc')
    std_filepath = os.path.join(clim_dir,std_filename)
    std_gridfile = nc.Dataset(std_filepath,mode='r')
    stddevheight = np.squeeze(std_gridfile.variables['H'][:])
    std_gridfile.close()
    #define extremes location
    max_filename = (f'Z500max_{month}{day}.nc')
    max_filepath = os.path.join(clim_dir,max_filename)
    max_gridfile = nc.Dataset(max_filepath,mode='r')
    maxheight = np.squeeze(max_gridfile.variables['H'][:])
    max_gridfile.close()
    min_filename = (f'Z500min_{month}{day}.nc')
    min_filepath = os.path.join(clim_dir,min_filename)
    min_gridfile = nc.Dataset(min_filepath,mode='r')
    minheight = np.squeeze(min_gridfile.variables['H'][:])
    min_gridfile.close()
    
    
    
    #REDUCE VARIABLES TO DESIRED AREA
    #define lat lon restrictions
    latmin = 24.5
    latmax = 75.5
    lonmin = -158.25
    lonmax = -90.75
    #reduce lat
    latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
    latind = np.where(latlims)[0]
    gridlatreduced = gridlat[latind]
    #reduce lon
    lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
    lonind = np.where(lonlims)[0]
    gridlonreduced = gridlon[lonind]
    #reduce height
    heightreduced = height[:,latind,:]
    heightreduced = heightreduced[:,:,lonind]
    #reduce climatology
    mean = meanheight[latind,:]
    mean = mean[:,lonind]
    stddev = stddevheight[latind,:]
    stddev = stddev[:,lonind]
    maximum = maxheight[latind,:]
    maximum = maximum[:,lonind]
    minimum = minheight[latind,:]
    minimum = minimum[:,lonind]    
    
    
    #LOOP THROUGH EACH TIMESTEP
    for timestep in np.arange(0,8,2): #plotting 3,12,21Z timesteps
        zulu = str(3*timestep)
        if len(zulu) < 2:
            zulu = '0' + zulu
    
        #REDUCE VARIABLES TO DESIRED TIMESTEP
        timenew = time[timestep]
        heightnew = heightreduced[timestep,:,:]
        
        #CALCULATE Z500 ANOMALY
        anomaly_meters = heightnew - mean
        anomaly_stddev = anomaly_meters / stddev
        
        #CALCULATE RECORD BREAKING REGIONS
        recordhigh = np.ma.masked_less(heightnew,maximum)
        recordlow = np.ma.masked_greater(heightnew,minimum)
    
        #MAP DESIRED VARIABLE
        #convert lat and lon into a 2D array
        lon, lat = np.meshgrid(gridlonreduced,gridlatreduced) 
        #define area threshold for basemap
        area_thresh = 1E4
        #create equidistant cylindrical projection basemap
        map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
                  urcrnrlon=lonmax,resolution='c',area_thresh=area_thresh)
        xi, yi = map(lon,lat)
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
            bbox=dict(facecolor='1', edgecolor='none', pad=1.5,alpha=0.85,zorder=10))

        #define z-axis bounds for colormap
        lowanom = -7
        highanom = 4
        #create colormap of MERRA2 data
        newmap = center_colormap(lowanom, highanom, center=0)
        colorm = map.pcolor(xi,yi,anomaly_stddev,shading='auto',cmap=newmap,vmin=lowanom,vmax=highanom,zorder=1)
        #define border color and thickness
        border_c = '0.4'
        border_w = 0.4
        map.drawcoastlines(color=border_c,linewidth=border_w,zorder=2)
        map.drawstates(color=border_c, linewidth=border_w,zorder=2)
        map.drawcountries(color=border_c, linewidth=border_w,zorder=2)
        gridlinefont = 9
        parallels = np.arange(35.,71.,15.)
        meridians = np.arange(-145.,-104.,20.)
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
        #add hatching for record high and lows
        mpl.rcParams['hatch.linewidth'] = 0.5
        mpl.rcParams['hatch.color'] = 'k'
        colormhigh = map.pcolor(xi,yi,recordhigh,hatch='///////',lw=0.1,alpha=0.,zorder=5)
        mpl.rcParams['hatch.linewidth'] = 0.8
        mpl.rcParams['hatch.color'] = 'y'
        colormlow = map.pcolor(xi,yi,recordlow,hatch='///////',lw=0.1,alpha=0.,zorder=5)
        
        
        i += 3
    n += 1
    
#CUSTOMIZE SUBPLOT SPACING
# fig.subplots_adjust(left=0.04,right=0.95,bottom=0.083, top=0.94,hspace=0.05, wspace=0.05) #bottom colorbar
fig.subplots_adjust(left=0.028,right=0.948,bottom=0.092, top=0.945,hspace=0.04, wspace=0.04) #bottom colorbar
#fig.add_axis([left,bottom, width,height])
# cbar_ax = fig.add_axes([0.15,0.04,0.7,0.0216]) #bottom colorbar
cbar_ax = fig.add_axes([0.05,0.045,0.876,0.025]) #bottom colorbar
cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(lowanom+1,highanom+1,2),orientation='horizontal')
cbar.ax.tick_params(labelsize=9)
# cbar.set_label(r'\boldmath$\sigma$',fontsize=9.5,labelpad=-0.5,fontweight='bold')
cbar.set_label(r'$\mathbf{\sigma}$',fontsize=9.5,labelpad=-0.5)

#show map
save_dir='I:\\Emma\\LaborDayWildfires\\Figures\\1FinalFigures'
os.chdir(save_dir)
#print(os.getcwd())
plt.savefig('F3_Z500anomalies_rev.png',dpi=300)
plt.show()