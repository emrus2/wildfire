# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 14:13:40 2022

@author: emrus2

This script is used to create a subplot of 500mb geopotential height
for 3 timesteps of each day for three days, colormap and contours

wind bounds are [0,90]
gph bounds are [9920,11070]

UPDATED ON 1/26/2023
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
    
#DEFINE DATES OF ANALYSIS 
#define hour, day, month, year of analysis
year = '2020'
month = '09'
#define integers of date
year_i,month_i = (int(year[2:4]),int(month))

#create subplot for mapping multiple timesteps
fig = plt.figure(figsize=(7,7.7))
fig.suptitle('V250 and Z250',fontsize=13,fontweight="bold",y=0.992)

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
    filename = (f'MERRA2_401.inst3_3d_asm_Np.{year}{month}{day}.SUB.nc')
    filename_wind = (f'MERRA2.inst3_3d_asm_Np.{year}{month}{day}.SUB.nc')
    base_dir = 'I:\\MERRA2\\Daily_and_Subdaily\\250_hPa_Geopotential_Height_3hourly'
    base_dir_wind = 'I:\\MERRA2\\Daily_and_Subdaily\\East_and_North_wind_components_at_250_hPa'
    filepath = os.path.join(base_dir,filename)
    filepath_wind = os.path.join(base_dir_wind,filename_wind)
    
    
    #COLLECT VARIABLE DATA FROM MERRA2 FILE
    #open the netcdf file in read mode
    #geopotential height file
    gridfile = nc.Dataset(filepath,mode='r')
    # print(gridfile)
    gridlat = gridfile.variables['lat'][:]
    gridlon = gridfile.variables['lon'][:]
    time = gridfile.variables['time'][:]
    height = gridfile.variables['H'][:] #time x height x lat x lon
    gridfile.close()
    #wind file
    gridfile_wind = nc.Dataset(filepath_wind,mode='r')
    ewind = gridfile_wind.variables['U'][:]
    nwind = gridfile_wind.variables['V'][:]
    gridfile_wind.close()
    
    
    #REDUCE VARIABLES TO DESIRED AREA
    #convert height to a 3D array
    heightsm = np.squeeze(height)
    #print(np.shape(heightsm)) #output shows 8x361x576, time,lat,lon
    #reduce wind variables to desired level (2:250hPa)
    # lev = 0
    ewindlev = ewind.squeeze()
    nwindlev = nwind.squeeze()
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
    heightreduced = heightsm[:,latind,:]
    heightreduced = heightreduced[:,:,lonind]
    #reduce winds
    ewindreduced = ewindlev[:,latind,:]
    ewindreduced = ewindreduced[:,:,lonind]
    nwindreduced = nwindlev[:,latind,:]
    nwindreduced = nwindreduced[:,:,lonind]
    
    
    #LOOP THROUGH EACH TIMESTEP
    #define zulu timestep index for mapping and title
    for timestep in np.arange(0,8,2): #plotting 3,12,21Z timesteps
        zulu = str(3*timestep)
        if len(zulu) < 2:
            zulu = '0' + zulu
    
        #REDUCE VARIABLES TO DESIRED TIMESTEP
        timenew = time[timestep]
        heightnew = heightreduced[timestep,:,:]
        ewindnew = ewindreduced[timestep,:,:]
        nwindnew = nwindreduced[timestep,:,:]
        
        #CALCULATE WIND MAGNITUDE USING EAST AND NORTH COMPONENTS
        wspeed = np.sqrt(ewindnew**2 + nwindnew**2)
        
        #MAP DESIRED VARIABLE
        #convert lat and lon into a 2D array
        lon, lat = np.meshgrid(gridlonreduced,gridlatreduced)  
        #define area threshold for basemap
        area_thresh = 1E4
        #create equidistant cylindrical projection basemap
        map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
                  urcrnrlon=lonmax,resolution='c',area_thresh=area_thresh)
        xi, yi = map(lon,lat)
        #add subplot (rows,columns,index)
        print(i)
        ax = fig.add_subplot(4,3,i)
        #add zulu labels on top row
        if day_i == 6:
            ax.set_ylabel(f'{zulu} UTC',fontsize=10,fontweight="bold",labelpad=0.3)
            ax.yaxis.set_label_coords(-0.01,0.5)
        #add sublabels to single plots
        sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
        ax.text(0.0, 1.0, sublabel[i], transform=ax.transAxes + sublabel_loc,
            fontsize=10, verticalalignment='top', fontweight = 'bold',
            bbox=dict(facecolor='1', edgecolor='none', pad=1.5, alpha = 0.85),zorder=10)
        #add date labels on left column
        if timestep == 0:
            ax.set_title(f'{day} SEP',fontsize=11,fontweight="bold",pad=2)
        #define z-axis bounds for colormap
        lowlim = 0
        highlim = 90
        lowlimgph = 9920
        highlimgph = 11070
        #create colormap of MERRA2 data
        colorm = map.pcolor(xi,yi,wspeed,shading='auto',cmap='YlOrRd',vmin=lowlim,vmax=highlim)
        #define border color and thickness
        border_c = '0.4'
        border_w = 0.4
        map.drawcoastlines(color=border_c,linewidth=border_w)
        map.drawstates(color=border_c, linewidth=border_w)
        map.drawcountries(color=border_c, linewidth=border_w)
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
        #define contour color and thickness
        contour_c = '0.1'
        contour_w = 0.5
        #create contour map
        contourm = map.contour(xi,yi,heightnew,levels=np.arange(lowlimgph,highlimgph+1,120),colors=contour_c,linewidths=contour_w,zorder=5)
        plt.clabel(contourm,levels=np.arange(lowlimgph,highlimgph+1,240),fontsize=6.5,inline_spacing=1,colors='b',zorder=5,manual=False)
        i += 3
    n += 1

#CUSTOMIZE SUBPLOT SPACING
# fig.subplots_adjust(left=0.04,right=0.95,bottom=0.083, top=0.94,hspace=0.05, wspace=0.05) #bottom colorbar
fig.subplots_adjust(left=0.028,right=0.948,bottom=0.092, top=0.945,hspace=0.04, wspace=0.04) #bottom colorbar
#fig.add_axis([left,bottom, width,height])
# cbar_ax = fig.add_axes([0.15,0.04,0.7,0.0216]) #bottom colorbar
cbar_ax = fig.add_axes([0.05,0.045,0.876,0.025]) #bottom colorbar
cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(lowlim+5,highlim+1,20),orientation='horizontal')
cbar.ax.tick_params(labelsize=9)
cbar.set_label('m/s',fontsize=9.5,labelpad=-0.5, fontweight='bold')

#SHOW MAP
save_dir='I:\\Emma\\LaborDayWildfires\\Figures\\1FinalFigures'
os.chdir(save_dir)
plt.savefig('F2_Z250_V250_rev.png',dpi=300)
plt.show()