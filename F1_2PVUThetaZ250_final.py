# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 13:48:42 2022

@author: emrus2

This script is used to plot Potential Temperature(Theta) data on the 2PVU surface
to identify Rossby wave breaking, daily subplot

theta bounds are [290,400]
gph bounds are [9770,11150]

UPDATED ON 1/26/2023
"""

#IMPORT MODULES
import netCDF4 as nc #installed using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap


#DEFINE DATES OF ANALYSIS 
year = '2020'
month = '09'
year_i,month_i= (int(year[2:4]),int(month))


#REDUCE VARIABLES TO DESIRED AREA
#define lat lon restrictions
latmin = 0
latmax = 90
lonmin = -180
lonmax = 180

#create subplot for mapping multiple timesteps
fig = plt.figure(figsize=(7,7.7))
#fig.suptitle(r'2PVU $\Theta$ and Z250',fontsize=10.5,fontweight="bold")
fig.suptitle(u'2PVU-\u03b8 and Z250',fontsize=13,fontweight="bold",y=0.992)

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
        
    #IMPORT 250MB GEOPOTENTIAL HEIGHT DATA
    filename_gph = (f'MERRA2_401.inst3_3d_asm_Np.{year}{month}{day}.SUB.nc')
    base_dir_gph = 'I:\\MERRA2\\Daily_and_Subdaily\\250_hPa_Geopotential_Height_3hourly'
    filepath_gph = os.path.join(base_dir_gph,filename_gph)
    gridfile_gph = nc.Dataset(filepath_gph,mode='r')
    gridlat = gridfile_gph.variables['lat'][:]
    gridlon = gridfile_gph.variables['lon'][:]
    time = gridfile_gph.variables['time'][:]
    height = gridfile_gph.variables['H'][:] #time x height x lat x lon
    gridfile_gph.close()
    
    #convert height to a 3D array
    heightsm = np.squeeze(height)
    
    #LOOP THROUGH EACH TIMESTEP
    for timestep in np.arange(0,8,2): #plotting 0,6,12,18Z timesteps
        zulu = str(3*timestep)
        if len(zulu) < 2:
            zulu = '0' + zulu
            
        #DEFINE 2PVU DATA LOCATION
        #set new directory and file name
        filename = (f'PVU2_Theta_{year}{month}{day}_{zulu}Z.nc')
        base_dir = 'I:\\Emma\\LaborDayWildfires\\Data\\2PVU_Data'
        filepath = os.path.join(base_dir,filename)
            
        
        #COLLECT VARIABLE DATA FROM MERRA2 FILE
        #open the netcdf file in read mode
        gridfile = nc.Dataset(filepath,mode='r')
       #gridlat = gridfile.variables['lat'][:]
       # gridlon = gridfile.variables['lon'][:]
        theta = gridfile.variables['THETA'][:] #lat x lon
        gridfile.close()
        
        
        #REDUCE VARIABLES TO DESIRED AREA    
        #reduce lat
        latlims = np.logical_and(gridlat >= latmin, gridlat <= latmax)
        latind = np.where(latlims)[0]
        gridlatreduced = gridlat[latind]
        #lonlims = np.logical_and(gridlon >= lonmin, gridlon <= lonmax)
        #lonind = np.where(lonlims)[0]
        #gridlonreduced = gridlon[lonind]
        gridlonreduced = gridlon
        #thetareduced = theta[latind,:]
        #thetanew = thetareduced[:,lonind]
        thetanew = theta[latind,:]
        heightreduced = heightsm[:,latind,:]
        #heightreduced = heightreduced[:,:,lonind]
        
        #REDUCE VARIABLES TO DESIRED TIMESTEP
        timenew = time[timestep]
        heightnew = heightreduced[timestep,:,:]
            
        #MAP DESIRED VARIABLE
        #convert lat and lon into a 2D array
        lon, lat = np.meshgrid(gridlonreduced,gridlatreduced)
        #define area threshold for basemap
        area_thresh = 1E4
        #create equidistant cylindrical projection basemap
        #map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax,resolution='c')
        map = Basemap(width=12000000,height=9000000,
                    rsphere=(6378137.00,6356752.3142),\
                    resolution='c',area_thresh=area_thresh,projection='lcc',\
                    lat_1=45.,lat_2=55,lat_0=50,lon_0=-115.)
        xi, yi = map(lon,lat)
        #add subplot (rows,columns,index)
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
            fontsize=10, verticalalignment='top', fontweight = 'bold',
            bbox=dict(facecolor='1', edgecolor='none', pad=1.5, alpha=0.85),zorder=10)
        #define z-axis bounds for colormap
        lowlim = 290
        highlim = 500
        lowlimgph = 9770
        highlimgph = 11150    
        #create colormap of MERRA2 data
        colorm = map.pcolor(xi,yi,thetanew,shading='auto',cmap='gist_rainbow_r',vmin=lowlim,vmax=highlim)
        #define border color and thickness
        border_c = '0.4'
        border_w = 0.4
        map.drawcoastlines(color=border_c,linewidth=border_w-0.1)
        map.drawstates(color=border_c, linewidth=border_w-0.1)
        map.drawcountries(color=border_c, linewidth=border_w-0.1)
        gridlinefont = 9
        parallels = np.arange(10.,71.,15.)
        meridians = np.arange(-140.,-80.,25.)
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
        contourm = map.contour(xi,yi,heightnew,levels=np.arange(lowlimgph+150,highlimgph+1,120),colors=contour_c,linewidths=contour_w,zorder=5)

        i += 3
    n += 1

#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.028,right=0.948,bottom=0.092, top=0.945,hspace=0.04, wspace=0.04) #bottom colorbar
#fig.add_axis([left,bottom, width,height])
cbar_ax = fig.add_axes([0.05,0.045,0.876,0.025]) #bottom colorbar
cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(lowlim+10,highlim+1,50),orientation='horizontal')
cbar.ax.tick_params(labelsize=9)
cbar.set_label('K',fontsize=9.5,labelpad=-0.5, fontweight='bold')

    
#SHOW MAP
save_dir='I:\\Emma\\LaborDayWildfires\\Figures\\1FinalFigures'
os.chdir(save_dir)
plt.savefig('F1_2PVUTheta_Z250_rev.png',dpi=300)
plt.show()