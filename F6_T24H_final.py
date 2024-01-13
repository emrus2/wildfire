# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 14:13:40 2022

@author: emrus2

This script is used to create a subplot of 500mb geopotential height
for 3 timesteps of each day for three days, colormap and contours

bounds = [-31,14]

UPDATED 08/15/2022
"""

#IMPORT MODULES
import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
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
fig.suptitle(u'\u0394T2M and Z500',fontsize=13,fontweight="bold",y=0.992)

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
    day_prev = str(day_i - 1)
    if len(day_prev) < 2:
        day_prev = '0' + day_prev
    #define MERRA2 location
    filename = (f'MERRA2.inst3_3d_asm_Np.{year}{month}{day}.SUB.nc')
    filename1 = (f'MERRA2.tavg1_2d_slv_Nx.{year}{month}{day_prev}.SUB.nc')
    filename2 = (f'MERRA2.tavg1_2d_slv_Nx.{year}{month}{day}.SUB.nc')
    base_dir = 'I:\\MERRA2\\Daily_and_Subdaily\\500_hPa_Geopotential_Height_3hourly'
    base_dir_temp = 'I:\\MERRA2\\Daily_and_Subdaily\\Temperature_at_2_meters_hourly'
    filepath = os.path.join(base_dir,filename)
    filepath1 = os.path.join(base_dir_temp,filename1)
    filepath2 = os.path.join(base_dir_temp,filename2)
    
    
    #COLLECT VARIABLE DATA FROM MERRA2 FILE
    #open the netcdf file in read mode
    #geopotential height file
    gridfile = nc.Dataset(filepath,mode='r')
    height = gridfile.variables['H'][:] #time x height x lat x lon
    gridfile.close()
    #temperature 1 file
    gridfile1 = nc.Dataset(filepath1,mode='r')
    gridlat = gridfile1.variables['lat'][:]
    gridlon = gridfile1.variables['lon'][:]
    time = gridfile1.variables['time'][:]
    temp1 = gridfile1.variables['T2M'][:] #time x height x lat x lon
    gridfile1.close()
    #temperature 2 file
    gridfile2 = nc.Dataset(filepath2,mode='r')
    temp2 = gridfile2.variables['T2M'][:] #time x height x lat x lon
    gridfile2.close()
    
    
    #REDUCE VARIABLES TO DESIRED AREA
    #convert height to a 3D array
    heightsm = np.squeeze(height)    
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
    #reduce height
    heightreduced = heightsm[:,latind,:]
    heightreduced = heightreduced[:,:,lonind]
    #reduce temps
    temp1reduced = temp1[:,latind,:]
    temp1reduced = temp1reduced[:,:,lonind]
    temp2reduced = temp2[:,latind,:]
    temp2reduced = temp2reduced[:,:,lonind]
    
    #REDUCE TEMPERATURE ARRAYS TO 8 TIMESTEPS TO ALIGN WITH GPH DATA 
    temp1sm = temp1reduced[np.arange(0,24,3),:,:]
    temp2sm = temp2reduced[np.arange(0,24,3),:,:]
    
    
    #LOOP THROUGH EACH TIMESTEP
    for timestep in np.arange(0,8,2): #plotting 3,12,21Z timesteps
        zulu = str(3*timestep)
        if len(zulu) < 2:
            zulu = '0' + zulu
    
        #REDUCE VARIABLES TO DESIRED TIMESTEP
        timenew = time[timestep]
        heightnew = heightreduced[timestep,:,:]
        temp1new = temp1sm[timestep,:,:]
        temp2new = temp2sm[timestep,:,:]
        
        #CALCULATE CHANGE IN TEMPERATURE
        tempchange = temp2new - temp1new
        
        #MAP DESIRED VARIABLE
        #convert lat and lon into a 2D array
        lon, lat = np.meshgrid(gridlonreduced,gridlatreduced)  
        #define area threshold for basemap
        area_thresh = 1E4
        #create equidistant cylindrical projection basemap
        map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
                  urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
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
            bbox=dict(facecolor='1', edgecolor='none', pad=1.5,alpha=0.85),zorder=10)
        
        #define z-axis bounds for colormap
        lowlim = -31
        highlim = 14
        lowlimgph = 5300
        highlimgph = 6050
        #create colormap of MERRA2 data
        newmap = center_colormap(lowlim, highlim, center=0)
        colorm = map.pcolor(xi,yi,tempchange,shading='auto',cmap=newmap,vmin=lowlim,vmax=highlim)
        #define border color and thickness
        border_c = '0.4'
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

        #define contour color and thickness
        contour_c = '0.1'
        contour_w = 0.5
        #create contour map
        contourm = map.contour(xi,yi,heightnew,levels=np.arange(lowlimgph,highlimgph+1,70),colors=contour_c,linewidths=contour_w,zorder=5)
        plt.clabel(contourm,levels=np.arange(lowlimgph,highlimgph+1,140),fontsize=6,inline_spacing=1,colors='k',zorder=5)

        i += 3
    n += 1

#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.028,right=0.948,bottom=0.092, top=0.945,hspace=0.04, wspace=0.04) #bottom colorbar
#fig.add_axis([left,bottom, width,height])
cbar_ax = fig.add_axes([0.05,0.045,0.876,0.025]) #bottom colorbar
cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(lowlim+1,highlim+1,10),orientation='horizontal')
cbar.ax.tick_params(labelsize=9)
cbar.set_label(u'\N{DEGREE SIGN}C',fontsize=9.5,labelpad=-0.5,fontweight='bold')

#SHOW MAP
save_dir='I:\\Emma\\LaborDayWildfires\\Figures\\1FinalFigures'
os.chdir(save_dir)
plt.savefig('F6_DeltaT2M_Z500_rev.png',dpi=300)
plt.show()