# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 14:00:14 2022

@author: emrus2

This script is used to plot HYSPLIT backward air parcel trajectory data from a
downloaded trajectory results .txt file using Basemap and high-res topography data
"""

"""
Adjust the min and max bounds for the height (which should be adjusted for terrain height)

Terrain heights can be run together but arrival heights need to be adjusted. Fix this...
"""
#%% IMPORTING PACKAGES
#IMPORT REQUIRED MODULES
import os
import pandas as pd
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.transforms as mtransforms
from matplotlib.lines import Line2D
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap

#%% UPLOADING TOPOGRAPHY DATA -only complete on initial HYSPLIT mapping, takes a long time
#define ETOPO data filepath
etopofile = 'ETOPO1_Bed_g_gmt4.grd'
etopodir = 'I:\\Emma\\LaborDayWildfires\\Data\\Topo_Data\\ETOPO1_Bedrock'
etopopath = os.path.join(etopodir,etopofile)


#COLLECT VARIABLE DATA FROM ETOPO FILE
gridfile = nc.Dataset(etopopath)

#grab etopo data from dataset
gridlon = gridfile.variables['x'][:]
gridlat = gridfile.variables['y'][:]
etopo = gridfile.variables['z'][:] #lat x lon
gridfile.close()

#define resolution of topography data
interval = 1 #the higher the interval, the lower the resolution
skip1 = slice(None,None,interval)
skip = (slice(None, None, interval), slice(None, None, interval))
gridlatreduced = gridlat[skip1]
gridlonreduced = gridlon[skip1]
etoporeduced = etopo[skip]

#%% REST OF ANALYSIS 

#create function for string spacing of locations
def string_spacing(string):
    new_word = ''
    i = 0
    for l in string:
        if l.isupper() and i != 0:
            new_word += ' '
        new_word += l    
        i += 1
    return new_word
                
#define various locations
locations = ['Detroit','Medford','Malden','GrandJunction','LincolnCity','Denver']
#locations = [locations[0]]

#create dictionary for state of location
d_states = {'Detroit':'OR','Medford':'OR','Malden':'WA','GrandJunction':'CO','LincolnCity':'OR','Denver':'CO'}
    
#create dictionaries of bounds for plotting meterological variables
    #mins, maxs, intervals, and tickstarts for each of the three plots below the map
    #[height,potential temp, ambient temp]
d_ymins = {'Detroit':[-300,-28,-28],'Medford':[-407,-45,-45],'Malden':[-220,-15,-15],'GrandJunction':[-300,-7,-7],'LincolnCity':[-387,-40,-40],'Denver':[-300,-4,-4]}
d_ymaxs = {'Detroit':[5500,43,43],'Medford':[7600,50,50],'Malden':[4000,45,45],'GrandJunction':[5600,55,55],'LincolnCity':[7200,40,40],'Denver':[5200,55,55]}
d_intervals = {'Detroit':[1500,20,20],'Medford':[2000,20,20],'Malden':[1000,15,15],'GrandJunction':[1500,15,15],'LincolnCity':[2000,20,20],'Denver':[1500,15,15]}
d_tickstarts = {'Detroit':[0,-20,-20],'Medford':[0,-40,-40],'Malden':[0,-15,-15],'GrandJunction':[0,0,0],'LincolnCity':[0,-40,-40],'Denver':[0,0,0]}
d_ylabelpad = {'Detroit':[3,6,6],'Medford':[3,6,6],'Malden':[3,6,6],'GrandJunction':[3,13,13],'LincolnCity':[3,6,6],'Denver':[3,13,13]}
#define subplot labels
sublabel = ['a)','b)','c)']

#determine day, month, year of analysis
year = '2020'
month = '09'
days = ['06','07','08']
#days = ['06']

#define trajectory colors, symbols, and arrival hours not including 0th index
colors = ['null','r','b','g','c','m','darkorange']
markers = ['null','>','x','d','s','.','p']
arriv_hours = ['null', '20', '16', '12', '08', '04', '00']
#define marker density and size
markdensity = 6
marksize = 4

#%% FIGURE MAKING

for location in locations:
    #create subplot, determine sizes
    fig = plt.figure(figsize=[14.2, 9])
    location_spaced = string_spacing(location)
    TITLE = f'{location_spaced}, {d_states[location]}'
    fig.suptitle(TITLE,fontweight='bold',size=16, y=0.9875)
    spec = fig.add_gridspec(ncols=3, nrows=6) #spec=row,column
    
    n_column = 0
    for day in days:
        n_row = 3
        n_var = 0
        #convert day into integer
        day_int = int(day)
    
        #define HYSPLIT variables of interest to plot on subplot 2
        variables = ['TerrainHeight','PotentialTemperature','AmbientTemperature']
        #variables = ['TerrainHeight']
        names = ['H (m)',u'\u03B8 (\N{DEGREE SIGN}C)',u'T (\N{DEGREE SIGN}C)']
    
        for var in variables:
            #define HYSPLIT data filepath
            filename = (f'hysplit_traj_{var}_{year}{month}{day}.txt')
            base_dir = f'I:\\Emma\\LaborDayWildfires\\Data\\Hysplit_Trajectory_Data\\{location}'
            filepath = os.path.join(base_dir,filename)
    
            df = pd.read_csv(filepath, sep='\s+', skiprows = 19,names = ['TrajNum','MetGrid',\
            'Year','Month','Day','Hour','Minute','ForecastHour','TrajAge','Lat','Lon','Height','Pressure','MetVar']) #separation is one or more spaces
            
            if var == 'TerrainHeight':
                #calculate height above sea level using terrain height and height
                df['MetVar'] += df['Height']
                #set ylabel padding
                ylabelpad = 3
                #adjust arrival times for trajectories
                arrival = np.arange(16,-1,-4) #(0,4,12,16) UTC
                for i in np.arange(1,6): #change arrival times for rows 2-6
                    df.at[i,'Day'] = day_int
                    df.at[i,'Hour'] = arrival[i-1]
        
                #CREATE MAP TO PLOT TRAJECTORIES
                #determine maximum and minimum lat/lon values of dataframe
                latmin = round((min(df['Lat']) - 3))
                latmax = round((max(df['Lat']) + 3))
                lonmin = round((min(df['Lon']) - 5))
                lonmax = round((max(df['Lon']) + 5))
        
                #determine if map is not square
                lats = latmax - latmin
                lons = lonmax - lonmin
                diff = abs(lons-lats)
                #adjust to create square map
                if lats > lons:
                    lonmin -= diff/2
                    lonmax += diff/2
                if lats < lons:
                    latmin -= diff/2
                    latmax += diff/2
                lats = latmax - latmin
                #print(latmin,latmax)
                #print(lonmax,lonmin)
                
                #adjust bounds for aesthetic of specific maps
                if location == "Malden" and n_column == 2 \
                or location == "GrandJunction" and n_column != 0 \
                or location == "Denver" and n_column == 1 \
                or location == "Medford" and n_column == 2 \
                or location == "LincolnCity" and n_column != 1 \
                or location == "Detroit":
                    latmin -= 1
                    latmax -= 1
                
                #reduce the latitude, longitude, and etopo to Western US
                latlims = np.logical_and(gridlatreduced > latmin, gridlatreduced < latmax)
                latind = np.where(latlims)[0]
                gridlatnew = gridlatreduced[latind]
                
                lonlims = np.logical_and(gridlonreduced > lonmin, gridlonreduced < lonmax)
                lonind = np.where(lonlims)[0]
                gridlonnew = gridlonreduced[lonind]
                
                etoponew = etoporeduced[latind,:]
                etoponew = etoponew[:,lonind]
    
                #convert lat and lon into a 2D array
                lon, lat = np.meshgrid(gridlonnew,gridlatnew)
    
                #SUBPLOT 1
                #create subplot 1, define gridspec
                ax1 = fig.add_subplot(spec[0:3,n_column])
                #set subplot title
                ax1.set_title(f'{days[n_column]} SEP',fontsize=12,fontweight='bold',pad=0.4)
                #create equidistant cylindrical projection basemap for subplot 1
                area_thresh = 10000
                map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh) #low resolution
                xi, yi = map(lon,lat)
    
                #create a colormap of topography data
                vmin,vmax = (-2500,3750)
                colorm = map.pcolor(xi,yi,etoponew,shading='auto',cmap='terrain',zorder=1,vmin=vmin,vmax=vmax)
                #plt.colorbar(colorm)

                #define border color and thickness
                border_c = '0.1'
                border_w = 0.5
    
                #map underlying map, with  parallel labels on the left, and meridian labels on the bottom
                map.drawlsmask(land_color='none',ocean_color='paleturquoise',zorder=2)
                try:
                    map.drawcoastlines(color=border_c,linewidth=border_w,zorder=3)
                except:
                    pass
                map.drawstates(color=border_c, linewidth=border_w,zorder=3)
                map.drawcountries(color=border_c, linewidth=border_w,zorder=3)
                                
                #define gridline interval for parallels and longitudes
                if lats > 21:
                    gridline_int = 10
                else:
                    gridline_int = 5
                map.drawparallels(np.arange(30,latmax,gridline_int), labels=[0,1,0,0], fontsize=10,color=border_c, linewidth=border_w,zorder=3)
                map.drawmeridians(np.arange(-95,lonmin+1,-gridline_int), labels=[0,0,0,1], fontsize=10,color=border_c, linewidth=border_w,zorder=3)

                #Plot trajectories on map
                #create legend elememt list
                legend_elements = []
                for traj, track in df.groupby(['TrajNum']):
                    latitude = track.Lat.values
                    longitude = track.Lon.values
                    xi,yi = map(longitude, latitude)
                    ax1.plot(xi,yi,'-',linewidth=1.5,color=colors[traj],marker=markers[traj],markevery=markdensity,markersize=marksize,zorder=4)
                    #create legend element for each trajectory
                    legend_elem = Line2D([0], [0], color=colors[traj], marker=markers[traj], label = f'{arriv_hours[traj]} UTC')
                    legend_elements.insert(0,legend_elem)
                #plot destination lat,lon
                dest_lon, dest_lat = (df.at[0,'Lon'],df.at[0,'Lat'])
                x, y = map(dest_lon,dest_lat)
                ax1.plot(x,y,marker='*', color='k',markersize=6,zorder=5)
                
                #add sublabels to single plots
                sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
                ax1.text(0.005, 0.995, sublabel[n_column], transform=ax1.transAxes + sublabel_loc,
                    fontsize=10, verticalalignment='top', 
                    bbox=dict(facecolor='1', edgecolor='none', pad=2),zorder=10)
                
                #add legend for arrival times
                if n_column == 0:
                    ax1.legend(handles=legend_elements, loc='lower left', title='Arrival Time',title_fontproperties = {'weight':'bold'})
                   
            #convert temperatures to celsius
            else:             
                df['MetVar'] -= 273.15
                #set ylabel padding
                ylabelpad = 4
                
            #convert a reduced 2-num year to full 4-num year
            year_prefix = year[0:2]     #first two numbers of year from year variable
            #concatenate first two numbers to last two numbers
            df['Year'] = year_prefix + df['Year'].astype(str)
            #convert date information to full singular date
            df['Dates'] = pd.to_datetime(df[['Year','Month','Day','Hour']])
                    
            #SUBPLOT 2
            #create subplot 2, define gridspec
            ax2 = fig.add_subplot(spec[n_row,n_column])
    
            #PLOT TRAJECTORY HEIGHT DATA ON SUBPLOT 2
            for traj, track in df.groupby(['TrajNum']):
                trackvar = track.MetVar.values
                dates = track.Dates.values
                ax2.plot(dates,trackvar,'-',linewidth=1,color=colors[traj],marker=markers[traj],markevery=markdensity,markersize=marksize)
                
            #plot trajectory arrivals
            for i in np.arange(6):
                arriv_time = df.at[i,'Dates']
                arriv_var = df.at[i,'MetVar']
                ax2.plot(arriv_time,arriv_var,'*',color='k',markersize=marksize)
        
            #add gridlines
            ax2.grid(which='major',axis='y',linestyle=':',color=border_c,linewidth=border_w)
            #define dictionaries for values
            ymins = d_ymins[location]
            ymaxs = d_ymaxs[location]
            intervals = d_intervals[location]
            tickstarts = d_tickstarts[location]
            ylabelpad = d_ylabelpad[location]
            #set yaxis label and ticks
            ax2.yaxis.set_label_position('right')
            ax2.yaxis.tick_right()
            bottom_lim, top_lim, interval,tickstart = (ymins[n_var],ymaxs[n_var],intervals[n_var],tickstarts[n_var])
            ax2.set_ylim(bottom=bottom_lim,top=top_lim)
            yticks = np.arange(tickstart,top_lim+1,interval)  
            ax2.set_yticks(yticks)
            
            #create ylabels on rightmost side only
            if n_column == 2:
                ax2.set_ylabel(names[n_var],fontweight='bold',labelpad=ylabelpad[n_var])
        
            #set x axis label and ticks
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H'))    
            
            #hide tick labels for interior plots
            if n_var != 2:
                plt.gca().tick_params(axis='x',label1On=False)
            else:
                ax2.set_xlabel(f'          {days[n_column]} SEP',fontsize=8)
            ax2.tick_params(axis='both',which='both',labelsize=8)
                
            #add to counters
            n_row += 1
            n_var += 1
            
        n_column += 1
        
    #adjust location of subplots, positions are fractions of the figure width and height
    fig.subplots_adjust(top=0.94, bottom = 0.0475, right = 0.957, left = 0.043, wspace=.12)
    
    #SAVE FIGURE
    save_dir='I:\\Emma\\LaborDayWildfires\\Figures\\HYSPLIT'
    os.chdir(save_dir)
    plt.savefig(f'TRAJECTORIES_{location}_legends.png',dpi=300)
    plt.show()