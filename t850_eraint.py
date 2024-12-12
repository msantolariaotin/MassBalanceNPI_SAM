'''
for i in DJF MAM JJA SON;do
python t850_eraint.py ${i}
done
'''
import sys
sys.path.insert(1, '/home/msantolaria/Documents/MyPythonLibrary/ClimAnag/')
import climbasis as climb
from climbasis import *
import domain as dom
import myplot
import glob
import cartopy
import cartopy.feature as cpf
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from shapely import geometry
from collections import namedtuple
from shapely.geometry.polygon import LinearRing

# In[3]:


source='/home/msantolaria/Documents/MyResearch/GabiSAM/UBSam/'
resultsDir=source + 'Results/'
plotsDir=source + 'Plots/'


#resultsDir='/media/maria/MARIAEXT2/WinterTrendsObs/'
sourceData='/home/msantolaria/Documents/Data/'


# In[5]:


iyr=1980
imon=1
fmon=12
fyr=2014
mapa=True

domain='SHWest'
season=sys.argv[1]
variable='t850'
units='K'
model=str('eraint')
fileName ='t850_mon_eraint_1979-2015.nc'
ds = xr.open_dataset(sourceData+model+'/'+fileName)['t']

ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180))
ds = ds.roll(longitude=int(len(ds['longitude']) / 2), roll_coords=True)
#field=dom.field_dom(ds,domain)

lat,lon=climb.latlon(ds)
ylat=ds.coords['latitude']
xlon=ds.coords['longitude']

field=ds.sel(latitude=slice(-20,-90),longitude=slice(-180,0))
latS,latN,lonW,lonE,latlim,lonlim=dom.coord_domain(domain)
#field=dom.field_sel_dom(ds,domain)
print(ds.units)

if season[0]=='m':
    rmon=int(season.split('mon')[1])
    vals,anoms=climb.monthly_selection(field,rmon,iyr,fyr)
else:
    if season[0]=='D':
        print('Winter',iyr,'-',iyr+1,fyr-1,'-',fyr)
        vals,anoms=climb.seasonal_selection(field,season,6,iyr,6,fyr)
    else:
        vals,anoms=climb.seasonal_selection(field,season,1,iyr,12,fyr)


anoms_detrend=climb.detrend_dim(vals, 'time', deg=1)

clim=vals.mean('time')
std=vals.std('time')
std_det=anoms_detrend.std('time')
#--------------------------------
subdomain=domain

'''
##Clim-----------------
#clevs=np.arange(0,14.25,0.25)
clevs=np.arange(-30,32,2)
color='bwr'
exp='obs'
decomp='original'
lat,lon=climb.latlon(clim)
figclim=myplot.oneplot_ds(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevs, color=color, units=units, subdomain=subdomain, ds=clim-273.15, xlon=clim[lon], ylat=clim[lat], season=season, iyr=iyr, fyr=fyr,extent=mapa)
ofileC='clim_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figclim.savefig(plotsDir+ofileC+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileC)
plt.show()
# In[21]:


##STD-----------------
clevs=np.arange(0,3.5,0.1)
color='rainbow'
exp='obs'
decomp='original'
lat,lon=climb.latlon(std)
figstd=myplot.oneplot_ds(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevs, color=color, units=units, subdomain=subdomain, ds=std, xlon=clim[lon], ylat=clim[lat], season=season, iyr=iyr, fyr=fyr,extent=mapa)
ofileS='anoms_std_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figstd.savefig(plotsDir+ofileS+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileS)

###MEAN DET
clevs=np.arange(0,3.5,0.1)
color='rainbow'
exp='obs'
decomp='original-detrend'
lat,lon=climb.latlon(std_det)
figstd_det=myplot.oneplot_ds(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevs, color=color, units=units, subdomain=subdomain, ds=std_det, xlon=clim[lon], ylat=clim[lat], season=season, iyr=iyr, fyr=fyr,extent=mapa)
ofileS='anoms_det_mean_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+subdomain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figstd_det.savefig(plotsDir+ofileS+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileS)
plt.show()
'''
#-----------------------------------
par=climb.trend_vect(anoms.time,anoms,'time')
trend=par[0]
intercept=par[1]
rvalue=par[2]
pvalue=par[3]
stderr=par[4]
#------------------------------------
clevsT=np.arange(-1.5,1.6,0.1)
color='seismic'
exp='obs'
decomp='original'
lat,lon=climb.latlon(trend)
#figtrend=myplot.oneplot_trend(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevsT, units=units, subdomain=subdomain, par=par, xlon=par[0][lon], ylat=par[0][lat], season=season, iyr=iyr, fyr=fyr,extent=False)

latS,latN,lonW,lonE,latlim,lonlim=dom.coord_domain(subdomain)
fig,axs= plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},figsize=(8,10))
lons, lats = np.meshgrid(par[0][lon],par[0][lat])
CS1=axs.contourf(lons,lats, 10*par[0][:,:],clevsT,
                transform=ccrs.PlateCarree(),
                cmap=color,extent='both')
levels=[0,0.1,1.0]
cs = axs.contourf(lons,lats,par[3][:,:], transform=ccrs.PlateCarree(),levels=levels,
                hatches=["+", ""], alpha=0.)
    # Draw the coastines for each subplot
axs.coastlines()
axs.add_feature(cfeature.BORDERS, linestyle=':', alpha=1)
    #axs.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='None', facecolor='None'))
axs.set_extent([lonW, lonE, latS,latN])
        # Longitude and latitude labels
axs.set_xticks(np.arange(lonW,lonE,30), crs=ccrs.PlateCarree())
axs.set_yticks(np.arange(latS,latN,10), crs=ccrs.PlateCarree())
    # Adjust the location of the subplots on the page to make room for the colorbar
geom = geometry.box(minx=-79,maxx=-70,miny=-44.7,maxy=-49.3)
axs.add_geometries([geom], facecolor='None',edgecolor='green',crs=cartopy.crs.PlateCarree(), alpha=0.3)
fig.subplots_adjust(bottom=0.35, top=0.7, left=0.20, right=0.80,
                    wspace=0.05, hspace=0.5)
    # Add a colorbar axis at the bottom of the graph
    #([xmin,ymin,dx,dy])
cbar_ax = fig.add_axes([0.2, 0.3, 0.6, 0.02])
    # Draw the colorbar
cbar=fig.colorbar(CS1, cax=cbar_ax,orientation='horizontal',label='%s $(10 \, yr)^{-1}$'%(units))
    ## Add a big title at the top
#plt.suptitle('%s %s %s %s %s %s %s-%s'%(variable,decomp,exp,model,subdomain,season,str(iyr),str(fyr)),y=0.73)

ofileT='spatialtrend_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
fig.savefig(plotsDir+ofileT+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileT)
#plt.show()

'''
ratio=10*abs(trend)/std_det

clevs=np.arange(0,1.5,0.01)
color='BuGn'
exp='obs'
decomp='original'
lat,lon=climb.latlon(trend)
figratio=myplot.oneplot_ds(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevs, color=color, units=units, subdomain=subdomain, ds=ratio, xlon=clim[lon], ylat=clim[lat], season=season, iyr=iyr, fyr=fyr,extent=False)

ofileR='ratio_trend-stddet_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figratio.savefig(plotsDir+ofileR+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileR)
plt.show()
'''
