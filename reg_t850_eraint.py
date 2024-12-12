'''
for i in DJF MAM JJA SON;do
python reg_t850_eraint.py ${i}
done
'''
import sys
sys.path.insert(1, '/home/msantolaria/Documents/MyPythonLibrary/ClimAnag/')
import climbasis as climb
from climbasis import *
import domain as dom
import myplot
import glob
from myplot import *

# In[3]:


source='/home/msantolaria/Documents/MyResearch/GabiSAM/UBSam/'
resultsDir=source + 'Results/'
plotsDir=source + 'Plots/'


#resultsDir='/media/maria/MARIAEXT2/WinterTrendsObs/'
sourceData='/home/msantolaria/Documents/Data/'


season=sys.argv[1]

if season=='DJF':
    iyr=1979
else:
    iyr=1980
imon=1
fmon=12
fyr=2014
mapa=True

domain='SHWest'
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

field=ds.sel(latitude=slice(-20,-90),longitude=slice(-170,-10))
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

print(anoms.shape)
anoms_detrend=climb.detrend_dim(vals, 'time', deg=1)

clim=vals.mean('time')
std=vals.std('time')
std_det=anoms_detrend.std('time')
#--------------------------------
subdomain=domain
#-----------------------------------------
#Regression
sourceIndex='/home/msantolaria/Documents/MyResearch/GabiSAM/UBSam/Index/'
mode='SAM'
if season=='DJF':
    indexName=mode+'_'+season+'_19791980_20132014.txt'
else:
    indexName=mode+'_'+season+'_1980_2014.txt'
index=np.loadtxt(sourceIndex+indexName)
print(indexName)
print(index.shape)
#-----------------------------------
par=climb.trend_vect(index,anoms,'time')
trend=par[0]
intercept=par[1]
rvalue=par[2]
pvalue=par[3]
stderr=par[4]
#------------------------------------
clevsT=np.arange(-1.0,1.1,0.1)
color='PiYG'
exp='obs'
decomp='original'
lat,lon=climb.latlon(trend)
#figtrend=myplot.oneplot_rvalue(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevsT, units=units, subdomain=subdomain, par=par, xlon=par[0][lon], ylat=par[0][lat], season=season, iyr=iyr, fyr=fyr,extent=mapa)

latS,latN,lonW,lonE,latlim,lonlim=dom.coord_domain(subdomain)
fig,axs= plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},figsize=(8,10))
lons, lats = np.meshgrid(par[0][lon],par[0][lat])
CS1=axs.contourf(lons,lats, par[2][:,:],clevsT,
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
axs.add_geometries([geom], facecolor='None',edgecolor='red',crs=ccrs.PlateCarree(), alpha=0.3)
fig.subplots_adjust(bottom=0.35, top=0.7, left=0.20, right=0.80,
                    wspace=0.05, hspace=0.5)
    # Add a colorbar axis at the bottom of the graph
    #([xmin,ymin,dx,dy])
cbar_ax = fig.add_axes([0.2, 0.3, 0.6, 0.02])
    # Draw the colorbar
cbar=fig.colorbar(CS1, cax=cbar_ax,orientation='horizontal',label='%s $(10 \, yr)^{-1}$'%(units))
    ## Add a big title at the top
#plt.suptitle('%s %s %s %s %s %s %s-%s'%(variable,decomp,exp,model,subdomain,season,str(iyr),str(fyr)),y=0.73)

ofileT='reg_'+mode+'_NOAA_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
fig.savefig(plotsDir+ofileT+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileT)
#plt.show()
