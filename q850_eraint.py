'''
for i in DJF MAM JJA SON;do
python q850_eraint.py ${i}
done
'''
import sys
sys.path.insert(1, '/home/msantolaria/Documents/MyPythonLibrary/ClimAnag/')
import climbasis as climb
from climbasis import *
import domain as dom
import myplot
import glob


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
variable='q850'
model=str('eraint')
fileName ='q850_mon_eraint_1979-2015.nc'
ds = xr.open_dataset(sourceData+model+'/'+fileName)['q']
#filename = 'q925.mon.eraint_197901_201512.nc'
ylat=ds.coords['latitude']
xlon=ds.coords['longitude']
units=ds.units

ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180))
ds = ds.roll(longitude=int(len(ds['longitude']) / 2), roll_coords=True)
#field=dom.field_dom(ds,domain)

lat,lon=climb.latlon(ds)
ylat=ds.coords['latitude']
xlon=ds.coords['longitude']

field=ds.sel(latitude=slice(-20,-90),longitude=slice(-180,0))
latS,latN,lonW,lonE,latlim,lonlim=dom.coord_domain(domain)
#field=dom.field_sel_dom(ds,domain)
if units=='kg kg**-1':
    print(units)
    field=field*1000
    units='g/kg'
    print('new unit',units)
else:
    print(units)

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


##Clim-----------------
#clevs=np.arange(0,14.25,0.25)
clevs=np.arange(0,16,2)
color='Blues'
exp='obs'
decomp='original'
lat,lon=climb.latlon(clim)
figclim=myplot.oneplot_ds(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevs, color=color, units=units, subdomain=subdomain, ds=clim, xlon=clim[lon], ylat=clim[lat], season=season, iyr=iyr, fyr=fyr,extent=mapa)
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
#-----------------------------------
par=climb.trend_vect(anoms.time,anoms,'time')
trend=par[0]
intercept=par[1]
rvalue=par[2]
pvalue=par[3]
stderr=par[4]
#------------------------------------
clevsT=np.arange(-1.5,1.6,0.1)
color='RdBu_r'
exp='obs'
decomp='original'
lat,lon=climb.latlon(trend)
figtrend=myplot.oneplot_trend(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevsT, units=units, subdomain=subdomain, par=par, xlon=par[0][lon], ylat=par[0][lat], season=season, iyr=iyr, fyr=fyr,extent=False)
ofileT='spatialtrend_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figtrend.savefig(plotsDir+ofileT+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileT)
plt.show()
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
