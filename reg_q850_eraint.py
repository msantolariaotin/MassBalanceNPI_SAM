'''
for i in DJF MAM JJA SON;do
python reg_q850_eraint.py ${i}
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

ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180))
ds = ds.roll(longitude=int(len(ds['longitude']) / 2), roll_coords=True)
#field=dom.field_dom(ds,domain)
units=ds.units
lat,lon=climb.latlon(ds)
ylat=ds.coords['latitude']
xlon=ds.coords['longitude']

field=ds.sel(latitude=slice(-20,-90),longitude=slice(-170,-10))

if units=='kg kg**-1':
    print(units)
    field=field*1000
    units='g/kg'
    print('new unit',units)
else:
    print(units)

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
mode='ELNINO34'
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
color='RdBu_r'
exp='obs'
decomp='original'
lat,lon=climb.latlon(trend)
figtrend=myplot.oneplot_rvalue(variable=variable, decomp=decomp, exp=exp, model=model, clevs=clevsT, units=units, subdomain=subdomain, par=par, xlon=par[0][lon], ylat=par[0][lat], season=season, iyr=iyr, fyr=fyr,extent=mapa)
ofileT='reg_'+mode+'_NOAA_'+variable+'_'+model+'_'+exp+'_'+decomp+'_'+domain+'_'+season+'_'+str(iyr)+'_'+str(fyr)
figtrend.savefig(plotsDir+ofileT+'.png',format='png')
print('Figure save at ',plotsDir, 'as',ofileT)
plt.show()
