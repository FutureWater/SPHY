# The Spatial Processes in HYdrology (SPHY) model:
# A spatially distributed hydrological model 
# Copyright (C) 2013-2019  FutureWater
# Email: sphy@futurewater.nl
#
# Authors (alphabetical order):
# P. Droogers, J. Eekhout, W. Immerzeel, S. Khanal, A. Lutz, G. Simons, W. Terink
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# This model uses the sphy_config.cfg as configuration file.

import time, shutil, os, glob, configparser, sys, csv, math
import pandas as pd
import pcraster as pcr
import pcraster.framework as pcrm
import numpy as np

tic = time.time()

# Read the model configuration file
config = configparser.RawConfigParser()
config.read(sys.argv[1])

class sphy(pcrm.DynamicModel):
	def __init__(self):
		# Print model info
		print('The Spatial Processes in HYdrology (SPHY) model is')
		print('developed and owned by FutureWater, Wageningen, The Netherlands')
		print('Version 3.0, released June 2019')
		print(' ')

		#-Missing value definition
		self.MV= -9999

		# Read the modules to be used
		self.GlacFLAG = config.getint('MODULES','GlacFLAG')
		self.SnowFLAG = config.getint('MODULES','SnowFLAG')
		self.RoutFLAG = config.getint('MODULES','RoutFLAG')
		self.ResFLAG = config.getint('MODULES','ResFLAG')
		self.LakeFLAG = config.getint('MODULES','LakeFLAG')
		self.DynVegFLAG = config.getint('MODULES','DynVegFLAG')
		self.GroundFLAG = config.getint('MODULES','GroundFLAG')
		self.SedFLAG = config.getint('MODULES','SedFLAG')
		self.SedTransFLAG = config.getint('MODULES','SedTransFLAG')

		# import the required modules
		import datetime, calendar, ET, rootzone, subzone
		import utilities.reporting as reporting
		import utilities.timecalc as timecalc
		import utilities.netcdf2PCraster as netcdf2PCraster
		from math import pi
		#-standard python modules
		self.datetime = datetime
		self.calendar = calendar
		self.pi = pi
		#-FW defined modules
		self.reporting = reporting
		self.timecalc = timecalc
		self.netcdf2PCraster = netcdf2PCraster
		self.ET = ET
		self.rootzone = rootzone
		self.subzone = subzone
		del datetime, calendar, pi, reporting, timecalc, ET, rootzone, subzone
		#-import additional modules if required
		if self.GlacFLAG == 1:
			self.SnowFLAG = 1
			self.GroundFLAG = 1

		#-read the input and output directories from the configuration file
		self.inpath = config.get('DIRS', 'inputdir')
		self.outpath = config.get('DIRS', 'outputdir')

		#-set the timing criteria
		sy1 = config.getint('TIMING', 'startyear_timestep1')
		sm1 = config.getint('TIMING', 'startmonth_timestep1')
		sd1 = config.getint('TIMING', 'startday_timestep1')
		sy = config.getint('TIMING', 'startyear')
		sm = config.getint('TIMING', 'startmonth')
		sd = config.getint('TIMING', 'startday')
		ey = config.getint('TIMING', 'endyear')
		em = config.getint('TIMING', 'endmonth')
		ed = config.getint('TIMING', 'endday')
		self.startdate = self.datetime.datetime(sy,sm,sd)
		self.enddate = self.datetime.datetime(ey,em,ed)
		self.ts1date = self.datetime.datetime(sy1,sm1,sd1)
		self.dateAfterUpdate = self.startdate - self.datetime.timedelta(days=1)  #-only required for glacier retreat (create dummy value here to introduce the variable)

		#-set date input for reporting
		self.startYear = sy
		self.endYear = ey
		self.spinUpYears = config.getint('TIMING', 'spinupyears')
		self.simYears = self.endYear - self.startYear - self.spinUpYears + 1

		#-set the 2000 julian date number
		self.julian_date_2000 = 2451545
		#-read name of reporting table
		self.RepTab = config.get('REPORTING','RepTab')
		#-set the option to calculate the fluxes in mm for the upstream area
		self.mm_rep_FLAG = config.getint('REPORTING','mm_rep_FLAG')

		#-set the option to calculate the fluxes per component in mm for the upstream area
		pars = ['Prec', 'ETa', 'GMelt', 'QSNOW', 'QROOTR', 'QROOTD', 'QRAIN', 'QGLAC', 'QBASE', 'QTOT', 'Seep']
		for i in pars:
			var = i + '_mm_FLAG'
			setattr(self, var, config.getint('REPORTING', var))

		#-set the option to calculate the timeseries of the water balance
		self.wbal_TSS_FLAG = config.getint('REPORTING','wbal_TSS_FLAG')

		#-setting clone map
		self.clonefile = self.inpath + config.get('GENERAL','mask')
		pcr.setclone(self.clonefile)
		self.clone = pcr.ifthen(pcr.readmap(self.clonefile), pcr.boolean(1))
		
		self.cellArea = pcr.cellvalue(pcr.cellarea(),1)[0]

		#-read general maps
		self.DEM = pcr.readmap(self.inpath + config.get('GENERAL','dem'))
		self.Slope = pcr.readmap(self.inpath + config.get('GENERAL','Slope'))
		self.Locations = pcr.readmap(self.inpath + config.get('GENERAL','locations'))

		#-read soil calibration fractions
		self.RootFieldFrac = config.getfloat('SOIL_CAL', 'RootFieldFrac')
		self.RootSatFrac = config.getfloat('SOIL_CAL', 'RootSatFrac')
		self.RootDryFrac = config.getfloat('SOIL_CAL', 'RootDryFrac')
		self.RootWiltFrac = config.getfloat('SOIL_CAL', 'RootWiltFrac')
		self.RootKsatFrac = config.getfloat('SOIL_CAL', 'RootKsatFrac')

		#-read soil maps
		#-check for PedotransferFLAG
		self.PedotransferFLAG = config.getint('PEDOTRANSFER', 'PedotransferFLAG')
		#-if pedotransfer functions are used read the sand, clay, organic matter and bulk density maps, otherwise read the soil hydraulic properties
		if self.PedotransferFLAG == 1:
			import utilities.pedotransfer
			self.pedotransfer = utilities.pedotransfer
			del utilities.pedotransfer

			#-read init processes pedotransfer
			self.pedotransfer.init(self, pcr, config, np)
		else:
			#self.Soil = pcr.readmap(self.inpath + config.get('SOIL','Soil'))
			self.RootFieldMap = pcr.readmap(self.inpath + config.get('SOIL','RootFieldMap')) * self.RootFieldFrac
			self.RootSatMap = pcr.readmap(self.inpath + config.get('SOIL','RootSatMap')) * self.RootSatFrac
			self.RootDryMap = pcr.readmap(self.inpath + config.get('SOIL','RootDryMap')) * self.RootDryFrac
			self.RootWiltMap = pcr.readmap(self.inpath + config.get('SOIL','RootWiltMap')) * self.RootWiltFrac
			self.RootKsat = pcr.readmap(self.inpath + config.get('SOIL','RootKsat')) * self.RootKsatFrac
			self.SubSatMap = pcr.readmap(self.inpath + config.get('SOIL','SubSatMap'))
			self.SubFieldMap = pcr.readmap(self.inpath + config.get('SOIL','SubFieldMap'))
			self.SubKsat = pcr.readmap(self.inpath + config.get('SOIL','SubKsat'))
			self.RootDrainVel = self.RootKsat * self.Slope

		#-Read and set the soil parameters
		pars = ['CapRiseMax','RootDepthFlat','SubDepthFlat']
		for i in pars:
			try:
				setattr(self, i, pcr.readmap(self.inpath + config.get('SOILPARS',i)))
			except:
				setattr(self, i, config.getfloat('SOILPARS',i))
		
		# groundwater storage as third storage layer. This is used instead of a fixed bottomflux
		if self.GroundFLAG == 1:
			import modules.groundwater 
			self.groundwater = modules.groundwater
			del modules.groundwater

			#-read init processes groundwater
			self.groundwater.init(self, pcr, config)

		else:
			# if groundwater module is not used, read seepage and gwl_base
			self.SeepStatFLAG = config.getint('SOILPARS','SeepStatic')
			if self.SeepStatFLAG == 0: # set the seepage map series
				self.Seepmaps = self.inpath + config.get('SOILPARS', 'SeePage')
			else: #-set a static map or value for seepage
				try:
					self.SeePage = pcr.readmap(self.inpath + config.get('SOILPARS','SeePage'))
				except:
					self.SeePage = config.getfloat('SOILPARS','SeePage')
			try:
				self.GWL_base = pcr.readmap(self.inpath + config.get('SOILPARS','GWL_base'))
			except:
				self.GWL_base = config.getfloat('SOILPARS','GWL_base')

			self.SubDrainVel = self.SubKsat * self.Slope

		#-calculate soil properties
		self.RootField = self.RootFieldMap * self.RootDepthFlat
		self.RootSat = self.RootSatMap * self.RootDepthFlat
		self.RootDry = self.RootDryMap * self.RootDepthFlat
		self.RootWilt = self.RootWiltMap * self.RootDepthFlat
		self.SubSat = self.SubSatMap * self.SubDepthFlat
		self.SubField = self.SubFieldMap * self.SubDepthFlat
		self.RootTT = pcr.max((self.RootSat - self.RootField) / self.RootKsat, 0.0001)
		self.SubTT = pcr.max((self.SubSat - self.SubField) / self.SubKsat, 0.0001)
		# soil max and soil min for scaling of gwl if groundwater module is not used
		if self.GroundFLAG == 0:
			self.SoilMax = self.RootSat + self.SubSat
			self.SoilMin = self.RootDry + self.SubField

		#-read land use map
		self.LandUse = pcr.readmap(self.inpath + config.get('LANDUSE','LandUse'))

		#-Use the dynamic vegetation module
		if self.DynVegFLAG == 1:
			#-import dynamic vegetation module
			import modules.dynamic_veg
			self.dynamic_veg = modules.dynamic_veg
			del modules.dynamic_veg

			#-read init processes dynamic vegetation
			self.dynamic_veg.init(self, pcr, config)
		#-read the crop coefficient table if the dynamic vegetation module is not used
		else:
			self.KcStatFLAG = config.getint('LANDUSE', 'KCstatic')
			if self.KcStatFLAG == 1:
				#-read land use map and kc table
				self.kc_table = self.inpath + config.get('LANDUSE','CropFac')
				self.Kc = pcr.lookupscalar(self.kc_table, self.LandUse)
			else:
				#-set the kc map series
				self.Kcmaps = self.inpath + config.get('LANDUSE', 'KC')

		#-read the p factor table if the plant water stress module is used
		self.PlantWaterStressFLAG = config.getint('PWS','PWS_FLAG')
		if self.PlantWaterStressFLAG == 1:
			PFactor = self.inpath + config.get('PWS', 'PFactor')
			self.PMap = pcr.lookupscalar(PFactor, self.LandUse)

		#-read and set glacier maps and parameters if glacier module is used
		if self.GlacFLAG:
			#-import glacier module
			import modules.glacier
			self.glacier = modules.glacier
			del modules.glacier

			#-read init processes glacier module
			self.glacier.init(self, pcr, config, pd, np, os)

		#-read and set snow maps and parameters if snow modules are used
		if self.SnowFLAG == 1:
			#-import snow module
			import modules.snow
			self.snow = modules.snow
			del modules.snow
		
			#-read init processes glacier module
			self.snow.init(self, pcr, config)

		#-read and set climate forcing and the calculation of etref

		#-read precipitation data
		#-read flag for precipitation forcing by netcdf
		self.precNetcdfFLAG = config.getint('CLIMATE', 'precNetcdfFLAG')
		if self.precNetcdfFLAG == 1:
			#-read configuration for forcing by netcdf
			self.netcdf2PCraster.getConfigNetcdf(self, config, 'Prec', 'CLIMATE')

			#-determine x,y-coordinates of netcdf file and model domain and indices of netcdf corresponding to model domain
			self.netcdf2PCraster.netcdf2pcrInit(self, pcr, 'Prec')
		else:
			#-read precipitation forcing folder
			self.Prec = self.inpath + config.get('CLIMATE','Prec')


		#-read precipitation data
		#-read flag for temperature forcing by netcdf
		self.tempNetcdfFLAG = config.getint('CLIMATE', 'tempNetcdfFLAG')
		if self.tempNetcdfFLAG == 1:
			#-read configuration for forcing by netcdf
			self.netcdf2PCraster.getConfigNetcdf(self, config, 'Temp', 'CLIMATE')

			#-determine x,y-coordinates of netcdf file and model domain and indices of netcdf corresponding to model domain
			self.netcdf2PCraster.netcdf2pcrInit(self, pcr, 'Temp')
		else:
			#-read temperature forcing folder
			self.Tair = self.inpath + config.get('CLIMATE','Tair')
		#-read flag for etref time series input
		self.ETREF_FLAG = config.getint('ETREF','ETREF_FLAG')
		#-determine the use of a given etref time-series or calculate etref using Hargreaves
		if self.ETREF_FLAG == 1:
			self.ETref = self.inpath + config.get('ETREF','ETref')
		else:
			self.Lat = pcr.readmap(self.inpath + config.get('ETREF','Lat'))
			#-read flag for minimum temperature forcing by netcdf
			self.TminNetcdfFLAG = config.getint('ETREF', 'TminNetcdfFLAG')
			if self.TminNetcdfFLAG == 1:
				#-read configuration for forcing by netcdf
				self.netcdf2PCraster.getConfigNetcdf(self, config, 'Tmin', 'ETREF')

				#-determine x,y-coordinates of netcdf file and model domain and indices of netcdf corresponding to model domain
				self.netcdf2PCraster.netcdf2pcrInit(self, pcr, 'Tmin')
			else:
				self.Tmin = self.inpath + config.get('ETREF','Tmin')
			#-read flag for maximum temperature forcing by netcdf
			self.TmaxNetcdfFLAG = config.getint('ETREF', 'TmaxNetcdfFLAG')
			if self.TmaxNetcdfFLAG == 1:
				#-read configuration for forcing by netcdf
				self.netcdf2PCraster.getConfigNetcdf(self, config, 'Tmax', 'ETREF')

				#-determine x,y-coordinates of netcdf file and model domain and indices of netcdf corresponding to model domain
				self.netcdf2PCraster.netcdf2pcrInit(self, pcr, 'Tmax')
			else:
				self.Tmax = self.inpath + config.get('ETREF','Tmax')
			self.Gsc = config.getfloat('ETREF', 'Gsc')
			import hargreaves
			self.Hargreaves = hargreaves
			del hargreaves
			
		#-read and set routing maps and parameters
		if self.RoutFLAG == 1:
			import modules.routing 
			self.routing = modules.routing
			del modules.routing
			
			#-read init processes routing
			self.routing.init(self, pcr, config)

		#-read and set routing maps and parameters
		if self.ResFLAG == 1 or self.LakeFLAG == 1:
			#-import advanced routing module
			import modules.advanced_routing
			self.advanced_routing = modules.advanced_routing
			del modules.advanced_routing
			
			#-read init processes advanced routing
			self.advanced_routing.init(self, pcr, config)

		#-read lake maps and parameters if lake module is used
		if self.LakeFLAG == 1:
			#-import lakes module
			import modules.lakes
			self.lakes = modules.lakes
			del modules.lakes

			#-read init processes lakes
			self.lakes.init(self, pcr, config)

		#-read reservior maps and parameters if reservoir module is used
		if self.ResFLAG == 1:
			#-import reservoirs module
			import modules.reservoirs
			self.reservoirs = modules.reservoirs
			del modules.reservoirs

			#-read init processes reservoirs
			self.reservoirs.init(self, pcr, config)

		#-read flag for calculation of ET in reservoirs
		self.ETOpenWaterFLAG = config.getint('OPENWATER', 'ETOpenWaterFLAG')
		if self.ETOpenWaterFLAG == 1:
			#-read kc value for open water
			self.kcOpenWater = config.getfloat('OPENWATER', 'kcOpenWater')
			#-read openwater fraction map
			self.openWaterFrac = pcr.readmap(self.inpath + config.get('OPENWATER', 'openWaterFrac'))
			#-determine openwater map with values of each reservoir/lake in the extent of the openwater
			self.openWater = pcr.ifthenelse(self.openWaterFrac > 0, pcr.scalar(1), pcr.scalar(0))
			self.openWaterNominal = pcr.clump(pcr.nominal(self.openWater))
			self.openWaterNominal = pcr.nominal(pcr.areamaximum(pcr.scalar(self.ResID), self.openWaterNominal))
		else:
			#-set all cells to 0 for openwater fraction map
			self.openWaterFrac = self.DEM * 0
			self.openWater = 0
			self.ETOpenWater = 0

		#-read maps and parameters for infiltration excess
		self.InfilFLAG = config.getfloat('INFILTRATION', 'Infil_excess')
		if self.InfilFLAG == 1:
			self.K_eff = config.getfloat('INFILTRATION', 'K_eff')
			try:
				self.Alpha = config.getfloat('INFILTRATION', 'Alpha')
			except:
				self.Alpha = pcr.readmap(self.inpath + config.get('INFILTRATION', 'Alpha'))
			try:
				self.Labda_Infil = config.getfloat('INFILTRATION', 'Labda_infil')
			except:
				self.Labda_Infil = pcr.readmap(self.inpath + config.get('INFILTRATION', 'Labda_infil'))
			try:
				self.paved_table = self.inpath + config.get('INFILTRATION','PavedFrac')
				self.pavedFrac = pcr.lookupscalar(self.paved_table, self.LandUse)
			except:
				self.pavedFrac = 0

		#-read maps and parameters for soil erosion
		if self.SedFLAG == 1:
			#-read soil erosion model selector (1 for MUSLE, 2 for MMF)
			self.SedModel = config.getfloat('SEDIMENT', 'SedModel')

			#-read rock fraction map
			self.RockFrac = pcr.readmap(self.inpath + config.get('SEDIMENT', 'RockFrac'))

			#-Read flag if channels should be excluded from the detachment by runoff calculation
			self.exclChannelsFLAG = config.getint('SEDIMENT', 'exclChannelsFLAG')
			
			#-determine hillslope map if channels should be excluded
			if self.exclChannelsFLAG == 1:
				#-determine upstream area map
				self.UpstreamArea = pcr.accuflux(self.FlowDir, 1) * pcr.cellarea() / 10**6

				#-determine upstream area larger than upstream_km2 and define hillslope cells based on upstream area
				self.Upstream_km2 = config.getfloat('SEDIMENT', 'upstream_km2')
				self.Hillslope = pcr.scalar(self.UpstreamArea <= self.Upstream_km2)

			#-read MUSLE input parameters
			if self.SedModel == 1:
				#-import musle module
				import modules.musle
				self.musle = modules.musle
				del modules.musle

				#-read init processes musle
				self.musle.init(self, pcr, config)

			#-read MMF input parameters
			if self.SedModel == 2:
				#-import mmf module
				import modules.mmf
				self.mmf = modules.mmf
				del modules.mmf

				#-read init processes mmf
				self.mmf.init(self, pcr, config)

			#-read INCA input parameters
			if self.SedModel == 3:
				#-import INCA module
				import modules.inca
				self.inca = modules.inca
				del modules.inca

				#-read init processes INCA
				self.inca.init(self, pcr, config)

			#-read SHETRAN input parameters
			if self.SedModel == 4:
				#-import SHETRAN module
				import modules.shetran
				self.shetran = modules.shetran
				del modules.shetran

				#-read init processes SHETRAN
				self.shetran.init(self, pcr, config)

			#-read DHSVM input parameters
			if self.SedModel == 5:
				#-import DHSVM module
				import modules.dhsvm
				self.dhsvm = modules.dhsvm
				del modules.dhsvm

				#-read init processes DHSVM
				self.dhsvm.init(self, pcr, config)

			#-read HSPF input parameters
			if self.SedModel == 6:
				#-import HSPF module
				import modules.hspf
				self.hspf = modules.hspf
				del modules.hspf

				#-read init processes HSPF
				self.hspf.init(self, pcr, config)

			#-read input parameters for sediment transport
			if self.SedTransFLAG == 1:
				#-import sediment transport module
				import modules.sediment_transport
				self.sediment_transport = modules.sediment_transport
				del modules.sediment_transport

				#-read init processes sediment transport
				self.sediment_transport.init(self, pcr, config, csv, np)

		#-set the global option for radians
		pcr.setglobaloption('radians')

	#-initial section
	def initial(self):
		#-timer
		self.counter = (self.startdate - self.ts1date).days
		#-initial date
		self.curdate = self.startdate

		#-initial soil properties
		#-initial rootwater content
		if not config.get('SOIL_INIT','RootWater'):
			self.RootWater = self.RootField
		else:
			try:
				self.RootWater = config.getfloat('SOIL_INIT','RootWater')
			except:
				self.RootWater = pcr.readmap(self.inpath + config.get('SOIL_INIT','RootWater'))
		#-initial water content in subsoil
		if not config.get('SOIL_INIT','SubWater'):
			self.SubWater = self.SubField
		else:
			try:
				self.SubWater = config.getfloat('SOIL_INIT','SubWater')
			except:
				self.SubWater = pcr.readmap(self.inpath + config.get('SOIL_INIT','SubWater'))
		#-initial water storage in rootzone + subsoil
		self.SoilWater = self.RootWater + self.SubWater
		#-initial capillary rise
		try:
			self.CapRise = config.getfloat('SOIL_INIT','CapRise')
		except:
			self.CapRise = pcr.readmap(self.inpath + config.get('SOIL_INIT','CapRise'))
		#-initial drainage from rootzone
		try:
			self.RootDrain = config.getfloat('SOIL_INIT','RootDrain')
		except:
			self.RootDrain = pcr.readmap(self.inpath + config.get('SOIL_INIT','RootDrain'))

		if self.DynVegFLAG == 1:
			#-read initial conditions dynamic vegetation
			self.dynamic_veg.initial(self, pcr)

		elif self.KcStatFLAG == 0:
			#-set initial kc value to one, if kc map is not available for first timestep
			self.KcOld = pcr.scalar(1)

		#-initial groundwater properties
		if self.GroundFLAG == 1:
			#-read initial conditions groundwater module
			self.groundwater.initial(self, pcr, config)
		else:
			#-initial drainage from subsoil
			try:
				self.SubDrain = config.getfloat('SOIL_INIT','SubDrain')
			except:
				self.SubDrain = pcr.readmap(self.inpath + config.get('SOIL_INIT','SubDrain'))
			#-initial seepage value if seepage map series is used
			if self.SeepStatFLAG == 0:
				self.SeepOld = pcr.scalar(0)

		#-initial snow properties
		if self.SnowFLAG:
			#-read initial conditions snow module
			self.snow.initial(self, pcr, config)
		else:
			self.SnowStore = pcr.scalar(0)

		#-initial glacier properties
		if self.GlacFLAG:
			#-read initial conditions glacier
			self.glacier.initial(self, pcr, pd)
		else:
			self.GlacFrac = pcr.scalar(0)

		#-initial routed total runoff and of individual components
		if self.RoutFLAG == 1:
			#-read init processes routing
			self.routing.initial(self, pcr, config)

		#-initial storage in lakes and reservoirs
		if self.LakeFLAG == 1 or self.ResFLAG == 1:
			#-Read initial storages from table/reservoir file
			if self.LakeFLAG == 1:
				#-read initial conditions lakes
				self.lakes.initial(self, pcr, config)
			if self.ResFLAG == 1:
				#-read initial conditions reservoirs
				self.reservoirs.initial(self, pcr, config)

			#-read init processes advanced routing
			self.advanced_routing.initial(self, pcr, config)

		#-Initial routed volume of sediment
		if self.SedFLAG == 1 and self.SedTransFLAG == 1:
			#-read init processes sediment transport
			self.sediment_transport.initial(self, pcr, config)

		#-Initial values for reporting and setting of time-series
		#-set time-series reporting for mm flux from upstream area for prec and eta
		pars = ['Prec','ETa','GMelt','QSNOW','QROOTR','QROOTD','QRAIN','QGLAC','QBASE','QTOT','Seep']
		for i in pars:
			if self.mm_rep_FLAG == 1 and eval('self.' + i + '_mm_FLAG'):
				setattr(self, i + 'SubBasinTSS', pcrm.TimeoutputTimeseries(i + 'SubBasinTSS', self, self.Locations, noHeader=False))

		#-WATER BALANCE
		self.oldRootWater = self.RootWater
		self.oldSubWater = self.SubWater
		if self.GroundFLAG:
			self.oldGw = self.Gw
		if self.wbal_TSS_FLAG:
			self.wbalTSS = pcrm.TimeoutputTimeseries("wbalTSS", self, self.Locations, noHeader=True)
			self.wbalTotTSS = pcrm.TimeoutputTimeseries("wbalTotTSS", self, self.Locations, noHeader=True)
		self.waterbalanceTot = pcr.scalar(0.)

		#-read reporting options from csv file
		self.reporting.initial(self, pcr, csv, pcrm)
		
		#-set reporting of water balances for lakes
		if self.LakeFLAG == 1 and config.getint('REPORTING', 'Lake_wbal') ==1:
			#-read initial conditions reporting lakes
			self.lakes.initial_reporting(self, pcr, pcrm)
		#-set reporting of water balances for reservoirs
		if self.ResFLAG == 1 and config.getint('REPORTING', 'Res_wbal') == 1:
			#-read initial conditions reporting reservoirs
			self.reservoirs.initial_reporting(self, pcr, pcrm)

	def dynamic(self):
		self.counter+=1
		print(str(self.curdate.day)+'-'+str(self.curdate.month)+'-'+str(self.curdate.year)+'  t = '+str(self.counter))
		
		#-Snow and rain fraction settings for non-glacier part of model cell
		SnowFrac = pcr.ifthenelse(self.SnowStore > 0, pcr.scalar(1 - self.GlacFrac), 0)
		RainFrac = pcr.ifthenelse(self.SnowStore == 0, pcr.scalar(1 - self.GlacFrac), 0)

		#-Read the precipitation time-series
		if self.precNetcdfFLAG == 1:
			#-read forcing by netcdf input
			Precip = self.netcdf2PCraster.netcdf2pcrDynamic(self, pcr, 'Prec')
		else:
			#-read forcing by map input
			Precip = pcr.readmap(pcrm.generateNameT(self.Prec, self.counter))
		PrecipTot = Precip
		#-Report Precip
		self.reporting.reporting(self, pcr, 'TotPrec', Precip)
		self.reporting.reporting(self, pcr, 'TotPrecF', Precip * (1-self.GlacFrac))

		#-Temperature and determine reference evapotranspiration
		if self.tempNetcdfFLAG == 1:
			#-read forcing by netcdf input
			Temp = self.netcdf2PCraster.netcdf2pcrDynamic(self, pcr, 'Temp')
		else:
			#-read forcing by map input
			Temp = pcr.readmap(pcrm.generateNameT(self.Tair, self.counter))

		if self.ETREF_FLAG == 0:
			if self.TminNetcdfFLAG == 1:
				#-read forcing by netcdf input
				TempMin = self.netcdf2PCraster.netcdf2pcrDynamic(self, pcr, 'Tmin')
			else:
				#-read forcing by map input
				TempMin = pcr.readmap(pcrm.generateNameT(self.Tmin, self.counter))
				# TempMin = pcr.readmap(pcrm.generateNameT(self.Tmin, self.curdate.timetuple().tm_yday))
			if self.TmaxNetcdfFLAG == 1:
				#-read forcing by netcdf input
				TempMax = self.netcdf2PCraster.netcdf2pcrDynamic(self, pcr, 'Tmax')
			else:
				#-read forcing by map input
				TempMax = pcr.readmap(pcrm.generateNameT(self.Tmax, self.counter))
				# TempMax = pcr.readmap(pcrm.generateNameT(self.Tmax, self.curdate.timetuple().tm_yday))
			ETref = self.Hargreaves.Hargreaves(pcr, self.Hargreaves.extrarad(self, pcr), Temp, TempMax, TempMin)
		else:
			ETref = pcr.readmap(pcrm.generateNameT(self.ETref, self.counter))
		self.reporting.reporting(self, pcr, 'TotETref', ETref)
		self.reporting.reporting(self, pcr, 'TotETrefF', ETref * (1-self.GlacFrac))

		#-Interception and effective precipitation
		if self.DynVegFLAG == 1:
			#-read dynamic processes dynamic vegetation
			Precip = self.dynamic_veg.dynamic(self, pcr, pcrm, np, Precip, ETref)

		elif self.KcStatFLAG == 0:
			#-Try to read the KC map series
			try:
				self.Kc = pcr.readmap(pcrm.generateNameT(self.Kcmaps, self.counter))
				self.KcOld = self.Kc
			except:
				self.Kc = self.KcOld
		#-report mm effective precipitation for sub-basin averages
		if self.mm_rep_FLAG == 1 and self.Prec_mm_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
			self.PrecSubBasinTSS.sample(pcr.catchmenttotal(Precip * (1-self.GlacFrac), self.FlowDir) / pcr.catchmenttotal(1, self.FlowDir))

		#-Snow, rain, and glacier calculations for glacier fraction of cell
		if self.GlacFLAG:
			#-read dynamic processes glacier
			Rain_GLAC, Snow_GLAC, ActSnowMelt_GLAC, SnowR_GLAC, GlacMelt, GlacPerc, self.GlacR = self.glacier.dynamic(self, pcr, pd, Temp, Precip)
		#-If glacier module is not used, then
		else:
			Rain_GLAC = 0
			Snow_GLAC = 0
			ActSnowMelt_GLAC = 0
			self.TotalSnowStore_GLAC = 0
			SnowR_GLAC = 0
			self.GlacR = 0
			GlacMelt = 0
			GlacPerc = 0

		# Calculate snow and rain for non-glacier part of cell
		if self.SnowFLAG == 1:
			#-read dynamic processes snow
			Rain, self.SnowR, OldTotalSnowStore = self.snow.dynamic(self, pcr, Temp, Precip, Snow_GLAC, ActSnowMelt_GLAC, SnowFrac, RainFrac, SnowR_GLAC)
		else:
			Rain = Precip
			self.SnowR = 0
			OldTotalSnowStore = 0
			self.TotalSnowStore = 0
		#-Report Rain
		self.reporting.reporting(self, pcr, 'TotRain', Rain)
		self.reporting.reporting(self, pcr, 'TotRainF', Rain * (1-self.GlacFrac) + Rain_GLAC)  # for entire cell

		#-Potential evapotranspiration
		ETpot = self.ET.ETpot(ETref, self.Kc)
		if self.ETOpenWaterFLAG == 1:
			self.ETOpenWater = self.ET.ETpot(ETref, self.kcOpenWater)
		#-Report ETpot
		self.reporting.reporting(self, pcr, 'TotETpot', ETpot)
		self.reporting.reporting(self, pcr, 'TotETpotF', ETpot * RainFrac)

		#-Rootzone calculations
		self.RootWater = self.RootWater + self.CapRise
		#-Calculate rootzone runoff
		tempvar = self.rootzone.RootRunoff(self, pcr, RainFrac, Rain)
		#-Rootzone runoff
		RootRunoff = tempvar[0]
		#-Infiltration
		Infil = tempvar[1]
		#-Report infiltration
		self.reporting.reporting(self, pcr, 'Infil', Infil)
	    #-Updated rootwater content
		self.RootWater = pcr.ifthenelse(RainFrac > 0, self.RootWater + Infil, self.RootWater)

		#-Actual evapotranspiration
		if self.PlantWaterStressFLAG == 1:
			etreddry = self.ET.ks(self, pcr, ETpot)
		else:
			etreddry = pcr.max(pcr.min((self.RootWater - self.RootDry) / (self.RootWilt - self.RootDry), 1), 0)
		self.reporting.reporting(self, pcr, 'PlantStress', 1 - etreddry)
		ETact = self.ET.ETact(pcr, ETpot, self.RootWater, self.RootSat, etreddry, RainFrac)
		#-Report the actual evapotranspiration
		self.reporting.reporting(self, pcr, 'TotETact', ETact * (1-self.openWaterFrac) + self.ETOpenWater * self.openWaterFrac)
		#-Actual evapotranspiration, corrected for rain fraction
		ActETact = ETact * RainFrac
		#-Report the actual evapotranspiration, corrected for rain fraction
		self.reporting.reporting(self, pcr, 'TotETactF', ActETact)
		if self.mm_rep_FLAG == 1 and self.ETa_mm_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
			self.ETaSubBasinTSS.sample(pcr.catchmenttotal(ActETact, self.FlowDir) / pcr.catchmenttotal(1, self.FlowDir))
		#-Update rootwater content
		self.RootWater = pcr.max(self.RootWater - ETact, 0)

		#-Calculate drainage
		temp_RootDrain = self.rootzone.RootDrainage(pcr, self.RootWater, self.RootDrain, self.RootField, self.RootSat, self.RootDrainVel, self.RootTT)
		#-Calculate percolation
		temp_rootperc = self.rootzone.RootPercolation(pcr, self.RootWater, self.SubWater, self.RootField, self.RootTT, self.SubSat)
		#-Total sum of water able to leave the soil
		RootOut = temp_RootDrain + temp_rootperc
		#-Calculate new values for drainage and percolation (to be used when RootOut > RootExcess)
		newdrain, newperc = self.rootzone.CalcFrac(pcr, self.RootWater, self.RootField, temp_RootDrain, temp_rootperc)
		#-Determine whether the new values need to be used
		rootexcess = pcr.max(self.RootWater - self.RootField, 0)
		self.RootDrain = pcr.ifthenelse(RootOut > rootexcess, newdrain, temp_RootDrain)
		rootperc = pcr.ifthenelse(RootOut > rootexcess, newperc, temp_rootperc)
		#-Update the RootWater content
		# Roottemp = self.RootWater
		self.RootWater = self.RootWater - (self.RootDrain + rootperc)

		#-Report rootzone percolation, corrected for fraction
		self.reporting.reporting(self, pcr, 'TotRootPF', rootperc * (1 - self.GlacFrac))
		#-Report rootwater content
		self.reporting.reporting(self, pcr, 'StorRootW', self.RootWater * (1-self.openWaterFrac))

		#-Sub soil calculations
		self.SubWater = self.SubWater + rootperc
		if self.GroundFLAG == 0:
			if self.SeepStatFLAG == 0:
				try:
					self.SeePage = pcr.readmap(pcrm.generateNameT(self.Seepmaps, self.counter))
					self.SeepOld = self.SeePage
				except:
					self.SeePage = self.SeepOld

			#-Report seepage
			self.reporting.reporting(self, pcr, 'TotSeepF', pcr.scalar(self.SeePage))
			self.SubWater = pcr.min(pcr.max(self.SubWater - self.SeePage, 0), self.SubSat)
			if self.mm_rep_FLAG == 1 and self.Seep_mm_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
				self.SeepSubBasinTSS.sample(pcr.catchmenttotal(self.SeePage, self.FlowDir) / pcr.catchmenttotal(1, self.FlowDir))
		#-Capillary rise
		self.CapRise = self.subzone.CapilRise(pcr, self.SubField, self.SubWater, self.CapRiseMax, self.RootWater, self.RootSat, self.RootField)
		#-Report capillary rise, corrected for fraction
		self.reporting.reporting(self, pcr, 'TotCapRF', self.CapRise * (1-self.GlacFrac))
		#-Update sub soil water content
		self.SubWater = self.SubWater - self.CapRise
		if self.GroundFLAG == 1:   # sub percolation will be calculated instead of subdrainage
			subperc = self.subzone.SubPercolation(pcr, self.SubWater, self.SubField, self.SubTT, self.Gw, self.GwSat)
			ActSubPerc = subperc * (1-self.GlacFrac)
			#-Report the subzone percolation, corrected for the fraction
			self.reporting.reporting(self, pcr, 'TotSubPF', ActSubPerc)
			#-Update sub soil water content
			self.SubWater = self.SubWater - subperc
		else: # sub drainage will be calculated instead of sub percolation
			self.SubDrain = self.subzone.SubDrainage(pcr, self.SubWater, self.SubField, self.SubSat, self.SubDrainVel, self.SubDrain, self.SubTT)
			#-Report drainage from subzone
			self.reporting.reporting(self, pcr, 'TotSubDF', self.SubDrain)
			#-Update sub soil water content
			self.SubWater = self.SubWater - self.SubDrain
		#-Report rootwater content
		self.reporting.reporting(self, pcr, 'StorSubW', self.SubWater * (1-self.openWaterFrac))

		#-Changes in soil water storage
		OldSoilWater = self.SoilWater
		self.SoilWater = (self.RootWater + self.SubWater) * (1-self.GlacFrac)

		#-Rootzone runoff
		self.RootRR = RootRunoff * RainFrac * (1-self.openWaterFrac)
		#-Report rootzone runoff, corrected for fraction
		self.reporting.reporting(self, pcr, 'TotRootRF', self.RootRR)
		#-Rootzone drainage
		self.RootDR = self.RootDrain * (1-self.GlacFrac) * (1-self.openWaterFrac)
		#-Report rootzone drainage, corrected for fraction
		self.reporting.reporting(self, pcr, 'TotRootDF', self.RootDR)
		#-Rain runoff
		self.RainR = self.RootRR + self.RootDR
		#-Report rain runoff
		self.reporting.reporting(self, pcr, 'TotRainRF', self.RainR)

		#-Groundwater calculations
		if self.GroundFLAG == 1:
			#-read dynamic processes groundwater
			self.groundwater.dynamic(self, pcr, ActSubPerc, GlacPerc)
		else:
			#-Use drainage from subsoil as baseflow
			self.BaseR = self.SubDrain
			#-Groundwater level as scaled between min and max measured gwl
			SoilAct = self.RootWater + self.SubWater
			SoilRel = (SoilAct - self.SoilMin) / (self.SoilMax - self.SoilMin) # scale between 0 (dry) and 1 (wet)
			GWL = self.GWL_base - (SoilRel-0.5) * self.GWL_base
			#-Report groundwater
			self.reporting.reporting(self, pcr, 'GWL', GWL)

		#-Report Total runoff
		TotR = self.BaseR + self.RainR + self.SnowR + self.GlacR
		self.reporting.reporting(self, pcr, 'TotRF', TotR)

		#-Routing for lake and/or reservoir modules
		if self.LakeFLAG == 1 or self.ResFLAG == 1:
			#-read dynamic processes advanced routing
			Q = self.advanced_routing.dynamic(self, pcr, pcrm, config, TotR, self.ETOpenWater, PrecipTot)

		#-Normal routing module
		elif self.RoutFLAG == 1:
			Q = self.routing.dynamic(self, pcr, TotR)
			
			if self.GlacFLAG:
				#-read dynamic reporting processes glacier
				self.glacier.dynamic_reporting(self, pcr, pd, np)

		#-Water balance
		if self.GlacFLAG and self.GlacRetreat == 1:
			GlacTable_MODid = self.GlacTable.loc[:,['FRAC_GLAC', 'ICE_DEPTH']]
			GlacTable_MODid['ICE_DEPTH'] = GlacTable_MODid['ICE_DEPTH'] * GlacTable_MODid['FRAC_GLAC']
			GlacTable_MODid = GlacTable_MODid.groupby(GlacTable_MODid.index).sum()
			GlacTable_MODid.fillna(0., inplace=True)
			#-Report pcraster map of glacier depth
			iceDepth = np.zeros(self.ModelID_1d.shape)
			iceDepth[self.GlacierKeys] = GlacTable_MODid['ICE_DEPTH']
			iceDepth = iceDepth.reshape(self.ModelID.shape)
			iceDepth = pcr.numpy2pcr(pcr.Scalar, iceDepth, self.MV)
			iceDepth = pcr.ifthen(self.clone, iceDepth)  #-only use values where clone is True
			iceDepth = iceDepth * 1000 # in mm
			#-change in storage
			dS = ((self.RootWater - self.oldRootWater) + (self.SubWater - self.oldSubWater)) * (1-self.GlacFrac) + (self.Gw - self.oldGw) + \
				(self.TotalSnowStore-OldTotalSnowStore) + (iceDepth - self.oldIceDepth)
			#-set old state variables for glacier
			self.oldIceDepth = iceDepth; iceDepth = None; del iceDepth;
			GlacTable_MODid = None; del GlacTable_MODid;
		elif self.GroundFLAG:
			#-change in storage
			dS = ((self.RootWater - self.oldRootWater) + (self.SubWater - self.oldSubWater)) * (1-self.GlacFrac) + (self.Gw - self.oldGw) + \
				(self.TotalSnowStore-OldTotalSnowStore)
			# set old state variables for groundwater
			self.oldGw = self.Gw
		else:
			#-change in storage
			dS = ((self.RootWater - self.oldRootWater) + (self.SubWater - self.oldSubWater)) * (1-self.GlacFrac) + (self.TotalSnowStore-OldTotalSnowStore)

		#-water balance per time step
		if self.GroundFLAG:
			waterbalance = Precip - ActETact - self.BaseR - self.RainR - self.SnowR - self.GlacR - dS
		else:
			waterbalance = Precip - ActETact - self.BaseR - self.RainR - self.SnowR - dS - self.SeePage
		self.reporting.reporting(self, pcr, 'wbal', waterbalance)

		#-total water balance
		self.waterbalanceTot = self.waterbalanceTot + waterbalance
		#-report water balance and accumulated water balance
		if self.wbal_TSS_FLAG and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
			self.wbalTSS.sample(pcr.catchmenttotal(waterbalance, self.FlowDir) / pcr.catchmenttotal(1., self.FlowDir))
			self.wbalTotTSS.sample(pcr.catchmenttotal(self.waterbalanceTot, self.FlowDir) / pcr.catchmenttotal(1., self.FlowDir))
		# set old state variables
		self.oldRootWater = self.RootWater
		self.oldSubWater = self.SubWater
		waterbalance = None; del waterbalance; dS = None; del dS;
		#-End of water balance calculations

		#-Sediment yield
		if self.SedFLAG == 1:
			#-determine runoff in mm per day
			if self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1:
				Runoff = (Q * 3600 * 24) / pcr.cellarea() * 1000
			else: 
				Runoff = TotR

			#-MUSLE
			if self.SedModel == 1:
				#-read dynamic processes musle
				self.musle.dynamic(self, pcr, Runoff)

				#-sediment transport
				if self.SedTransFLAG == 1:
					#-read dynamic sediment transport processes musle
					self.sediment_transport.dynamic_musle(self, pcr)

			#-Modified Morgan-Morgan-Finney model
			if self.SedModel == 2:
				#-determine soil erosion in transport (G)
				G = self.mmf.dynamic(self, pcr, Precip, Runoff)
				
				#-sediment transport
				if self.SedTransFLAG == 1:
					#-read dynamic sediment transport processes mmf
					self.sediment_transport.dynamic_mmf(self, pcr, Runoff, np, G)

			#-INCA
			if self.SedModel == 3:
				#-determine soil erosion
				Sed = self.inca.dynamic(self, pcr, Precip, Q)

			#-SHETRAN
			if self.SedModel == 4:
				#-determine soil erosion
				Sed = self.shetran.dynamic(self, pcr, np, Precip, Q)

			#-DHSVM
			if self.SedModel == 5:
				#-determine soil erosion
				Sed = self.dhsvm.dynamic(self, pcr, np, Precip, Q)

			#-HSPF
			if self.SedModel == 6:
				#-determine soil erosion
				Sed = self.hspf.dynamic(self, pcr, np, Precip, Runoff)


		#-update current date
		self.curdate = self.curdate + self.datetime.timedelta(days=1)

# END OF SPHY CLASS

SPHY = sphy()
timesteps = SPHY.timecalc.timesteps(SPHY)
RunSPHY= pcrm.DynamicFramework(SPHY,lastTimeStep= timesteps,firstTimestep= 1)
RunSPHY.run()

# move tss files to output directory
tssfiles = glob.glob('*.tss')
for i in tssfiles:
	if os.path.exists(SPHY.outpath + i):
		os.remove(SPHY.outpath + i)
	shutil.move(i, SPHY.outpath)

toc = time.time()
dt = toc - tic
print('Simulation succesfully completed in '+str(dt)+' seconds!')
