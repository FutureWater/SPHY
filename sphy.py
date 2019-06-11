# The Spatial Processes in HYdrology (SPHY) model:
# A spatially distributed hydrological model that calculates soil-water and
# cryosphere processes on a cell-by-cell basis.
#
# Copyright (C) 2013  FutureWater
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
#
# Email: info@futurewater.nl

#-Authorship information-###################################################################
__authors__ = "W. Terink, A. Lutz, G. Simons, W. Immerzeel and P. Droogers"
__copyright__ = "FutureWater"
__license__ = "GPL"
__version__ = "2.2"
__email__ = "info@futurewater.nl"
__date__ ='23 December 2017'
############################################################################################

# This model uses the sphy_config.cfg as configuration file.

import time, shutil, os, glob, ConfigParser
import pandas as pd
import pcraster as pcr
import pcraster.framework as pcrm
from pcraster._pcraster import Scalar
import numpy as np

tic = time.clock()

# Read the model configuration file
config = ConfigParser.RawConfigParser()
config.read('sphy_config.cfg')

class sphy(pcrm.DynamicModel):
	def __init__(self):
		# Print model info
		print 	'The Spatial Processes in HYdrology (SPHY) model is ' \
				'developed by Wilco Terink, Wageningen, The Netherlands'
		print   'Version 2.2'
		print ' '

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
		
		# import the required modules
		import datetime, calendar, reporting, timecalc, ET, rootzone, subzone
		from math import pi
		#-standard python modules
		self.datetime = datetime
		self.calendar = calendar
		self.pi = pi
		#-FW defined modules
		self.reporting = reporting
		self.timecalc = timecalc
		self.ET = ET
		self.rootzone = rootzone
		self.subzone = subzone
		del datetime, calendar, pi, reporting, timecalc, ET, rootzone, subzone
		#-import additional modules if required
		if self.GlacFLAG == 1:
			self.SnowFLAG = 1
			self.GroundFLAG = 1
		if self.SnowFLAG == 1:
			import snow # snow melt processes
			self.snow = snow
			del snow
		if self.RoutFLAG == 1:
			import routing # simple routing scheme
			self.routing = routing
			del routing
		if self.LakeFLAG == 1:
			import lakes # import lake module
			self.lakes = lakes
			del lakes
		if self.ResFLAG == 1:
			import reservoirs # import reservoir module
			self.reservoirs = reservoirs
			del reservoirs
		if self.LakeFLAG == 1 or self.ResFLAG == 1:
			import advanced_routing  # overwrite the simple routing scheme
			self.routing = advanced_routing
			del advanced_routing
			self.RoutFLAG = 0
		if self.DynVegFLAG == 1:
			import dynamic_veg # dynamic crop growth using ndvi or kc time-series
			self.dynamic_veg = dynamic_veg
			del dynamic_veg
		if self.GroundFLAG == 1:
			import groundwater # groundwater storage as third storage layer. This is used instead of a fixed bottomflux
			self.groundwater = groundwater
			del groundwater
			
		#-read the input and output directories from the configuration file
		self.inpath = config.get('DIRS', 'inputdir')
		self.outpath = config.get('DIRS', 'outputdir')

		#-set the timing criteria
		sy = config.getint('TIMING', 'startyear')
		sm = config.getint('TIMING', 'startmonth')
		sd = config.getint('TIMING', 'startday')
		ey = config.getint('TIMING', 'endyear')
		em = config.getint('TIMING', 'endmonth')
		ed = config.getint('TIMING', 'endday')
		self.startdate = self.datetime.datetime(sy,sm,sd)
		self.enddate = self.datetime.datetime(ey,em,ed)
		self.dateAfterUpdate = self.startdate - self.datetime.timedelta(days=1)  #-only required for glacier retreat (create dummy value here to introduce the variable) 
		
		#-get start date of first forcing file in forcing directory
		syF = config.getint('TIMING', 'startyear_F')
		smF = config.getint('TIMING', 'startmonth_F')
		sdF = config.getint('TIMING', 'startday_F')
		self.startdateF = self.datetime.datetime(syF, smF, sdF)
		
		#-set the global options
		pcr.setglobaloption('radians')
		#-set the 2000 julian date number
		self.julian_date_2000 = 2451545
		#-set the option to calculate the fluxes in mm for the upstream area
		self.mm_rep_FLAG = config.getint('REPORTING','mm_rep_FLAG')
			
		#-setting clone map
		clonemap = self.inpath + config.get('GENERAL','mask')
		pcr.setclone(clonemap)
		self.clone = pcr.ifthen(pcr.readmap(clonemap), pcr.boolean(1))
		
		self.cellArea = pcr.cellvalue(pcr.cellarea(),1)[0]
		
		#-read general maps
		self.DEM = pcr.readmap(self.inpath + config.get('GENERAL','dem'))
		self.Slope = pcr.readmap(self.inpath + config.get('GENERAL','Slope'))
		self.Locations = pcr.readmap(self.inpath + config.get('GENERAL','locations'))
		
		#-read soil maps
		self.RootFieldMap = pcr.readmap(self.inpath + config.get('SOIL','RootFieldMap'))
		self.RootSatMap = pcr.readmap(self.inpath + config.get('SOIL','RootSatMap'))
		self.RootDryMap = pcr.readmap(self.inpath + config.get('SOIL','RootDryMap'))
		self.RootWiltMap = pcr.readmap(self.inpath + config.get('SOIL','RootWiltMap'))
		self.RootKsat = pcr.readmap(self.inpath + config.get('SOIL','RootKsat'))
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
		if self.GroundFLAG == 0:   # if groundwater module is not used, read seepage and gwl_base
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
		else: # if groundwater module is used, then read the groundwater soil parameters
			pars = ['GwDepth','GwSat','deltaGw','BaseThresh','alphaGw','YieldGw']
			for i in pars:
				try:
					setattr(self, i, pcr.readmap(self.inpath + config.get('GROUNDW_PARS',i)))
				except:
					setattr(self, i, config.getfloat('GROUNDW_PARS',i))
					
		#-calculate soil properties
		self.RootField = self.RootFieldMap * self.RootDepthFlat
		self.RootSat = self.RootSatMap * self.RootDepthFlat
		self.RootDry = self.RootDryMap * self.RootDepthFlat
		self.RootWilt = self.RootWiltMap * self.RootDepthFlat
		self.SubSat = self.SubSatMap * self.SubDepthFlat
		self.SubField = self.SubFieldMap * self.SubDepthFlat
		self.RootTT = (self.RootSat - self.RootField) / self.RootKsat
		self.SubTT = (self.SubSat - self.SubField) / self.SubKsat
		# soil max and soil min for scaling of gwl if groundwater module is not used
		if self.GroundFLAG == 0:
			self.SoilMax = self.RootSat + self.SubSat
			self.SoilMin = self.RootDry + self.SubField

		#-read the crop coefficient table if the dynamic vegetation module is not used
		if self.DynVegFLAG == 0:
			self.KcStatFLAG = config.getint('LANDUSE', 'KCstatic')
			if self.KcStatFLAG == 1:
				#-read land use map and kc table
				self.LandUse = pcr.readmap(self.inpath + config.get('LANDUSE','LandUse'))
				self.kc_table = self.inpath + config.get('LANDUSE','CropFac')
				self.Kc = pcr.lookupscalar(self.kc_table, self.LandUse)
			else:
				#-set the kc map series
				self.Kcmaps = self.inpath + config.get('LANDUSE', 'KC')
		#-Use the dynamic vegetation module
		else:
			#-set the ndvi map series to be read
			self.ndvi = self.inpath + config.get('DYNVEG', 'NDVI')
			#-read the vegetation parameters
			pars = ['NDVImax','NDVImin','NDVIbase','KCmax','KCmin','LAImax','FPARmax','FPARmin']
			for i in pars:
				try:
					setattr(self, i, pcr.readmap(self.inpath + config.get('DYNVEG', i)))
				except:
					setattr(self, i, config.getfloat('DYNVEG', i))
			
		#-read and set glacier maps and parameters if glacier module is used
		if self.GlacFLAG:
			#-Table with glaciers properties (U_ID,MOD_ID,GLAC_ID,MOD_H,GLAC_H,DEBRIS,FRAC_GLAC)
			self.GlacTable = pd.read_csv(os.path.join(self.inpath, config.get('GLACIER', 'GlacTable')))
			cols = pd.DataFrame(columns=['MOD_T', 'GLAC_T', 'Prec_GLAC', 'Rain_GLAC', 'Snow_GLAC', 'PotSnowMelt_GLAC',\
										'ActSnowMelt_GLAC', 'OldSnowStore_GLAC', 'SnowStore_GLAC', 'OldSnowWatStore_GLAC', 'SnowWatStore_GLAC',\
										'MaxSnowWatStore_GLAC', 'OldTotalSnowStore_GLAC', 'TotalSnowStore_GLAC',\
										'SnowR_GLAC', 'GlacMelt', 'GlacR', 'GlacPerc'])
			self.GlacTable = pd.concat([self.GlacTable, cols], axis=1).fillna(0)
			#-sort on MOD_ID column
			self.GlacTable.sort_values(by='MOD_ID', inplace=True)
			#-Read the glacier maps
			self.ModelID = pcr.pcr2numpy(pcr.readmap(os.path.join(self.inpath, config.get('GLACIER', 'ModelID'))), self.MV) #-Nominal
			self.ModelID_1d = self.ModelID.flatten()	#-1 dim array with model cell IDs
			SelModelID = pd.unique(self.GlacTable['MOD_ID'])  #-model id cells for which to extract temperature, precip, etc. (=cells that have glaciers)
			#-Create keys for glacier cells (index in ModelID_1d where cell has glacier)
			self.GlacierKeys= pcr.numpy.ones(self.ModelID_1d.shape)* self.MV
			n = pcr.numpy.arange(self.ModelID_1d.size)
			iCnt= 0
			for ID in SelModelID:
				if ID in self.ModelID_1d:
					key= n[self.ModelID_1d == ID]
					self.GlacierKeys[iCnt] = key
				iCnt += 1
			self.GlacierKeys = self.GlacierKeys[self.GlacierKeys != self.MV]
			self.GlacierKeys = [int(float(x)) for x in self.GlacierKeys] #-make the index an integer
			#-Read the glacier parameters
			pars = ['DDFG','DDFDG','GlacF']
			for	i in pars:
				try:
					setattr(self, i, pcr.readmap(self.inpath + config.get('GLACIER',i)))
				except:
					setattr(self, i, config.getfloat('GLACIER',i))
			#-Lapse rate for temperature
			self.TLapse_table = pd.read_csv(self.inpath + config.get('CLIMATE','TLapse'), header=None, index_col=0, sep=' ', skipinitialspace =True)
			#-Map with glacier IDs
			self.GlacID = pcr.readmap(self.inpath + config.get('GLACIER','GlacID'))
			#-Check if glacier retreat should be calculated
			self.GlacRetreat = config.getint('GLACIER','GlacRetreat')
			if self.GlacRetreat == 1:
				#-Get date for updating the glacier fraction (once a year)
				GlacUpdate = config.get('GLACIER', 'GlacUpdate').split(',')
				self.GlacUpdate = {}
				self.GlacUpdate['month'] = int(GlacUpdate[1])
				self.GlacUpdate['day'] = int(GlacUpdate[0])
					
		#-read and set snow maps and parameters if snow modules are used
		if self.SnowFLAG:
			pars = ['Tcrit','SnowSC','DDFS']
			for	i in pars:
				try:
					setattr(self, i, pcr.readmap(self.inpath + config.get('SNOW',i)))
				except:
					setattr(self, i, config.getfloat('SNOW',i))
		
		#-read and set climate forcing and the calculation of etref
		self.Prec = self.inpath + config.get('CLIMATE','Prec')
		self.Tair = self.inpath + config.get('CLIMATE','Tair')
		self.ETREF_FLAG = config.getint('ETREF','ETREF_FLAG')
		#-determine the use of a given etref time-series or calculate etref using Hargreaves
		if self.ETREF_FLAG == 1:
			self.ETref = self.inpath + config.get('ETREF','ETref')
		else:
			self.Lat = pcr.readmap(self.inpath + config.get('ETREF','Lat'))
			self.Tmax = self.inpath + config.get('ETREF','Tmax')
			self.Tmin = self.inpath + config.get('ETREF','Tmin')
			self.Gsc = config.getfloat('ETREF', 'Gsc')
			import hargreaves
			self.Hargreaves = hargreaves
			del hargreaves
		
		#-read and set routing maps and parameters
		if self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1:
			self.FlowDir = pcr.readmap(self.inpath + config.get('ROUTING','flowdir'))
			try:
				self.kx = pcr.readmap(self.inpath + config.get('ROUTING','kx'))
			except:
				self.kx = config.getfloat('ROUTING','kx')
		
		pcr.setglobaloption('matrixtable')
		#-read lake maps and parameters if lake module is used
		if self.LakeFLAG == 1:
			# nominal map with lake IDs
			self.LakeID = pcr.cover(pcr.readmap(self.inpath + config.get('LAKE','LakeId')), 0)
			# lookup table with function for each lake (exp, 1-order poly, 2-order poly, 3-order poly)
			LakeFunc_Tab = self.inpath + config.get('LAKE', 'LakeFunc')
			# lookup table with Qh-coeficients for each lake
			LakeQH_Tab = self.inpath + config.get('LAKE', 'LakeQH')
			# lookup table with Sh-coeficients for each lake
			LakeSH_Tab = self.inpath + config.get('LAKE', 'LakeSH')
			# lookup table with hS-coeficients for each lake
			LakeHS_Tab = self.inpath + config.get('LAKE', 'LakeHS')
			# create lake coefficient maps
			self.LakeQH_Func = pcr.lookupnominal(LakeFunc_Tab, 1, self.LakeID)
			self.LakeSH_Func = pcr.lookupnominal(LakeFunc_Tab, 2, self.LakeID)
			self.LakeHS_Func = pcr.lookupnominal(LakeFunc_Tab, 3, self.LakeID)
			# Read QH coefficients
			self.LakeQH_exp_a = pcr.lookupscalar(LakeQH_Tab, 1, self.LakeID) 
			self.LakeQH_exp_b = pcr.lookupscalar(LakeQH_Tab, 2, self.LakeID)
			self.LakeQH_pol_b = pcr.lookupscalar(LakeQH_Tab, 3, self.LakeID)
			self.LakeQH_pol_a1 = pcr.lookupscalar(LakeQH_Tab, 4, self.LakeID)
			self.LakeQH_pol_a2 = pcr.lookupscalar(LakeQH_Tab, 5, self.LakeID)
			self.LakeQH_pol_a3 = pcr.lookupscalar(LakeQH_Tab, 6, self.LakeID)
			# Read SH coefficients
			self.LakeSH_exp_a = pcr.lookupscalar(LakeSH_Tab, 1, self.LakeID) 
			self.LakeSH_exp_b = pcr.lookupscalar(LakeSH_Tab, 2, self.LakeID)
			self.LakeSH_pol_b = pcr.lookupscalar(LakeSH_Tab, 3, self.LakeID)
			self.LakeSH_pol_a1 = pcr.lookupscalar(LakeSH_Tab, 4, self.LakeID)
			self.LakeSH_pol_a2 = pcr.lookupscalar(LakeSH_Tab, 5, self.LakeID)
			self.LakeSH_pol_a3 = pcr.lookupscalar(LakeSH_Tab, 6, self.LakeID)
			# Read HS coefficients
			self.LakeHS_exp_a = pcr.lookupscalar(LakeHS_Tab, 1, self.LakeID) 
			self.LakeHS_exp_b = pcr.lookupscalar(LakeHS_Tab, 2, self.LakeID)
			self.LakeHS_pol_b = pcr.lookupscalar(LakeHS_Tab, 3, self.LakeID)
			self.LakeHS_pol_a1 = pcr.lookupscalar(LakeHS_Tab, 4, self.LakeID)
			self.LakeHS_pol_a2 = pcr.lookupscalar(LakeHS_Tab, 5, self.LakeID)
			self.LakeHS_pol_a3 = pcr.lookupscalar(LakeHS_Tab, 6, self.LakeID)
			#-read water level maps and parameters if available
			try:
				self.UpdateLakeLevel = pcr.readmap(self.inpath + config.get('LAKE','updatelakelevel'))
				self.LLevel = self.inpath + config.get('LAKE','LakeFile')
				print 'measured lake levels will be used to update lake storage'
			except:
				pass
			
		#-read reservior maps and parameters if reservoir module is used
		if self.ResFLAG == 1:
			# nominal map with reservoir IDs
			self.ResID = pcr.cover(pcr.readmap(self.inpath + config.get('RESERVOIR','ResId')), 0)
			# lookup table with operational scheme to use (simple or advanced)
			ResFunc_Tab = self.inpath + config.get('RESERVOIR', 'ResFuncStor')
			# Reservoir function
			self.ResFunc = pcr.cover(pcr.lookupscalar(ResFunc_Tab, 1, self.ResID), 0)
			try:
				# lookup table with coefficients for simple reservoirs
				ResSimple_Tab = self.inpath + config.get('RESERVOIR', 'ResSimple')
				# Read coefficients for simple reservoirs
				self.ResKr = pcr.lookupscalar(ResSimple_Tab, 1, self.ResID)
				self.ResSmax = pcr.lookupscalar(ResSimple_Tab, 2, self.ResID) * 10**6 # convert to m3
				self.ResSimple = True
			except:
				self.ResSimple = False
			try:
				# lookup table with coefficients for advanced reservoirs
				ResAdvanced_Tab = self.inpath + config.get('RESERVOIR', 'ResAdv')
				# Read coefficients for advanced reservoirs
				self.ResEVOL = pcr.lookupscalar(ResAdvanced_Tab, 1, self.ResID) * 10**6 # convert to m3
				self.ResPVOL = pcr.lookupscalar(ResAdvanced_Tab, 2, self.ResID) * 10**6 # convert to m3
				self.ResMaxFl = pcr.lookupscalar(ResAdvanced_Tab, 3, self.ResID) * 10**6 # convert to m3/d
				self.ResDemFl = pcr.lookupscalar(ResAdvanced_Tab, 4, self.ResID) * 10**6 # convert to m3/d
				self.ResFlStart = pcr.lookupscalar(ResAdvanced_Tab, 5, self.ResID)
				self.ResFlEnd = pcr.lookupscalar(ResAdvanced_Tab, 6, self.ResID)
				self.ResAdvanced = True
			except:
				self.ResAdvanced = False
	
	def initial(self):
		#-initial section
		#-get the correct forcing file number, depending on the start date of your simulation
		#-and the start date of the first forcing file in your forcing directory.
		self.counter = (self.startdate - self.startdateF).days
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
			#-initial canopy storage
			self.Scanopy = 0
			#-initial ndvi if first map is not provided
			self.ndviOld = pcr.scalar((self.NDVImax + self.NDVImin)/2)
		elif self.KcStatFLAG == 0:
			#-set initial kc value to one, if kc map is not available for first timestep
			self.KcOld = pcr.scalar(1)
		
		#-initial groundwater properties
		if self.GroundFLAG == 1:
			#-initial groundwater recharge
			try:
				self.GwRecharge = config.getfloat('GROUNDW_INIT','GwRecharge')
			except:
				self.GwRecharge = pcr.readmap(self.inpath + config.get('GROUNDW_INIT','GwRecharge'))
			#-initial baseflow
			try:
				self.BaseR = config.getfloat('GROUNDW_INIT','BaseR')
			except:
				self.BaseR = pcr.readmap(self.inpath + config.get('GROUNDW_INIT','BaseR'))
			#-initial groundwater storage
			try:
				self.Gw = config.getfloat('GROUNDW_INIT','Gw')
			except:
				self.Gw = pcr.readmap(self.inpath + config.get('GROUNDW_INIT','Gw'))
			#-initial groundwater level
			try:
				self.H_gw = config.getfloat('GROUNDW_INIT','H_gw')
			except:
				self.H_gw = pcr.readmap(self.inpath + config.get('GROUNDW_INIT','H_gw'))
			self.H_gw = pcr.max((self.RootDepthFlat + self.SubDepthFlat + self.GwDepth)/1000 - self.H_gw, 0)	
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
			try:
				self.SnowStore = config.getfloat('SNOW_INIT','SnowIni')
			except:
				self.SnowStore = pcr.readmap(self.inpath + config.get('SNOW_INIT','SnowIni'))			
			#-initial water stored in snowpack
			try:
				self.SnowWatStore = config.getfloat('SNOW_INIT','SnowWatStore')
			except:
				self.SnowWatStore = pcr.readmap(self.inpath + config.get('SNOW_INIT','SnowWatStore'))
			self.TotalSnowStore = self.SnowStore + self.SnowWatStore
		else:
			self.SnowStore = pcr.scalar(0)

		#-initial glacier properties
		if self.GlacFLAG:
			self.GlacTable.set_index('MOD_ID', inplace=True)
			#-initial SnowStore on glacier
			try:
				SnowStore_1d = pcr.pcr2numpy(self.SnowStore, self.MV).flatten()
				S = pd.DataFrame(data={'SnowStore_GLAC': SnowStore_1d[self.GlacierKeys]}, index=self.ModelID_1d[self.GlacierKeys])
				self.GlacTable.update(S)
				S = None; SnowStore_1d = None; del S, SnowStore_1d
			except:
				self.GlacTable['SnowStore_GLAC'] = self.SnowStore
			#-initial SnowWatstore on glacier
			try:
				SnowWatStore_1d = pcr.pcr2numpy(self.SnowWatStore, self.MV).flatten()
				S = pd.DataFrame(data={'SnowWatStore_GLAC': SnowWatStore_1d[self.GlacierKeys]}, index=self.ModelID_1d[self.GlacierKeys])
				self.GlacTable.update(S)
				S = None; SnowWatStore_1d = None; del S, SnowWatStore_1d
			except:
				self.GlacTable['SnowWatStore_GLAC'] = self.SnowWatStore
			self.GlacTable['TotalSnowStore_GLAC']  = self.GlacTable['SnowStore_GLAC'] + self.GlacTable['SnowWatStore_GLAC']
			#-Create table for accumulating glacier melt over certain period
			self.GlacTable['AccuGlacMelt'] = 0.
			#-initialize the glacier fraction
			self.GlacFrac = pcr.numpy.zeros(self.ModelID_1d.shape)
			self.GlacFrac[self.GlacierKeys] = self.GlacTable.groupby(self.GlacTable.index).sum()['FRAC_GLAC']
			self.GlacFrac = self.GlacFrac.reshape(self.ModelID.shape)
			self.GlacFrac = pcr.numpy2pcr(Scalar, self.GlacFrac, self.MV)
			self.GlacFrac = pcr.ifthen(self.clone, self.GlacFrac)  #-only use values where clone is True
			pcr.report(self.GlacFrac, self.outpath + 'GlacFrac_' + self.curdate.strftime('%Y%m%d') + '.map')
			# 1-D Masks for debris and clean ice
			self.CImask = self.GlacTable['DEBRIS'] == 0
			self.DBmask = pcr.numpy.invert(self.CImask)
			#-If dynamic glacier retreat is on, then make spatial maps of initial glacier area, volume, and depth
			if self.GlacRetreat == 1:
				#-Calculate the  area and depth per model ID
				GlacTable_MODid = self.GlacTable.loc[:,['FRAC_GLAC', 'ICE_DEPTH']]
				GlacTable_MODid['ICE_DEPTH'] = GlacTable_MODid['ICE_DEPTH'] * GlacTable_MODid['FRAC_GLAC']
				GlacTable_MODid['AREA'] = GlacTable_MODid['FRAC_GLAC'] * self.cellArea
				GlacTable_MODid = GlacTable_MODid.groupby(GlacTable_MODid.index).sum()
				#-Report pcraster map of glacier area
				glacArea = pcr.numpy.zeros(self.ModelID_1d.shape)
				glacArea[self.GlacierKeys] = GlacTable_MODid['AREA']
				glacArea = glacArea.reshape(self.ModelID.shape)
				glacArea = pcr.numpy2pcr(Scalar, glacArea, self.MV)
				glacArea = pcr.ifthen(self.clone, glacArea)  #-only use values where clone is True
				pcr.report(glacArea, self.outpath + 'GlacArea_' + self.curdate.strftime('%Y%m%d') + '.map')
				#-Report pcraster map of glacier depth
				iceDepth = pcr.numpy.zeros(self.ModelID_1d.shape)
				iceDepth[self.GlacierKeys] = GlacTable_MODid['ICE_DEPTH']
				iceDepth = iceDepth.reshape(self.ModelID.shape)
				iceDepth = pcr.numpy2pcr(Scalar, iceDepth, self.MV)
				iceDepth = pcr.ifthen(self.clone, iceDepth)  #-only use values where clone is True
				pcr.report(iceDepth, self.outpath + 'iceDepth_' + self.curdate.strftime('%Y%m%d') + '.map')
				#-Report pcraster map of glacier volume
				pcr.report(iceDepth * glacArea, self.outpath + 'vIce_' + self.curdate.strftime('%Y%m%d') + '.map')
				#-Delete variables that are not needed
				glacArea = None; del glacArea; #iceDepth = None; del iceDepth; 
				GlacTable_MODid = None; del GlacTable_MODid
		else:
			self.GlacFrac = 0
			
		#-initial routed total runoff and of individual components
		if self.RoutFLAG == 1 or self.LakeFLAG==1 or self.ResFLAG==1:
			#-initial routed total runoff
			try:
				self.QRAold = config.getfloat('ROUT_INIT','QRA_init')
			except:
				try:
					self.QRAold = pcr.readmap(self.inpath + config.get('ROUT_INIT','QRA_init'))
				except:
					self.QRAold = 0
			#-initial routed runoff	for the individual components
			pars = ['RainRA','SnowRA','GlacRA','BaseRA']
			for i in pars:
				try:
					setattr(self, i + 'old', pcr.readmap(self.inpath + config.get('ROUT_INIT', i + '_init')))
					setattr(self, i + '_FLAG', True)
				except:
					try:
						setattr(self, i + 'old', config.getfloat('ROUT_INIT', i + '_init'))
						setattr(self, i + '_FLAG', True)
					except:
						setattr(self, i + '_FLAG', False)
		
		#-initial storage in lakes and reservoirs
		if self.LakeFLAG == 1 or self.ResFLAG == 1:
			#-Read initial storages from table/reservoir file
			if self.LakeFLAG == 1:
				LakeStor_Tab = self.inpath + config.get('LAKE', 'LakeStor')
				self.StorRES = pcr.cover(pcr.lookupscalar(LakeStor_Tab, 1, self.LakeID), 0) * 10**6  # convert to m3
				#-Qfrac for lake cells should be zero, else 1
				self.QFRAC = pcr.ifthenelse(self.LakeID != 0, pcr.scalar(0), 1)
			if self.ResFLAG == 1:
				ResStor_Tab = self.inpath + config.get('RESERVOIR', 'ResFuncStor')
				ResStor = pcr.cover(pcr.lookupscalar(ResStor_Tab, 2, self.ResID), 0) * 10**6  # convert to m3
				try:
					self.StorRES = self.StorRES + ResStor
					#-Qfrac for reservoir cells should be zero, else 1
					self.QFRAC = pcr.ifthenelse(self.ResID != 0, pcr.scalar(0), self.QFRAC)
				except:
					self.StorRES = ResStor
					#-Qfrac for reservoir cells should be zero, else 1
					self.QFRAC = pcr.ifthenelse(self.ResID != 0, pcr.scalar(0), 1)
			
			#-initial storage in lakes/reservoirs of individual flow components
			pars = ['RainRA','SnowRA','GlacRA','BaseRA']
			for i in pars:
				column = pars.index(i)  # identify column to be read from lake or reservoir table
				try: #-try to sum the storages read from the lake and reservoir tables if both thse modules are used
					setattr(self, i + 'stor', (pcr.cover(pcr.lookupscalar(LakeStor_Tab, column + 2, self.LakeID), 0) + \
											pcr.cover(pcr.lookupscalar(ResStor_Tab, column + 3, self.ResID), 0)) * 10**6)
					if eval('self.' + i + '_FLAG'):
						setattr(self, i + '_FLAG', True)
					else:
						setattr(self, i + '_FLAG', False)
				except:
					try: #-try to read the storages from the lake table
						setattr(self, i + 'stor', pcr.cover(pcr.lookupscalar(LakeStor_Tab, column + 2, self.LakeID), 0) * 10**6)
						if eval('self.' + i + '_FLAG'):
							setattr(self, i + '_FLAG', True)
						else:
							setattr(self, i + '_FLAG', False)
					except: #-try to read the storages from the reservoir table
						try:
							setattr(self, i + 'stor', pcr.cover(pcr.lookupscalar(ResStor_Tab, column + 3, self.ResID), 0) * 10**6)
							if eval('self.' + i + '_FLAG'):
								setattr(self, i + '_FLAG', True)
							else:
								setattr(self, i + '_FLAG', False)
						except:
							setattr(self, i + '_FLAG', False)
		
		#-Initial values for reporting and setting of time-series
		#-set time-series reporting for mm flux from upstream area for prec and eta 
		if self.mm_rep_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1): 
			self.PrecSubBasinTSS = pcrm.TimeoutputTimeseries("PrecSubBasinTSS", self, self.Locations, noHeader=False)
			self.ETaSubBasinTSS = pcrm.TimeoutputTimeseries("ETaSubBasinTSS", self, self.Locations, noHeader=False)
		if self.GlacFLAG == 1:
			pars = ['wbal','GWL','TotPrec','TotPrecE','TotInt','TotRain','TotETpot','TotETact','TotSnow','TotSnowMelt','TotGlacMelt','TotGlacMelt','TotRootR','TotRootD','TotRootP',\
				'TotSubP','TotCapR','TotGlacPerc','TotGwRecharge','TotRainR','TotBaseR','TotSnowR','TotGlacR','TotR','RainRAtot','SnowRAtot','GlacRAtot','BaseRAtot','QallRAtot']
			#-set time-series reporting for mm fluxes from upstream area if glacier and routing/reservoir modules are used
			if self.mm_rep_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
				self.GMeltSubBasinTSS = pcrm.TimeoutputTimeseries("GMeltSubBasinTSS", self, self.Locations, noHeader=False)
				self.QSNOWSubBasinTSS = pcrm.TimeoutputTimeseries("QSNOWSubBasinTSS", self, self.Locations, noHeader=False)
				self.QRAINSubBasinTSS = pcrm.TimeoutputTimeseries("QRAINSubBasinTSS", self, self.Locations, noHeader=False)
				self.QGLACSubBasinTSS = pcrm.TimeoutputTimeseries("QGLACSubBasinTSS", self, self.Locations, noHeader=False)
				self.QBASFSubBasinTSS = pcrm.TimeoutputTimeseries("QBASFSubBasinTSS", self, self.Locations, noHeader=False)
				self.QTOTSubBasinTSS = pcrm.TimeoutputTimeseries("QTOTSubBasinTSS", self, self.Locations, noHeader=False)
			#-Initialize reporting per glacier ID
			self.GlacID_flag = config.getint('REPORTING', 'GlacID_flag')
			if self.GlacID_flag:
				self.GlacVars = config.get('REPORTING', 'GlacID_report').split(',')  #-get the variables to report
				self.glacid = sorted(self.GlacTable['GLAC_ID'].unique())  #-get the unique glacier ids
				self.GlacID_memerror = config.getint('REPORTING', 'GlacID_memerror') #-check if a memory error occured while runnning the model
				if self.GlacID_memerror == 0:  #-with no memory error we can store the pandas dataframes in the computer's memory
					drange = pd.date_range(self.startdate, self.enddate, freq='D')
					for p in self.GlacVars: #-make panda dataframes for each variable to report
						setattr(self, p + '_Table', pd.DataFrame(index = drange, columns=self.glacid,dtype=np.float32))  #-create table for each variable to report
		elif self.SnowFLAG == 1:
			if self.GroundFLAG == 1:		
				pars = ['wbal','GWL','TotPrec','TotPrecE','TotInt','TotRain','TotETpot','TotETact','TotSnow','TotSnowMelt','TotRootR','TotRootD','TotRootP',\
					'TotSubP','TotCapR','TotGwRecharge','TotRainR','TotBaseR','TotSnowR','TotR','RainRAtot','SnowRAtot','BaseRAtot','QallRAtot']
				#-set time-series reporting for mm fluxes from upstream area if snow, groundwater and routing/reservoir modules are used
				if self.mm_rep_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
					self.QSNOWSubBasinTSS = pcrm.TimeoutputTimeseries("QSNOWSubBasinTSS", self, self.Locations, noHeader=False)
					self.QRAINSubBasinTSS = pcrm.TimeoutputTimeseries("QRAINSubBasinTSS", self, self.Locations, noHeader=False)
					self.QBASFSubBasinTSS = pcrm.TimeoutputTimeseries("QBASFSubBasinTSS", self, self.Locations, noHeader=False)
					self.QTOTSubBasinTSS = pcrm.TimeoutputTimeseries("QTOTSubBasinTSS", self, self.Locations, noHeader=False)
			else:
				pars = ['wbal','GWL','TotPrec','TotPrecE','TotInt','TotRain','TotETpot','TotETact','TotSnow','TotSnowMelt','TotRootR','TotRootD','TotRootP',\
					'TotSubD','TotCapR','TotSeep','TotRainR','TotSnowR','TotR','RainRAtot','SnowRAtot','BaseRAtot','QallRAtot']
				#-set time-series reporting for mm fluxes from upstream area if snow and routing/reservoir modules are used
				if self.mm_rep_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
					self.SeepSubBasinTSS = pcrm.TimeoutputTimeseries("SeepSubBasinTSS", self, self.Locations, noHeader=False)
					self.QSNOWSubBasinTSS = pcrm.TimeoutputTimeseries("QSNOWSubBasinTSS", self, self.Locations, noHeader=False)
					self.QRAINSubBasinTSS = pcrm.TimeoutputTimeseries("QRAINSubBasinTSS", self, self.Locations, noHeader=False)
					self.QBASFSubBasinTSS = pcrm.TimeoutputTimeseries("QBASFSubBasinTSS", self, self.Locations, noHeader=False)
					self.QTOTSubBasinTSS = pcrm.TimeoutputTimeseries("QTOTSubBasinTSS", self, self.Locations, noHeader=False)
		else:
			if self.GroundFLAG == 1:
				pars = ['wbal','GWL','TotPrec','TotPrecE','TotInt','TotRain','TotETpot','TotETact','TotRootR','TotRootD','TotRootP',\
					'TotSubP','TotCapR','TotGwRecharge','TotRainR','TotBaseR','TotR','RainRAtot','BaseRAtot','QallRAtot']
				#-set time-series reporting for mm fluxes from upstream area if groundwater and routing/reservoir modules are used
				if self.mm_rep_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
					self.QRAINSubBasinTSS = pcrm.TimeoutputTimeseries("QRAINSubBasinTSS", self, self.Locations, noHeader=False)
					self.QBASFSubBasinTSS = pcrm.TimeoutputTimeseries("QBASFSubBasinTSS", self, self.Locations, noHeader=False)
					self.QTOTSubBasinTSS = pcrm.TimeoutputTimeseries("QTOTSubBasinTSS", self, self.Locations, noHeader=False)
			else:
				pars = ['wbal','GWL','TotPrec','TotPrecE','TotInt','TotRain','TotETpot','TotETact','TotRootR','TotRootD','TotRootP',\
					'TotSubD','TotCapR','TotSeep','TotRainR','TotR','RainRAtot','BaseRAtot','QallRAtot']
				#-set time-series reporting for mm fluxes from upstream area if routing/reservoir modules are used
				if self.mm_rep_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
					self.SeepSubBasinTSS = pcrm.TimeoutputTimeseries("SeepSubBasinTSS", self, self.Locations, noHeader=False)
					self.QRAINSubBasinTSS = pcrm.TimeoutputTimeseries("QRAINSubBasinTSS", self, self.Locations, noHeader=False)
					self.QBASFSubBasinTSS = pcrm.TimeoutputTimeseries("QBASFSubBasinTSS", self, self.Locations, noHeader=False)
					self.QTOTSubBasinTSS = pcrm.TimeoutputTimeseries("QTOTSubBasinTSS", self, self.Locations, noHeader=False)
		#-remove routing output from reported list of parameters if these modules are not used			
		if self.RoutFLAG == 0 and self.ResFLAG == 0 and self.LakeFLAG == 0:
				rpars = ['RainRAtot','SnowRAtot','GlacRAtot','BaseRAtot','QallRAtot']
				for i in rpars:
					try:
						j = pars.index(i)
						del pars[j]
					except:
						pass
		#-set reporting options and read initial values
		for i in pars:
			mapoutops = config.get('REPORTING', i+'_mapoutput')
			TSoutops = config.get('REPORTING', i+'_TSoutput')
			if mapoutops == 'NONE' and TSoutops == 'NONE':
				print i + ' will NOT be reported'
			else:
				print i + ' will be reported'
				fname = config.get('REPORTING', i+'_fname')
				setattr(self, i+'_fname', fname)
				setattr(self, i, 0.)
				if mapoutops != 'NONE':
					mapoutops = mapoutops.split(",")
					for j in mapoutops:
						if j == 'D':
							setattr(self, i+'_Day', eval('self.'+i))
							setattr(self, i+'_Day_map', 1)
						elif j == 'M':
							setattr(self, i+'_Month', eval('self.'+i))
							setattr(self, i+'_Month_map', 1)
						elif j == 'Y':
							setattr(self, i+'_Year', eval('self.'+i))
							setattr(self, i+'_Year_map', 1)
						else:
							setattr(self, i+'_Final', eval('self.'+i))
							setattr(self, i+'_Final_map', 1)
				if TSoutops != 'NONE':
					TSoutops = TSoutops.split(",")
					for j in TSoutops:
						if j == 'D':
							setattr(self, i+'_Day', eval('self.'+i))
							setattr(self, i+'_DayTS', eval('pcrm.TimeoutputTimeseries("'+fname+'DTS'+'", self, self.Locations, noHeader=False)'))
						elif j == 'M':
							setattr(self, i+'_Month', eval('self.'+i))
							setattr(self, i+'_MonthTS', eval('pcrm.TimeoutputTimeseries("'+fname+'MTS'+'", self, self.Locations, noHeader=False)'))
						elif j == 'Y':
							setattr(self, i+'_Year', eval('self.'+i))
							setattr(self, i+'_YearTS', eval('pcrm.TimeoutputTimeseries("'+fname+'YTS'+'", self, self.Locations, noHeader=False)'))
		
		#-set reporting of water balances for lakes
		if self.LakeFLAG == 1 and config.getint('REPORTING', 'Lake_wbal') ==1:
			self.LakeInTSS = pcrm.TimeoutputTimeseries("LakeInTSS", self, self.LakeID, noHeader=True)
			self.LakeOutTSS = pcrm.TimeoutputTimeseries("LakeOutTSS", self, self.LakeID, noHeader=True)
			self.LakeStorTSS = pcrm.TimeoutputTimeseries("LakeStorTSS", self, self.LakeID, noHeader=True)
			if self.RainRA_FLAG==1: #-set reporting of water balances for individual components
				self.LakeRainInTSS = pcrm.TimeoutputTimeseries("LakeRainInTSS", self, self.LakeID, noHeader=True)
				self.LakeRainOutTSS = pcrm.TimeoutputTimeseries("LakeRainOutTSS", self, self.LakeID, noHeader=True)
				self.LakeRainStorTSS = pcrm.TimeoutputTimeseries("LakeRainStorTSS", self, self.LakeID, noHeader=True)
			if self.SnowRA_FLAG==1:
				self.LakeSnowInTSS = pcrm.TimeoutputTimeseries("LakeSnowInTSS", self, self.LakeID, noHeader=True)
				self.LakeSnowOutTSS = pcrm.TimeoutputTimeseries("LakeSnowOutTSS", self, self.LakeID, noHeader=True)
				self.LakeSnowStorTSS = pcrm.TimeoutputTimeseries("LakeSnowStorTSS", self, self.LakeID, noHeader=True)
			if self.GlacRA_FLAG==1:	
				self.LakeGlacInTSS = pcrm.TimeoutputTimeseries("LakeGlacInTSS", self, self.LakeID, noHeader=True)
				self.LakeGlacOutTSS = pcrm.TimeoutputTimeseries("LakeGlacOutTSS", self, self.LakeID, noHeader=True)
				self.LakeGlacStorTSS = pcrm.TimeoutputTimeseries("LakeGlacStorTSS", self, self.LakeID, noHeader=True)
			if self.BaseRA_FLAG==1:
				self.LakeBaseInTSS = pcrm.TimeoutputTimeseries("LakeBaseInTSS", self, self.LakeID, noHeader=True)
				self.LakeBaseOutTSS = pcrm.TimeoutputTimeseries("LakeBaseOutTSS", self, self.LakeID, noHeader=True)
				self.LakeBaseStorTSS = pcrm.TimeoutputTimeseries("LakeBaseStorTSS", self, self.LakeID, noHeader=True)
		#-set reporting of water balances for reservoirs
		if self.ResFLAG == 1 and config.getint('REPORTING', 'Res_wbal') == 1:
			self.ResInTSS = pcrm.TimeoutputTimeseries("ResInTSS", self, self.ResID, noHeader=True)
			self.ResOutTSS = pcrm.TimeoutputTimeseries("ResOutTSS", self, self.ResID, noHeader=True)
			self.ResStorTSS = pcrm.TimeoutputTimeseries("ResStorTSS", self, self.ResID, noHeader=True)
			if self.RainRA_FLAG==1: #-set reporting of water balances for individual components
				self.ResRainInTSS = pcrm.TimeoutputTimeseries("ResRainInTSS", self, self.ResID, noHeader=True)
				self.ResRainOutTSS = pcrm.TimeoutputTimeseries("ResRainOutTSS", self, self.ResID, noHeader=True)
				self.ResRainStorTSS = pcrm.TimeoutputTimeseries("ResRainStorTSS", self, self.ResID, noHeader=True)
			if self.SnowRA_FLAG==1:
				self.ResSnowInTSS = pcrm.TimeoutputTimeseries("ResSnowInTSS", self, self.ResID, noHeader=True)
				self.ResSnowOutTSS = pcrm.TimeoutputTimeseries("ResSnowOutTSS", self, self.ResID, noHeader=True)
				self.ResSnowStorTSS = pcrm.TimeoutputTimeseries("ResSnowStorTSS", self, self.ResID, noHeader=True)
			if self.GlacRA_FLAG==1:	
				self.ResGlacInTSS = pcrm.TimeoutputTimeseries("ResGlacInTSS", self, self.ResID, noHeader=True)
				self.ResGlacOutTSS = pcrm.TimeoutputTimeseries("ResGlacOutTSS", self, self.ResID, noHeader=True)
				self.ResGlacStorTSS = pcrm.TimeoutputTimeseries("ResGlacStorTSS", self, self.ResID, noHeader=True)
			if self.BaseRA_FLAG==1:
				self.ResBaseInTSS = pcrm.TimeoutputTimeseries("ResBaseInTSS", self, self.ResID, noHeader=True)
				self.ResBaseOutTSS = pcrm.TimeoutputTimeseries("ResBaseOutTSS", self, self.ResID, noHeader=True)
				self.ResBaseStorTSS = pcrm.TimeoutputTimeseries("ResBaseStorTSS", self, self.ResID, noHeader=True)
				
		#-WATER BALANCE
		self.oldRootWater = self.RootWater
		self.oldSubWater = self.SubWater
		if self.GroundFLAG:
			self.oldGw = self.Gw
		if self.GlacFLAG and self.GlacRetreat==1:
			self.oldIceDepth = iceDepth * 1000 # in mm
		self.wbalTSS = pcrm.TimeoutputTimeseries("wbalTSS", self, self.Locations, noHeader=True)
		self.wbalTotTSS = pcrm.TimeoutputTimeseries("wbalTotTSS", self, self.Locations, noHeader=True)
		self.waterbalanceTot = pcr.scalar(0.)

	def dynamic(self):
		self.counter+=1
		print str(self.curdate.day)+'-'+str(self.curdate.month)+'-'+str(self.curdate.year)+'  t = '+str(self.counter)
		
		# Snow and rain fraction settings for non-glacier part of model cel
		SnowFrac = pcr.ifthenelse(self.SnowStore > 0, pcr.scalar(1 - self.GlacFrac), 0)
		RainFrac = pcr.ifthenelse(self.SnowStore == 0, pcr.scalar(1 - self.GlacFrac), 0)
		
		#-Read the precipitation time-series
		Precip = pcr.readmap(pcrm.generateNameT(self.Prec, self.counter))
		#-Report Precip  
		self.reporting.reporting(self, pcr, 'TotPrec', Precip)
		
		#-Temperature and reference evapotranspiration
		Temp = pcr.readmap(pcrm.generateNameT(self.Tair, self.counter))
		if self.ETREF_FLAG == 0:
			TempMax = pcr.readmap(pcrm.generateNameT(self.Tmax, self.counter))
			TempMin = pcr.readmap(pcrm.generateNameT(self.Tmin, self.counter))
			ETref = self.Hargreaves.Hargreaves(pcr, self.Hargreaves.extrarad(self, pcr), Temp, TempMax, TempMin)
		else:
			ETref = pcr.readmap(pcrm.generateNameT(self.ETref, self.counter))

		#-Interception and effective precipitation
		#-Update canopy storage
		if self.DynVegFLAG == 1:
			#-try to read the ndvi map series. If not available, then use ndvi old
			try:
				ndvi = pcr.readmap(pcrm.generateNameT(self.ndvi, self.counter))
				self.ndviOld = ndvi
			except:
				ndvi = self.ndviOld
			#-fill missing ndvi values with ndvi base
			ndvi = pcr.ifthenelse(pcr.defined(ndvi) == 1, ndvi, self.NDVIbase)
			#-calculate the vegetation parameters
			vegoutput = self.dynamic_veg.Veg_function(pcr, ndvi, self.FPARmax, self.FPARmin, self.LAImax, self.NDVImin, self.NDVImax, self.KCmin, self.KCmax)
			#-Kc
			self.Kc = vegoutput[0]
			#-Update canopy storage
			self.Scanopy = self.Scanopy + Precip
			#-interception and effective precipitation
			intercep = self.dynamic_veg.Inter_function(pcr, self.Scanopy, vegoutput[1], ETref)
			#-interception
			Int = intercep[0]
			#-report interception
			self.reporting.reporting(self, pcr, 'TotInt', Int * (1-self.GlacFrac))
			#-effective precipitation for model cel
			Precip = intercep[1] * (1-self.GlacFrac) + Precip * self.GlacFrac
			#-Report effective precipitation
			self.reporting.reporting(self, pcr, 'TotPrecE', Precip)
			#-canopy storage
			self.Scanopy = intercep[2]
		elif self.KcStatFLAG == 0:
			#-Try to read the KC map series
			try:
				self.Kc = pcr.readmap(pcrm.generateNameT(self.Kcmaps, self.counter))
				self.KcOld = self.Kc
			except:
				self.Kc = self.KcOld
		#-report mm effective precipitation for sub-basin averages		
		if self.mm_rep_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
			self.PrecSubBasinTSS.sample(pcr.catchmenttotal(Precip, self.FlowDir) / pcr.catchmenttotal(1, self.FlowDir))

		#-Snow, rain, and glacier calculations for glacier fraction of cell
		if self.GlacFLAG:
			#-1 dim array of Tavg map
			T_1d = pcr.pcr2numpy(Temp, self.MV).flatten()  	
			T = pd.DataFrame(data={'MOD_T': T_1d[self.GlacierKeys]}, index=self.ModelID_1d[self.GlacierKeys])
			#-update table with model cel temperature
			self.GlacTable.update(T)
			T = None; T_1d = None; del T, T_1d
			#-lapse temperature for glaciers
			self.GlacTable['GLAC_T'] = self.GlacTable['MOD_T'] - (self.GlacTable['MOD_H'] - self.GlacTable['GLAC_H']) * float(self.TLapse_table.loc[self.curdate.month])
			#-1 dim array of Precip map
			P_1d = pcr.pcr2numpy(Precip, self.MV).flatten()  	
			P = pd.DataFrame(data={'Prec_GLAC': P_1d[self.GlacierKeys]}, index=self.ModelID_1d[self.GlacierKeys])
			#-update table with model cel precipitation
			self.GlacTable.update(P)
			P = None; P_1d = None; del P, P_1d
			#-Snow and rain differentiation
			self.GlacTable['Rain_GLAC'] = 0; self.GlacTable['Snow_GLAC'] = 0; 
			mask = self.GlacTable['GLAC_T'] >= self.Tcrit
			self.GlacTable.loc[mask, 'Rain_GLAC'] = self.GlacTable.loc[mask, 'Prec_GLAC'] 
			self.GlacTable.loc[pcr.numpy.invert(mask), 'Snow_GLAC'] = self.GlacTable.loc[pcr.numpy.invert(mask), 'Prec_GLAC']
			#-Set the melting temperature (=>0)
			Tmelt = pcr.numpy.maximum(self.GlacTable['GLAC_T'], 0)
			#-Snow melt
			self.GlacTable['PotSnowMelt_GLAC'] = Tmelt * self.DDFS
			self.GlacTable['ActSnowMelt_GLAC'] = pcr.numpy.minimum(self.GlacTable['SnowStore_GLAC'], self.GlacTable['PotSnowMelt_GLAC'])
			#-Update snow store
			self.GlacTable['OldSnowStore_GLAC'] = self.GlacTable['SnowStore_GLAC'] 
			self.GlacTable['SnowStore_GLAC'] = self.GlacTable['SnowStore_GLAC'] + self.GlacTable['Snow_GLAC'] - self.GlacTable['ActSnowMelt_GLAC']
			self.GlacTable.loc[self.GlacTable['GLAC_T'] < 0., 'SnowStore_GLAC'] = self.GlacTable.loc[self.GlacTable['GLAC_T'] < 0., 'SnowStore_GLAC'] + \
				self.GlacTable.loc[self.GlacTable['GLAC_T'] < 0., 'SnowWatStore_GLAC']
			#-Caclulate the maximum amount of water that can be stored in snowwatstore
			self.GlacTable['MaxSnowWatStore_GLAC'] = self.SnowSC * self.GlacTable['SnowStore_GLAC']
			self.GlacTable['OldSnowWatStore_GLAC'] = self.GlacTable['SnowWatStore_GLAC']
			#-Calculate the actual amount of water stored in snowwatstore
			self.GlacTable['SnowWatStore_GLAC'] = pcr.numpy.minimum(self.GlacTable['MaxSnowWatStore_GLAC'], self.GlacTable['SnowWatStore_GLAC'] +\
				self.GlacTable['ActSnowMelt_GLAC'] + self.GlacTable['Rain_GLAC'])
			self.GlacTable.loc[self.GlacTable['GLAC_T'] < 0., 'SnowWatStore_GLAC'] = 0
			#-Changes in total water storage in snow (SnowStore and SnowWatStore)
			self.GlacTable['OldTotalSnowStore_GLAC'] = self.GlacTable['TotalSnowStore_GLAC']
			self.GlacTable['TotalSnowStore_GLAC'] = self.GlacTable['SnowStore_GLAC'] + self.GlacTable['SnowWatStore_GLAC']
			#-Snow runoff
			mask = (self.GlacTable['SnowWatStore_GLAC'] == self.GlacTable['MaxSnowWatStore_GLAC']) & (self.GlacTable['OldSnowStore_GLAC'] > 0.) #-mask with true where SnowWatStore_GLAC == MaxSnowWatStore_GLAC 
			self.GlacTable.loc[mask,'SnowR_GLAC'] = self.GlacTable.loc[mask, 'ActSnowMelt_GLAC'] + self.GlacTable.loc[mask, 'Rain_GLAC'] - \
				(self.GlacTable.loc[mask, 'SnowWatStore_GLAC'] - self.GlacTable.loc[mask, 'OldSnowWatStore_GLAC'])
			self.GlacTable.loc[pcr.numpy.invert(mask), 'SnowR_GLAC'] = 0
			mask = None; del mask

			#-Glacier melt
			self.GlacTable['GlacMelt'] = 0   #-first set to 0 then update hereafter
			#-Masks for full glacier melt (=no snow melt in timestep) and partial glacier melt (=where snowpack is fully melted within timestep)
			partialMelt = (self.GlacTable['OldSnowStore_GLAC'] > 0.) & (self.GlacTable['SnowStore_GLAC'] == 0)  #-mask with true where snowpack has melted within timestep: for these cells also partial glacier melt
			fullMelt = (self.GlacTable['OldSnowStore_GLAC'] == 0) & (self.GlacTable['SnowStore_GLAC'] == 0)
			#-Melt from Clean Ice Glaciers
			mask = (partialMelt & self.CImask)  #-mask for Clean Ice glacier and partial melt
			self.GlacTable.loc[mask, 'GlacMelt'] = pcr.numpy.maximum(self.GlacTable.loc[mask, 'GLAC_T'] - (self.GlacTable.loc[mask, 'OldSnowStore_GLAC'] / self.DDFS), 0) * self.DDFG
			mask = (fullMelt & self.CImask)  #-mask for Clean Ice glacier and full melt
			self.GlacTable.loc[mask, 'GlacMelt'] = Tmelt.loc[mask] * self.DDFG
			#-Melt from Debris Covered Glaciers
			mask = (partialMelt & self.DBmask)  #-mask for Debris covered Glacier and partial melt
			self.GlacTable.loc[mask, 'GlacMelt'] = pcr.numpy.maximum(self.GlacTable.loc[mask, 'GLAC_T'] - (self.GlacTable.loc[mask, 'OldSnowStore_GLAC'] / self.DDFS), 0) * self.DDFDG
			mask = (fullMelt & self.DBmask)  #-mask for Debris covered Glacier and full melt
			self.GlacTable.loc[mask, 'GlacMelt'] = Tmelt.loc[mask] * self.DDFDG
			#-Accumulate glacier melt
			self.GlacTable['AccuGlacMelt'] = self.GlacTable['AccuGlacMelt'] + self.GlacTable['GlacMelt']
			#-Glacier runoff
			mask = self.GlacTable['OldSnowStore_GLAC'] == 0  #-only add rain to glacmelt when there was no snowpack at beginning to time-step
			#-Glacier runoff consisting of melt and rainfall
			self.GlacTable.loc[mask, 'GlacR'] = self.GlacF * (self.GlacTable.loc[mask, 'GlacMelt'] + self.GlacTable.loc[mask, 'Rain_GLAC'])
			#-Glacier percolation consisting of melt and rainfall
			self.GlacTable.loc[mask, 'GlacPerc'] = (1-self.GlacF) * (self.GlacTable.loc[mask, 'GlacMelt'] + self.GlacTable.loc[mask, 'Rain_GLAC'])
			#-Glacier runoff consisting of melt only
			self.GlacTable.loc[pcr.numpy.invert(mask), 'GlacR'] = self.GlacF * self.GlacTable.loc[pcr.numpy.invert(mask), 'GlacMelt']
			#-Glacier percolation consisting of melt only
			self.GlacTable.loc[pcr.numpy.invert(mask), 'GlacPerc'] = (1-self.GlacF) * self.GlacTable.loc[pcr.numpy.invert(mask), 'GlacMelt']
			mask = None; del mask
					
			#-Aggregation for model grid ID 
			GlacTable_MODid = self.GlacTable.loc[:,['Rain_GLAC', 'Snow_GLAC', 'ActSnowMelt_GLAC', 'SnowStore_GLAC',\
									'SnowWatStore_GLAC', 'TotalSnowStore_GLAC', 'SnowR_GLAC', 'GlacMelt', 'GlacR', 'GlacPerc']]
			GlacTable_MODid = GlacTable_MODid.multiply(self.GlacTable['FRAC_GLAC'], axis='index') #-Multiply with the glacier fraction
			GlacTable_MODid = GlacTable_MODid.groupby(GlacTable_MODid.index).sum() #-Summarize by model ID
			#-report back to model ID
			#-Rainfall on glacier
			Rain_GLAC = pcr.numpy.zeros(self.ModelID_1d.shape)
			Rain_GLAC[self.GlacierKeys] = GlacTable_MODid['Rain_GLAC']
			Rain_GLAC = Rain_GLAC.reshape(self.ModelID.shape)
			Rain_GLAC = pcr.numpy2pcr(Scalar, Rain_GLAC, self.MV)
			#-Snowfall on glacier
			Snow_GLAC = pcr.numpy.zeros(self.ModelID_1d.shape)
			Snow_GLAC[self.GlacierKeys] = GlacTable_MODid['Snow_GLAC']
			Snow_GLAC = Snow_GLAC.reshape(self.ModelID.shape)
			Snow_GLAC = pcr.numpy2pcr(Scalar, Snow_GLAC, self.MV)
			#-Act snowmelt from glacier
			ActSnowMelt_GLAC = pcr.numpy.zeros(self.ModelID_1d.shape)
			ActSnowMelt_GLAC[self.GlacierKeys] = GlacTable_MODid['ActSnowMelt_GLAC']
			ActSnowMelt_GLAC = ActSnowMelt_GLAC.reshape(self.ModelID.shape)
			ActSnowMelt_GLAC = pcr.numpy2pcr(Scalar, ActSnowMelt_GLAC, self.MV)
# 			#-Snowstore on glacier
# 			SnowStore_GLAC = pcr.numpy.zeros(self.ModelID_1d.shape)
# 			SnowStore_GLAC[self.GlacierKeys] = GlacTable_MODid['SnowStore_GLAC']
# 			SnowStore_GLAC = SnowStore_GLAC.reshape(self.ModelID.shape)
# 			SnowStore_GLAC = pcr.numpy2pcr(Scalar, SnowStore_GLAC, self.MV)
# 			#-SnowWatStore on glacier
# 			SnowWatStore_GLAC = pcr.numpy.zeros(self.ModelID_1d.shape)
# 			SnowWatStore_GLAC[self.GlacierKeys] = GlacTable_MODid['SnowWatStore_GLAC']
# 			SnowWatStore_GLAC = SnowWatStore_GLAC.reshape(self.ModelID.shape)
# 			SnowWatStore_GLAC = pcr.numpy2pcr(Scalar, SnowWatStore_GLAC, self.MV)
			#-TotalSnowStore on glacier
			TotalSnowStore_GLAC = pcr.numpy.zeros(self.ModelID_1d.shape)
			TotalSnowStore_GLAC[self.GlacierKeys] = GlacTable_MODid['TotalSnowStore_GLAC']
			TotalSnowStore_GLAC = TotalSnowStore_GLAC.reshape(self.ModelID.shape)
			TotalSnowStore_GLAC = pcr.numpy2pcr(Scalar, TotalSnowStore_GLAC, self.MV)
			#-SnowR from glacier
			SnowR_GLAC = pcr.numpy.zeros(self.ModelID_1d.shape)
			SnowR_GLAC[self.GlacierKeys] = GlacTable_MODid['SnowR_GLAC']
			SnowR_GLAC = SnowR_GLAC.reshape(self.ModelID.shape)
			SnowR_GLAC = pcr.numpy2pcr(Scalar, SnowR_GLAC, self.MV)
			#-Glacier melt
			GlacMelt = pcr.numpy.zeros(self.ModelID_1d.shape)
			GlacMelt[self.GlacierKeys] = GlacTable_MODid['GlacMelt']
			GlacMelt = GlacMelt.reshape(self.ModelID.shape)
			GlacMelt = pcr.numpy2pcr(Scalar, GlacMelt, self.MV)
			#-Report glacier melt
			self.reporting.reporting(self, pcr, 'TotGlacMelt', GlacMelt)
			if self.mm_rep_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
				self.GMeltSubBasinTSS.sample(pcr.catchmenttotal(GlacMelt, self.FlowDir) / pcr.catchmenttotal(1, self.FlowDir))
			#-Glacier runoff
			GlacR = pcr.numpy.zeros(self.ModelID_1d.shape)
			GlacR[self.GlacierKeys] = GlacTable_MODid['GlacR']
			GlacR = GlacR.reshape(self.ModelID.shape)
			GlacR = pcr.numpy2pcr(Scalar, GlacR, self.MV)
			#-Report glacier runoff
			self.reporting.reporting(self, pcr, 'TotGlacR', GlacR)
			#-Glacier percolation 
			GlacPerc = pcr.numpy.zeros(self.ModelID_1d.shape)
			GlacPerc[self.GlacierKeys] = GlacTable_MODid['GlacPerc']
			GlacPerc = GlacPerc.reshape(self.ModelID.shape)
			GlacPerc = pcr.numpy2pcr(Scalar, GlacPerc, self.MV)
			#-Report glacier percolation to groundwater
			self.reporting.reporting(self, pcr, 'TotGlacPerc', GlacPerc)
			GlacTable_MODid = None; del GlacTable_MODid
		#-If glacier module is not used, then
		else:
			Rain_GLAC = 0
			Snow_GLAC = 0
			ActSnowMelt_GLAC = 0
			TotalSnowStore_GLAC = 0
			SnowR_GLAC = 0
			GlacR = 0
			GlacMelt = 0
			GlacPerc = 0

		# Calculate snow and rain for non-glacier part of cell
		if self.SnowFLAG == 1:
			#-Snow and rain differentiation
			Snow = pcr.ifthenelse(Temp >= self.Tcrit, 0, Precip)
			Rain = pcr.ifthenelse(Temp < self.Tcrit, 0, Precip)
			#-Report Snow for entire cell (snow+glacier fraction)
			self.reporting.reporting(self, pcr, 'TotSnow', Snow * (1-self.GlacFrac) + Snow_GLAC)
			#-Snow melt
			PotSnowMelt = self.snow.PotSnowMelt(pcr, Temp, self.DDFS)
			ActSnowMelt = self.snow.ActSnowMelt(pcr, self.SnowStore, PotSnowMelt)
			#-Report snow melt for entire cell (snow+glacier fraction)
			self.reporting.reporting(self, pcr, 'TotSnowMelt', ActSnowMelt * (1-self.GlacFrac) + ActSnowMelt_GLAC)
			#-Update snow store
			self.SnowStore = self.snow.SnowStoreUpdate(pcr, self.SnowStore, Snow, ActSnowMelt, Temp, self.SnowWatStore)
			#-Caclulate the maximum amount of water that can be stored in snowwatstore
			MaxSnowWatStore = self.snow.MaxSnowWatStorage(self.SnowSC, self.SnowStore)
			OldSnowWatStore = self.SnowWatStore
			#-Calculate the actual amount of water stored in snowwatstore
			self.SnowWatStore = self.snow.SnowWatStorage(pcr, Temp, MaxSnowWatStore, self.SnowWatStore, ActSnowMelt, Rain)
			#-Changes in total water storage in snow (SnowStore and SnowWatStore)
			OldTotalSnowStore = self.TotalSnowStore
			self.TotalSnowStore = self.snow.TotSnowStorage(self.SnowStore, self.SnowWatStore, SnowFrac, RainFrac) + TotalSnowStore_GLAC  # for entire cell
			#-Snow runoff
			SnowR = self.snow.SnowR(pcr, self.SnowWatStore, MaxSnowWatStore, ActSnowMelt, Rain, OldSnowWatStore, SnowFrac) + SnowR_GLAC  # for entire cell
			#-Report Snow runoff
			self.reporting.reporting(self, pcr, 'TotSnowR', SnowR)
		else:
			Rain = Precip
			SnowR = 0
			OldTotalSnowStore = 0
			self.TotalSnowStore = 0
		#-Report Rain
		self.reporting.reporting(self, pcr, 'TotRain', Rain * (1-self.GlacFrac) + Rain_GLAC)  # for entire cell

		#-Potential evapotranspiration
		ETpot = self.ET.ETpot(ETref, self.Kc) 
		#-Report ETpot
		self.reporting.reporting(self, pcr, 'TotETpot', ETpot * RainFrac)
				
		#-Rootzone calculations
		self.RootWater = self.RootWater + pcr.ifthenelse(RainFrac > 0, Rain, 0) + self.CapRise
		#-Rootzone runoff
		RootRunoff = self.rootzone.RootRunoff(pcr, RainFrac, self.RootWater, self.RootSat)
		self.RootWater = self.RootWater - RootRunoff
		#-Actual evapotranspiration
		etreddry = pcr.max(pcr.min((self.RootWater - self.RootDry) / (self.RootWilt - self.RootDry), 1), 0)
		ETact = self.ET.ETact(pcr, ETpot, self.RootWater, self.RootSat, etreddry, RainFrac)
		#-Actual evapotranspiration, corrected for rain fraction
		ActETact = ETact * RainFrac	
		#-Report the actual evapotranspiration
		self.reporting.reporting(self, pcr, 'TotETact', ActETact)
		if self.mm_rep_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
			self.ETaSubBasinTSS.sample(pcr.catchmenttotal(ActETact, self.FlowDir) / pcr.catchmenttotal(1, self.FlowDir))
		#-Update rootwater content
		self.RootWater = pcr.max(self.RootWater - ETact, 0)
		#-Rootwater drainage
		self.RootDrain = self.rootzone.RootDrainage(pcr, self.RootWater, self.RootDrain, self.RootField, self.RootSat, self.RootDrainVel, self.RootTT)
		#-Update rootwater content
		self.RootWater = self.RootWater - self.RootDrain
		#-Rootwater percolation
		rootperc = self.rootzone.RootPercolation(pcr, self.RootWater, self.SubWater, self.RootField, self.RootTT, self.SubSat)
		#-Report rootzone percolation
		self.reporting.reporting(self, pcr, 'TotRootP', rootperc * (1 - self.GlacFrac))
		#-Update rootwater content
		self.RootWater = self.RootWater - rootperc
		
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
			self.reporting.reporting(self, pcr, 'TotSeep', pcr.scalar(self.SeePage))
			self.SubWater = pcr.min(pcr.max(self.SubWater - self.SeePage, 0), self.SubSat)
			if self.mm_rep_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
				self.SeepSubBasinTSS.sample(pcr.catchmenttotal(self.SeePage, self.FlowDir) / pcr.catchmenttotal(1, self.FlowDir))
		#-Capillary rise
		self.CapRise = self.subzone.CapilRise(pcr, self.SubField, self.SubWater, self.CapRiseMax, self.RootWater, self.RootSat, self.RootField)
		#-Report capillary rise
		self.reporting.reporting(self, pcr, 'TotCapR', self.CapRise * (1-self.GlacFrac))
		#-Update sub soil water content
		self.SubWater = self.SubWater - self.CapRise
		if self.GroundFLAG == 1:   # sub percolation will be calculated instead of subdrainage
			subperc = self.subzone.SubPercolation(pcr, self.SubWater, self.SubField, self.SubTT, self.Gw, self.GwSat)
			ActSubPerc = subperc * (1-self.GlacFrac)
			#-Report the subzone percolation
			self.reporting.reporting(self, pcr, 'TotSubP', ActSubPerc) 
			#-Update sub soil water content
			self.SubWater = self.SubWater - subperc
		else: # sub drainage will be calculated instead of sub percolation
			self.SubDrain = self.subzone.SubDrainage(pcr, self.SubWater, self.SubField, self.SubSat, self.SubDrainVel, self.SubDrain, self.SubTT)
			#-Report drainage from subzone
			self.reporting.reporting(self, pcr, 'TotSubD', self.SubDrain)
			#-Update sub soil water content
			self.SubWater = self.SubWater - self.SubDrain
			
		#-Rootzone runoff
		RootR = RootRunoff * RainFrac
		#-Report rootzone runoff
		self.reporting.reporting(self, pcr, 'TotRootR', RootR)
		#-Rootzone drainage
		RootD = self.RootDrain * (1-self.GlacFrac)
		#-Report rootzone drainage
		self.reporting.reporting(self, pcr, 'TotRootD', RootD)
		#-Rain runoff
		RainR = RootR + RootD
		#-Report rain runoff
		self.reporting.reporting(self, pcr, 'TotRainR', RainR)
		
		#-Groundwater calculations
		if self.GroundFLAG == 1:
			#-Groundwater recharge
			self.GwRecharge = self.groundwater.GroundWaterRecharge(pcr,	self.deltaGw, self.GwRecharge, ActSubPerc, GlacPerc)
			#-Report groundwater recharge
			self.reporting.reporting(self, pcr, 'TotGwRecharge', self.GwRecharge)
			#-Update groundwater storage
			self.Gw = self.Gw + self.GwRecharge
			#-Baseflow
			self.BaseR = self.groundwater.BaseFlow(pcr, self.Gw, self.BaseR, self.GwRecharge, self.BaseThresh, self.alphaGw)
			#-Report Baseflow
			self.reporting.reporting(self, pcr, 'TotBaseR', self.BaseR)
			#-Update groundwater storage
			self.Gw = self.Gw - self.BaseR
			#-Calculate groundwater level
			self.H_gw = self.groundwater.HLevel(pcr, self.H_gw, self.alphaGw, self.GwRecharge, self.YieldGw)
			#-Report groundwater
			self.reporting.reporting(self, pcr, 'GWL', ((self.SubDepthFlat + self.RootDepthFlat + self.GwDepth)/1000 - self.H_gw)*-1)
		else:
			#-Use drainage from subsoil as baseflow
			self.BaseR = self.SubDrain
			#-Groundwater level as scaled between min and max measured gwl
			SoilAct = self.RootWater + self.SubWater;
			SoilRel = (SoilAct - self.SoilMin) / (self.SoilMax - self.SoilMin) # scale between 0 (dry) and 1 (wet)
			GWL = self.GWL_base - (SoilRel-0.5) * self.GWL_base
			#-Report groundwater
			self.reporting.reporting(self, pcr, 'GWL', GWL)
			
		#-Report Total runoff
		self.reporting.reporting(self, pcr, 'TotR', self.BaseR + RainR + SnowR + GlacR)

		#-Routing for lake and/or reservoir modules
		if self.LakeFLAG == 1 or self.ResFLAG == 1:
			#-Update storage in lakes/reservoirs (m3) with specific runoff
			self.StorRES = self.StorRES + pcr.ifthenelse(self.QFRAC==0, 0.001 * pcr.cellarea() * (self.BaseR + RainR + GlacR + SnowR), 0)
			OldStorage = self.StorRES
			#-Calculate lake/reservoir outflow volumes
			if self.LakeFLAG ==1 and self.ResFLAG ==1:
				tempvar = self.lakes.UpdateLakeHStore(self, pcr, pcrm)
				LakeLevel = tempvar[0]
				self.StorRES = tempvar[1]
				LakeQ = self.lakes.QLake(self, pcr, LakeLevel)
				ResQ = self.reservoirs.QRes(self, pcr)
				Qout = pcr.ifthenelse(self.ResID != 0, ResQ, pcr.ifthenelse(self.LakeID!=0, LakeQ, 0))
			elif self.LakeFLAG ==1:
				tempvar = self.lakes.UpdateLakeHStore(self, pcr, pcrm)
				LakeLevel = tempvar[0]
				self.StorRES = tempvar[1]
				Qout = self.lakes.QLake(self, pcr, LakeLevel)
			else:
				Qout = self.reservoirs.QRes(self, pcr)

			#-Calculate volume available for routing (=outflow lakes/reservoir + cell specific runoff)
			RunoffVolume = pcr.upstream(self.FlowDir, Qout) + pcr.ifthenelse(self.QFRAC==0, 0, 0.001 * pcr.cellarea() * (self.BaseR + RainR + GlacR + SnowR))
			#-Routing of total flow
			tempvar = self.routing.ROUT(self, pcr, RunoffVolume, self.QRAold, Qout, self.StorRES)
			self.StorRES = tempvar[0]
			Q = tempvar[1]
			Qin = tempvar[2]
			self.QRAold = Q
			self.reporting.reporting(self, pcr, 'QallRAtot', Q)
			#-report flux in mm
			if self.mm_rep_FLAG == 1:
				self.QTOTSubBasinTSS.sample(((Q * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) * 1000)
			#-report lake and reservoir waterbalance
			if self.LakeFLAG == 1 and config.getint('REPORTING', 'Lake_wbal') ==1:
				self.LakeInTSS.sample(Qin)
				self.LakeOutTSS.sample(Qout) 
				self.LakeStorTSS.sample(self.StorRES)
			if self.ResFLAG == 1 and config.getint('REPORTING', 'Res_wbal') ==1:
				self.ResInTSS.sample(Qin)
				self.ResOutTSS.sample(Qout) 
				self.ResStorTSS.sample(self.StorRES)

			#-Routing of individual contributers
			#-Snow routing
			if self.SnowRA_FLAG == 1 and self.SnowFLAG == 1:
				self.SnowRAstor = self.SnowRAstor + pcr.ifthenelse(self.QFRAC==0, SnowR * 0.001 * pcr.cellarea(), 0)
				cQfrac = pcr.cover(self.SnowRAstor / OldStorage, 0)
				cQout = cQfrac * Qout
				cRunoffVolume = pcr.upstream(self.FlowDir, cQout) + pcr.ifthenelse(self.QFRAC==0, 0, 0.001 * pcr.cellarea() * SnowR)
				tempvar = self.routing.ROUT(self, pcr, cRunoffVolume, self.SnowRAold, cQout, self.SnowRAstor)
				self.SnowRAstor = tempvar[0]
				SnowRA = tempvar[1]
				cQin = tempvar[2]
				self.SnowRAold = SnowRA
				self.reporting.reporting(self, pcr, 'SnowRAtot', SnowRA)
				if self.mm_rep_FLAG == 1:
					self.QSNOWSubBasinTSS.sample(((SnowRA * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) *1000)
				#-report lake and reservoir waterbalance
				if self.LakeFLAG == 1 and config.getint('REPORTING', 'Lake_wbal') ==1:
					self.LakeSnowInTSS.sample(cQin)
					self.LakeSnowOutTSS.sample(cQout) 
					self.LakeSnowStorTSS.sample(self.SnowRAstor)
				if self.ResFLAG == 1 and config.getint('REPORTING', 'Res_wbal') ==1:
					self.ResSnowInTSS.sample(cQin)
					self.ResSnowOutTSS.sample(cQout) 
					self.ResSnowStorTSS.sample(self.SnowRAstor)
			#-Rain routing
			if self.RainRA_FLAG == 1:
				self.RainRAstor = self.RainRAstor + pcr.ifthenelse(self.QFRAC==0, RainR * 0.001 * pcr.cellarea(), 0)
				cQfrac = pcr.cover(self.RainRAstor / OldStorage, 0)
				cQout = cQfrac * Qout
				cRunoffVolume = pcr.upstream(self.FlowDir, cQout) + pcr.ifthenelse(self.QFRAC==0, 0, 0.001 * pcr.cellarea() * RainR)
				tempvar = self.routing.ROUT(self, pcr, cRunoffVolume, self.RainRAold, cQout, self.RainRAstor)
				self.RainRAstor = tempvar[0]
				RainRA = tempvar[1]
				cQin = tempvar[2]
				self.RainRAold = RainRA
				self.reporting.reporting(self, pcr, 'RainRAtot', RainRA)
				if self.mm_rep_FLAG == 1:
					self.QRAINSubBasinTSS.sample(((RainRA * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) *1000)
				#-report lake and reservoir waterbalance
				if self.LakeFLAG == 1 and config.getint('REPORTING', 'Lake_wbal') ==1:
					self.LakeRainInTSS.sample(cQin)
					self.LakeRainOutTSS.sample(cQout) 
					self.LakeRainStorTSS.sample(self.RainRAstor)
				if self.ResFLAG == 1 and config.getint('REPORTING', 'Res_wbal') ==1:
					self.ResRainInTSS.sample(cQin)
					self.ResRainOutTSS.sample(cQout) 
					self.ResRainStorTSS.sample(self.RainRAstor)
			#-Glacier routing
			if self.GlacRA_FLAG == 1 and self.GlacFLAG == 1:
				self.GlacRAstor = self.GlacRAstor + pcr.ifthenelse(self.QFRAC==0, GlacR * 0.001 * pcr.cellarea(), 0)
				cQfrac = pcr.cover(self.GlacRAstor / OldStorage, 0)
				cQout = cQfrac * Qout
				cRunoffVolume = pcr.upstream(self.FlowDir, cQout) + pcr.ifthenelse(self.QFRAC==0, 0, 0.001 * pcr.cellarea() * GlacR)
				tempvar = self.routing.ROUT(self, pcr, cRunoffVolume, self.GlacRAold, cQout, self.GlacRAstor)
				self.GlacRAstor = tempvar[0]
				GlacRA = tempvar[1]
				cQin = tempvar[2]
				self.GlacRAold = GlacRA
				self.reporting.reporting(self, pcr, 'GlacRAtot', GlacRA)
				if self.mm_rep_FLAG == 1:
					self.QGLACSubBasinTSS.sample(((GlacRA * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) *1000)
				#-report lake and reservoir waterbalance
				if self.LakeFLAG == 1 and config.getint('REPORTING', 'Lake_wbal') ==1:
					self.LakeGlacInTSS.sample(cQin)
					self.LakeGlacOutTSS.sample(cQout) 
					self.LakeGlacStorTSS.sample(self.GlacRAstor)
				if self.ResFLAG == 1 and config.getint('REPORTING', 'Res_wbal') ==1:
					self.ResGlacInTSS.sample(cQin)
					self.ResGlacOutTSS.sample(cQout) 
					self.ResGlacStorTSS.sample(self.GlacRAstor)
			#-Baseflow routing
			if self.BaseRA_FLAG == 1:
				self.BaseRAstor = self.BaseRAstor + pcr.ifthenelse(self.QFRAC==0, self.BaseR * 0.001 * pcr.cellarea(), 0)
				cQfrac = pcr.cover(self.BaseRAstor / OldStorage, 0)
				cQout = cQfrac * Qout
				cRunoffVolume = pcr.upstream(self.FlowDir, cQout) + pcr.ifthenelse(self.QFRAC==0, 0, 0.001 * pcr.cellarea() * self.BaseR)
				tempvar = self.routing.ROUT(self, pcr, cRunoffVolume, self.BaseRAold, cQout, self.BaseRAstor)
				self.BaseRAstor = tempvar[0]
				BaseRA = tempvar[1]
				cQin = tempvar[2]
				self.BaseRAold = BaseRA
				self.reporting.reporting(self, pcr, 'BaseRAtot', BaseRA)
				if self.mm_rep_FLAG == 1:
					self.QBASFSubBasinTSS.sample(((BaseRA * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) *1000)
				#-report lake and reservoir waterbalance
				if self.LakeFLAG == 1 and config.getint('REPORTING', 'Lake_wbal') ==1:
					self.LakeBaseInTSS.sample(cQin)
					self.LakeBaseOutTSS.sample(cQout) 
					self.LakeBaseStorTSS.sample(self.BaseRAstor)
				if self.ResFLAG == 1 and config.getint('REPORTING', 'Res_wbal') ==1:
					self.ResBaseInTSS.sample(cQin)
					self.ResBaseOutTSS.sample(cQout) 
					self.ResBaseStorTSS.sample(self.BaseRAstor)
						
		#-Normal routing module	
		elif self.RoutFLAG == 1:
			#-Rout total runoff
			Q = self.routing.ROUT(pcr, self.BaseR + RainR + GlacR + SnowR, self.QRAold, self.FlowDir, self.kx)
			self.QRAold = Q
			self.reporting.reporting(self, pcr, 'QallRAtot', Q)
			if self.mm_rep_FLAG == 1:
				self.QTOTSubBasinTSS.sample(((Q * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) * 1000)
			#-Snow routing
			if self.SnowRA_FLAG == 1 and self.SnowFLAG == 1:
				SnowRA = self.routing.ROUT(pcr, SnowR, self.SnowRAold, self.FlowDir, self.kx)
				self.SnowRAold = SnowRA
				self.reporting.reporting(self, pcr, 'SnowRAtot', SnowRA)
				if self.mm_rep_FLAG == 1:
					self.QSNOWSubBasinTSS.sample(((SnowRA * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) *1000)
			#-Rain routing
			if self.RainRA_FLAG == 1:
				RainRA = self.routing.ROUT(pcr, RainR, self.RainRAold, self.FlowDir, self.kx)
				self.RainRAold = RainRA
				self.reporting.reporting(self, pcr, 'RainRAtot', RainRA)
				if self.mm_rep_FLAG == 1:
					self.QRAINSubBasinTSS.sample(((RainRA * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) *1000)
			#-Glacier routing
			if self.GlacRA_FLAG == 1 and self.GlacFLAG == 1:
				GlacRA = self.routing.ROUT(pcr, GlacR, self.GlacRAold, self.FlowDir, self.kx)
				self.GlacRAold = GlacRA
				self.reporting.reporting(self, pcr, 'GlacRAtot', GlacRA)
				if self.mm_rep_FLAG == 1:
					self.QGLACSubBasinTSS.sample(((GlacRA * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) *1000)
			#-Baseflow routing
			if self.BaseRA_FLAG == 1:
				BaseRA = self.routing.ROUT(pcr, self.BaseR, self.BaseRAold, self.FlowDir, self.kx)
				self.BaseRAold = BaseRA
				self.reporting.reporting(self, pcr, 'BaseRAtot', BaseRA)
				if self.mm_rep_FLAG == 1:
					self.QBASFSubBasinTSS.sample(((BaseRA * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) *1000)
					
			if self.GlacFLAG:
				#-Reporting for each glacier ID
				if self.GlacID_flag:
					GlacTable_GLACid = self.GlacTable.loc[:, self.GlacVars]

					#-Muliply with fraction, summarize, and divide by total fraction to get glacier ID average 
					GlacTable_GLACid = GlacTable_GLACid.multiply(self.GlacTable['FRAC_GLAC'], axis='index')  #-Multiply with fraction
					GlacTable_GLACid['GLAC_ID'] = self.GlacTable['GLAC_ID']; GlacTable_GLACid['FRAC_GLAC'] = self.GlacTable['FRAC_GLAC']  #-Add GLAC_ID column and FRAC_GLAC column
					GlacTable_GLACid = GlacTable_GLACid.groupby('GLAC_ID').sum() #-Summarize by Glacier ID
					FracSum = GlacTable_GLACid['FRAC_GLAC'] #-Get summed glacier fraction					
					GlacTable_GLACid = GlacTable_GLACid.div(FracSum, axis='index') #-Divide by total fraction for glacier weighted average (frac==1)
					GlacTable_GLACid['FRAC_GLAC'] = FracSum; FracSum = None; del FracSum  
					GlacTable_GLACid = GlacTable_GLACid.transpose()  #-Transpose the glacier id table (ID as columns and vars as index)
					GlacTable_GLACid.fillna(0., inplace=True);
					
					#-Fill Glacier variable tables for reporting
					if self.GlacID_memerror == 0:
						for v in self.GlacVars:
							vv = getattr(self, v + '_Table'); vv.loc[self.curdate,:] = GlacTable_GLACid.loc[v,:]
						v = None; vv = None; del v, vv; GlacTable_GLACid = None; del GlacTable_GLACid
						if self.curdate == self.enddate: #-do the reporting at the final model time-step
							for v in self.GlacVars:
								eval('self.' + v + '_Table.to_csv("'  + self.outpath + v + '.csv")')
					else:
						df = pd.DataFrame(columns=self.glacid, dtype=np.float32)
						for v in self.GlacVars:
							#df = pd.DataFrame(columns=self.glacid, dtype=float)
							df.loc[self.curdate,:] = GlacTable_GLACid.loc[v,:]
							if self.curdate == self.startdate:
								df.to_csv(self.outpath + v + '.csv', mode='w')
							else:
								df.to_csv(self.outpath + v + '.csv', mode='a', header=False) #-no header for time-step > 1
						df = None; del df; GlacTable_GLACid = None; del GlacTable_GLACid

				#-Check if glacier retreat should be calculated
				if self.GlacRetreat == 1 and self.curdate.month == self.GlacUpdate['month'] and self.curdate.day == self.GlacUpdate['day']:
					self.dateAfterUpdate = self.curdate + self.datetime.timedelta(days=1)
					#-Create a table with fields used for updating the glacier fraction at defined update date
					GlacFracTable = self.GlacTable.loc[:,['U_ID', 'GLAC_ID','FRAC_GLAC','ICE_DEPTH']]
					#-Set the initial ice volumes and snow store
					GlacFracTable['V_ice_t0'] = GlacFracTable['FRAC_GLAC'] * GlacFracTable['ICE_DEPTH'] * self.cellArea
					GlacFracTable['TotalSnowStore_GLAC'] = self.GlacTable['TotalSnowStore_GLAC']
					GlacFracTable['AccuGlacMelt'] = self.GlacTable['AccuGlacMelt']
					GlacFracTable['dMelt'] = GlacFracTable['TotalSnowStore_GLAC'] - GlacFracTable['AccuGlacMelt']
					GlacFracTable['dMelt'] = GlacFracTable['dMelt'] / 1000 * GlacFracTable['FRAC_GLAC'] * self.cellArea  #-convert to m3
					#-Drop unnecessary columns
					GlacFracTable.drop(['TotalSnowStore_GLAC', 'AccuGlacMelt'], axis=1, inplace=True)
					
					#-Mask to determine ablation and accumulation UIDs
					ablMask = GlacFracTable['dMelt'] < 0.
					#-Set the ablation and accumulation in the corresponding fields
					GlacFracTable['Accumulation'] = 0.
					GlacFracTable['Ablation'] = 0.
					GlacFracTable.loc[pcr.numpy.invert(ablMask),'Accumulation'] = GlacFracTable.loc[pcr.numpy.invert(ablMask), 'dMelt']
					GlacFracTable.loc[ablMask,'Ablation'] = GlacFracTable.loc[ablMask, 'dMelt']
					#-Set the ice volumes for the ablation UIDs
					GlacFracTable['V_ice_ablation'] = 0.
					GlacFracTable.loc[ablMask,'V_ice_ablation'] = GlacFracTable.loc[ablMask,'V_ice_t0']
					#-Calculate totals per Glacier ID
					GlacID_grouped = GlacFracTable.groupby('GLAC_ID').sum()
					GlacID_grouped = GlacID_grouped.loc[:,['dMelt', 'Accumulation', 'V_ice_ablation']]

					#-Calculate total delta Melt (dMelt), accumulation, and ice volumes (of ablation cells) for each glacier ID
					GlacFracTable['dMelt_group'] = 0.
					GlacFracTable['Accumulation_group'] = 0.
					GlacFracTable['V_ice_ablation_group'] = 0.

					for index, row in GlacID_grouped.iterrows():
						GlacFracTable.loc[GlacFracTable['GLAC_ID'] == index,'dMelt_group'] = row['dMelt']
						GlacFracTable.loc[GlacFracTable['GLAC_ID'] == index,'Accumulation_group'] = row['Accumulation']
						GlacFracTable.loc[GlacFracTable['GLAC_ID'] == index,'V_ice_ablation_group'] = row['V_ice_ablation']
					#-Remove GlacID_grouped table
					GlacID_grouped = None; del GlacID_grouped
					#-Mask for determining if redistribution is negative (remove ice from accumulation cells) or positive (add ice to ablation cells)
					negDistMask = (GlacFracTable['dMelt_group'] < 0.) & (GlacFracTable['dMelt'] >= 0.)
					posDistMask = (GlacFracTable['dMelt_group'] < 0.) & (GlacFracTable['dMelt'] < 0.)
					#-Calculate the ice redistribution
					GlacFracTable['Ice_redist'] = 0.
					GlacFracTable.loc[negDistMask,'Ice_redist'] = -GlacFracTable.loc[negDistMask, 'Accumulation']
					GlacFracTable.loc[posDistMask,'Ice_redist'] = GlacFracTable.loc[posDistMask,'V_ice_ablation'] / GlacFracTable.loc[posDistMask,'V_ice_ablation_group'] *\
						GlacFracTable.loc[posDistMask,'Accumulation_group']
					#-Update ice volume
					GlacFracTable['V_ice_t1'] = GlacFracTable['V_ice_t0'] + GlacFracTable['Accumulation'] + GlacFracTable['Ablation'] + GlacFracTable['Ice_redist']
					#-Remove distribution masks
					negDistMask = None; del negDistMask; posDistMask = None; del posDistMask
					#-Remove unnecessary columns
					GlacFracTable.drop(['dMelt', 'Accumulation', 'Ablation','V_ice_ablation','dMelt_group',\
							'Accumulation_group', 'V_ice_ablation_group', 'Ice_redist'], axis=1, inplace=True)
					#-Calculate where updated ice volume becomes negative and postitive
					GlacFracTable['V_ice_negative'] = pcr.numpy.minimum(0., GlacFracTable['V_ice_t1'])
					GlacFracTable['V_ice_positive'] = pcr.numpy.maximum(0., GlacFracTable['V_ice_t1'])
					#-Calculate totals per Glacier ID
					GlacID_grouped = GlacFracTable.groupby('GLAC_ID').sum()
					GlacID_grouped = GlacID_grouped.loc[:,['V_ice_negative', 'V_ice_positive']]
					GlacFracTable['V_ice_negative_group'] = 0.
					GlacFracTable['V_ice_positive_group'] = 0.
					for index, row in GlacID_grouped.iterrows():
						GlacFracTable.loc[GlacFracTable['GLAC_ID'] == index,'V_ice_negative_group'] = row['V_ice_negative']
						GlacFracTable.loc[GlacFracTable['GLAC_ID'] == index,'V_ice_positive_group'] = row['V_ice_positive']
					#-Remove GlacID_grouped table
					GlacID_grouped = None; del GlacID_grouped
					#-Calculate the ice redistribution
					GlacFracTable['Ice_redist'] = GlacFracTable['V_ice_positive'] / GlacFracTable['V_ice_positive_group'] * GlacFracTable['V_ice_negative_group']
					GlacFracTable['Ice_redist'].fillna(0., inplace=True)
					#-Remove unnecessary columns
					GlacFracTable.drop(['V_ice_negative', 'V_ice_positive', 'V_ice_negative_group', 'V_ice_positive_group'], axis=1, inplace=True)
					#-Update ice volume
					GlacFracTable['V_ice_t2'] = pcr.numpy.maximum(0., GlacFracTable['V_ice_t1'] + GlacFracTable['Ice_redist'])
					#-Update ice thickness
					GlacFracTable['ICE_DEPTH_new'] = GlacFracTable['V_ice_t2'] / (GlacFracTable['FRAC_GLAC'] * self.cellArea)
					noIceMask = (GlacFracTable['ICE_DEPTH_new'] <=0.)  # it can be that melt is greater than total available ice: that results in balance error
					#-Update glacier fraction
					GlacFracTable['FRAC_GLAC_new'] = GlacFracTable['FRAC_GLAC']
					GlacFracTable.loc[noIceMask,'FRAC_GLAC_new'] = 0.
					noIceMask = None; del noIceMask

					#-Update the glactable
					self.GlacTable['FRAC_GLAC'] = GlacFracTable['FRAC_GLAC_new']
					self.GlacTable['ICE_DEPTH'] = GlacFracTable['ICE_DEPTH_new']
					self.GlacTable['SnowStore_GLAC'] = 0.
					self.GlacTable['SnowWatStore_GLAC'] = 0.
					self.GlacTable['TotalSnowStore_GLAC'] = 0.
					
					#-remove SnowWatStore_GLAC from total snowstore
					self.TotalSnowStore = self.TotalSnowStore - TotalSnowStore_GLAC
					
					#-Set accumulated glacier melt to zero as initial condition for next period
					self.GlacTable['AccuGlacMelt'] = 0.
					#-Remove the GlacFracTable
					GlacFracTable = None; del GlacFracTable
					
				#-Write glacier table, spatial maps of glacier fraction, ice depth, and average ice depth per glacier on day after update day
				if self.GlacRetreat == 1 and self.curdate == self.dateAfterUpdate:
					#-Write glac table to csv (can be used as initial setting for new run)
					glac_csv = self.GlacTable.loc[:,['U_ID','GLAC_ID','MOD_H','GLAC_H','DEBRIS','FRAC_GLAC','ICE_DEPTH']]
					glac_csv.insert(1, 'MOD_ID', self.GlacTable.index)
					glac_csv.to_csv(self.outpath + 'glacTable_' + self.curdate.strftime('%Y%m%d') + '.csv', index=False)
					glac_csv = None; del glac_csv
										
					#-Calculate average model ID ice depth and total model glacier fraction
					GlacTable_MODid = self.GlacTable.loc[:,['FRAC_GLAC', 'ICE_DEPTH']]
					GlacTable_MODid['ICE_DEPTH'] = GlacTable_MODid['ICE_DEPTH'] * GlacTable_MODid['FRAC_GLAC']
					GlacTable_MODid['AREA'] = GlacTable_MODid['FRAC_GLAC'] * self.cellArea
					GlacTable_MODid = GlacTable_MODid.groupby(GlacTable_MODid.index).sum()
					GlacTable_MODid.fillna(0., inplace=True)
					
					#-Report updated glacier fraction map
					self.GlacFrac = pcr.numpy.zeros(self.ModelID_1d.shape)
					self.GlacFrac[self.GlacierKeys] = GlacTable_MODid['FRAC_GLAC']
					self.GlacFrac = self.GlacFrac.reshape(self.ModelID.shape)
					self.GlacFrac = pcr.numpy2pcr(Scalar, self.GlacFrac, self.MV)
					self.GlacFrac = pcr.ifthen(self.clone, self.GlacFrac)  #-only use values where clone is True
					pcr.report(self.GlacFrac, self.outpath + 'GlacFrac_' + self.curdate.strftime('%Y%m%d') + '.map')
					
					#-Report pcraster map of glacier area
					glacArea = pcr.numpy.zeros(self.ModelID_1d.shape)
					glacArea[self.GlacierKeys] = GlacTable_MODid['AREA']
					glacArea = glacArea.reshape(self.ModelID.shape)
					glacArea = pcr.numpy2pcr(Scalar, glacArea, self.MV)
					glacArea = pcr.ifthen(self.clone, glacArea)  #-only use values where clone is True
					pcr.report(glacArea, self.outpath + 'GlacArea_' + self.curdate.strftime('%Y%m%d') + '.map')
					#-Report pcraster map of glacier depth
					iceDepth = pcr.numpy.zeros(self.ModelID_1d.shape)
					iceDepth[self.GlacierKeys] = GlacTable_MODid['ICE_DEPTH']
					iceDepth = iceDepth.reshape(self.ModelID.shape)
					iceDepth = pcr.numpy2pcr(Scalar, iceDepth, self.MV)
					iceDepth = pcr.ifthen(self.clone, iceDepth)  #-only use values where clone is True
					pcr.report(iceDepth, self.outpath + 'iceDepth_' + self.curdate.strftime('%Y%m%d') + '.map')
					#-Report pcraster map of glacier volume
					pcr.report(iceDepth * glacArea, self.outpath + 'vIce_' + self.curdate.strftime('%Y%m%d') + '.map')
					#-Delete variables that are not needed
					glacArea = None; del glacArea; iceDepth = None; del iceDepth; GlacTable_MODid = None; del GlacTable_MODid

		#-Water balance
		if self.GlacFLAG and self.GlacRetreat == 1:
			GlacTable_MODid = self.GlacTable.loc[:,['FRAC_GLAC', 'ICE_DEPTH']]
			GlacTable_MODid['ICE_DEPTH'] = GlacTable_MODid['ICE_DEPTH'] * GlacTable_MODid['FRAC_GLAC']
			GlacTable_MODid = GlacTable_MODid.groupby(GlacTable_MODid.index).sum()
			GlacTable_MODid.fillna(0., inplace=True)
			#-Report pcraster map of glacier depth
			iceDepth = pcr.numpy.zeros(self.ModelID_1d.shape)
			iceDepth[self.GlacierKeys] = GlacTable_MODid['ICE_DEPTH']
			iceDepth = iceDepth.reshape(self.ModelID.shape)
			iceDepth = pcr.numpy2pcr(Scalar, iceDepth, self.MV)
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
			waterbalance = Precip - ActETact - self.BaseR - RainR - SnowR - GlacR - dS
		else:
			waterbalance = Precip - ActETact - self.BaseR - RainR - SnowR - dS - self.SeePage
		self.report(waterbalance, self.outpath + 'wbal')
		#-total water balance
		self.waterbalanceTot = self.waterbalanceTot + waterbalance
		self.report(self.waterbalanceTot, self.outpath + 'wbalT')
		#-report water balance and accumulated water balance
		self.wbalTSS.sample(pcr.catchmenttotal(waterbalance, self.FlowDir) / pcr.catchmenttotal(1., self.FlowDir))
		self.wbalTotTSS.sample(pcr.catchmenttotal(self.waterbalanceTot, self.FlowDir) / pcr.catchmenttotal(1., self.FlowDir))
		# set old state variables
		self.oldRootWater = self.RootWater
		self.oldSubWater = self.SubWater
		waterbalance = None; del waterbalance; dS = None; del dS;
		#-End of wate balance calculations
			
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
	
toc = time.clock()
dt = toc - tic
print 'Simulation succesfully completed in '+str(dt)+' seconds!'
