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

print('irrigation module imported')

#-init processes irrigation
def init(self, pcr, config):
    
#    #Irrigation scenarios
#    self.IrrigationMethod = config.getint('IRRIGATION', 'IrrigationMethod')
    # # Irrigated land use 
    # self.IrrigatedLand = pcr.lookupscalar(self.IrrigationTable,self.LandUse)

    #-read table with irrigation input parameters per landuse class
    pcr.setglobaloption('matrixtable')
    IrrigationTable = self.inpath + config.get('IRRIGATION', 'IrrigationTable')
    self.IrrigationStart = pcr.lookupscalar(IrrigationTable, 1, self.LandUse)
    self.IrrigationEnd = pcr.lookupscalar(IrrigationTable, 2, self.LandUse)
    self.MAD = pcr.lookupscalar(IrrigationTable, 3, self.LandUse)
    pcr.setglobaloption('columntable')

    # Irrigated land use 
    self.IrrigatedLand = pcr.ifthenelse(self.IrrigationStart > 0, pcr.scalar(1), pcr.scalar(0))

#    #Baresoil parameter = 0.15  
#    self.KcBaresoil = config.getfloat('IRRIGATION','KC_baresoil')
#    #Readily available water (RAW) replenishment rate 
#    self.Replenishment  = config.getfloat('IRRIGATION','QUOTA')
#    #Crops water stress coefficient 
#    self.Ks = self.inpath + config.get('IRRIGATION' , 'Ks')
#    self.KsTable = pcr.lookupscalar(self.Ks, self.LandUse)

#-initial conditions irrigation
def initial(self, pcr):
    # Irrigation input as model runs    
    self.IrrigationWater  = 0


# dynamic processes irrigation
def dynamic(self, pcr, RAW): 
    # Irrigation to be applied only when Crops are into the field. 
    # UnderCultivation = pcr.ifthenelse(self.Kc != self.KcBaresoil, self.IrrigatedLand ,0)

    UnderCultivation = self.ones
    UnderCultivation = pcr.ifthenelse(self.IrrigationEnd < self.IrrigationStart, pcr.ifthenelse(pcr.pcrand(self.IrrigationEnd < self.curdate.timetuple().tm_yday, self.IrrigationStart > self.curdate.timetuple().tm_yday), 0, UnderCultivation), UnderCultivation)
    UnderCultivation = pcr.ifthenelse(self.IrrigationEnd > self.IrrigationStart, pcr.ifthenelse(pcr.pcror(self.curdate.timetuple().tm_yday > self.IrrigationEnd, self.curdate.timetuple().tm_yday < self.IrrigationStart), 0, UnderCultivation), UnderCultivation)
    UnderCultivation = pcr.ifthenelse(self.IrrigationEnd == 0, 0, UnderCultivation)

    # #Crop Available water within RootZone 
    # Taw = (self.RootField-self.RootDry)
    # #Depletion Factor 
    # d   = pcr.max(pcr.min(self.PMap + 0.04 * (5 - ETpot), 0.8), 0.1) 
    # # Readily available water 
    # Raw = (Taw * d ) 
    # #Plant and soil specific soil moisture content from which plant water stress starts to occur
    # RootPWS = (self.RootField - Raw)
    # # Crop water Stress coefficient
    # Ks = pcr.max(pcr.min((self.RootWater - self.RootDry) / (RootPWS - self.RootDry),1),0)
    # #Actual Evapotranspiration horticulture area Rambla del Albujon
    # CropETa = ETact * self.IrrigatedLand
    # #Potential Evapotranspiration horticulture area Rambla del Albujon
    # CropETp = ETpot *self.IrrigatedLand

    # #BASELINE SCENARIO 
    # if self.IrrigationMethod == 1:
    # #RAW depletion threshold.  
    # threshold = (self.RootField - Raw*self.MAD ) * self.IrrigatedLand
    # #Delivering Irrigation 
    # IrrigationScheduling = pcr.ifthenelse( self.RootWater <= threshold, Raw*self.Replenishment,0) * self.IrrigatedLand

#    #EVAPOTRANSPIRATION DEFICIT
#    elif self.IrrigationMethod == 2
#       IrrigationScheduling= pcr.ifthenelse(ETpot- ETact >0,ETpot- ETact,0) * self.IrrigatedLand
     
#    #STATIC CROP WATER STRESS THRESHOLD
#    elif self.IrrigationMethod == 3 :

#       IrrigationScheduling = pcr.ifthenelse(Ks <= self.KsTable, Raw*self.Replenishment,0) *self.IrrigatedLand
      
    #-Computation net irrigation water 
    irrigationThreshold = self.RootField - RAW * self.MAD
    IrrigationWater = pcr.ifthenelse(self.RootWater <= irrigationThreshold, irrigationThreshold - self.RootWater, 0) * UnderCultivation
    self.reporting.reporting(self, pcr, 'TotIrr', IrrigationWater)

    # pcr.report(IrrigationWater, self.outpath + "IrrigationWater_" + str(self.counter).zfill(3) + ".map")

#    #Time series irrigation water samount within the Rambla del Albujon 
#    IrrigationTOT= pcr.areatotal((IrrigationWater/1000)*pcr.cellarea(),self.clone)
#    #Time series actual evapotranspiration withimn the Rambla del Albujon 
#    HorticultureETa = pcr.areatotal((CropETa/1000)*pcr.cellarea(),self.clone) 
#    #Time series potential evapotranspiration within the Rambòa del Albujon 
#    HorticultureETp = pcr.areatotal((CropETp/1000)*pcr.cellarea(),self.clone)
    return IrrigationWater

