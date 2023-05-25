# The Spatial Processes in HYdrology (SPHY) model:
# A spatially distributed hydrological model 
# Copyright (C) 2013-2023  FutureWater
# Email: sphy@futurewater.nl
#
# Authors (alphabetical order):
# S. Contreras, P. Droogers, J. Eekhout, W. Immerzeel, S. Khanal, A. Lutz, G. Simons, W. Terink
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
    #-read table with irrigation input parameters per landuse class
    pcr.setglobaloption('matrixtable')
    IrrigationTable = self.inpath + config.get('IRRIGATION', 'IrrigationTable')
    self.IrrigationStart = pcr.lookupscalar(IrrigationTable, 1, self.LandUse)
    self.IrrigationEnd = pcr.lookupscalar(IrrigationTable, 2, self.LandUse)
    self.MAD = pcr.lookupscalar(IrrigationTable, 3, self.LandUse)
    pcr.setglobaloption('columntable')

    # Irrigated land use 
    self.IrrigatedLand = pcr.ifthenelse(self.IrrigationStart > 0, pcr.scalar(1), pcr.scalar(0))


#-initial conditions irrigation
def initial(self, pcr):
    # Irrigation input as model runs    
    self.IrrigationWater  = 0


# dynamic processes irrigation
def dynamic(self, pcr, RAW):
    #-determine cells under cultivation
    UnderCultivation = self.ones
    UnderCultivation = pcr.ifthenelse(self.IrrigationEnd < self.IrrigationStart, pcr.ifthenelse(pcr.pcrand(self.IrrigationEnd < self.curdate.timetuple().tm_yday, self.IrrigationStart > self.curdate.timetuple().tm_yday), 0, UnderCultivation), UnderCultivation)
    UnderCultivation = pcr.ifthenelse(self.IrrigationEnd > self.IrrigationStart, pcr.ifthenelse(pcr.pcror(self.curdate.timetuple().tm_yday > self.IrrigationEnd, self.curdate.timetuple().tm_yday < self.IrrigationStart), 0, UnderCultivation), UnderCultivation)
    UnderCultivation = pcr.ifthenelse(self.IrrigationEnd == 0, 0, UnderCultivation)

    #-Computation net irrigation water 
    irrigationThreshold = self.RootField - RAW * self.MAD
    IrrigationWater = pcr.ifthenelse(self.RootWater <= irrigationThreshold, irrigationThreshold - self.RootWater, 0) * UnderCultivation
    self.reporting.reporting(self, pcr, 'TotIrr', IrrigationWater)

    return IrrigationWater

