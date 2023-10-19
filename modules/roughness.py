# Roughness module
# Copyright (C) 2017-2019 Joris Eekhout / Spanish National Research Council (CEBAS-CSIC)
# Email: jeekhout@cebas.csic.es
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

from modules.travel_time_routing import waterDepth


print('roughness module imported')

#-Manning for tilled conditions (manningTilled; s/m1/3)
def manningTillage(self, pcr):
    manningTilled = pcr.exp(-2.1132 + 0.0349 * self.RFR)
    return manningTilled

#-Manning for vegetated conditions (manningVegetated; s/m1/3)
def manningVegetation(waterDepth, diameter, noElements):
    manningVegetated = (waterDepth**(0.67)) / ((2 * 9.81) / (diameter * noElements))**0.5
    return manningVegetated

#-Determine roughness
def manningField(self, pcr, waterDepth):
    #-Determine Manning's roughness for in-field deposition
    manningHillslopeVegetation = self.roughness.manningVegetation(waterDepth, self.Diameter, self.NoElements)
    manningHillslopeVegetation = pcr.ifthenelse(self.NoVegetation == 1, 0, manningHillslopeVegetation)
    manningHillslopeVegetation = pcr.ifthenelse(self.NoErosion == 1, 0, manningHillslopeVegetation)
    manningHillslopeVegetation = pcr.ifthenelse(self.n_table > 0, self.n_table, manningHillslopeVegetation)
    manningHillslope = (self.n_soil**2 + manningHillslopeVegetation**2)**0.5

    #-Determine Manning's roughness after harvest
    if self.harvest_FLAG:
        manningHillslopeHarvestVegetation = self.roughness.manningVegetation(waterDepth, self.Diameter_harvest, self.NoElements_harvest)
        manningHillslopeHarvestVegetation = pcr.ifthenelse(self.Tillage_harvest == 1, 0, manningHillslopeHarvestVegetation)
        manningHillslopeHarvest = (self.n_soil**2 + manningHillslopeHarvestVegetation**2)**0.5
        manningHillslopeHarvest = pcr.cover(manningHillslopeHarvest, manningHillslope)
    
    return manningHillslope, manningHillslopeHarvest

#-init processes roughness module
def init(self, pcr, config):
    #-read input parameters
    self.n_bare = config.getfloat('EROSION', 'manningBare')
    self.RFR = config.getfloat('EROSION', 'RFR')
    self.n_tilled = pcr.scalar(config.getfloat('EROSION', 'manningTillage'))
    self.n_soil = pcr.ifthenelse(self.Tillage == 1, self.n_tilled, self.n_bare)

#-dynamic processes roughness module
def dynamic(self, pcr):
    #-determine areas that have been harvested
    if self.harvest_FLAG:
        self.Harvested = self.ones * 0
        self.Harvested = pcr.ifthenelse(self.Harvest < self.Sowing, pcr.ifthenelse(pcr.pcrand(self.Harvest < self.curdate.timetuple().tm_yday, self.Sowing > self.curdate.timetuple().tm_yday), 1, self.Harvested), self.Harvested)
        self.Harvested = pcr.ifthenelse(self.Harvest > self.Sowing, pcr.ifthenelse(pcr.pcror(self.curdate.timetuple().tm_yday > self.Harvest, self.curdate.timetuple().tm_yday < self.Sowing), 1, self.Harvested), self.Harvested)
        self.Harvested = pcr.ifthenelse(self.Harvest == 0, 0, self.Harvested)

    #-Determine Manning's roughness for in-field deposition
    self.n_field, self.n_field_harvest = self.roughness.manningField(self, pcr, self.waterDepth)

    #-replace roughness for vegetated conditions for tilled soil conditions in case of harvested areas
    if self.harvest_FLAG:
        self.manningHillslope = pcr.ifthenelse(self.Harvested == 1, self.n_field_harvest, self.n_field)
    else:
        self.manningHillslope = self.n_field
