# Re-infiltration module
# Copyright (C) 2021-2023 Joris Eekhout / Spanish National Research Council (CEBAS-CSIC)
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

print('re-infiltration module imported')

#-Determine the flow area (m2) based on the channel storage
def area(self, pcr):
    area = self.channelStorage / pcr.celllength()
    return area

#-Determine bankfull channel area (m2)
def channel_area(self):
    channelArea = self.channelDepth * self.channelWidth
    return channelArea

#-function to determine water height (m) from channel dimensions and channel storage
def water_height(self, pcr, area, channelArea, floodplainFRAC):
    waterHeight = pcr.ifthenelse(floodplainFRAC, (area - channelArea) / self.floodplainWidth + self.channelDepth, area / self.channelWidth)
    return waterHeight

#-function to determine infiltration
def infiltration(pcr, infilPotential, waterHeight):
    Infil = pcr.ifthenelse(infilPotential > waterHeight, waterHeight, infilPotential)
    return Infil 


#-dynamic re-infiltration processes
def dynamic(self, pcr):
    #-Determine the flow area (m2) based on the channel storage
    area = self.reinfiltration.area(self, pcr)

    #-Determine bankfull channel area (m2)
    channelArea = channel_area(self)

    #-Test if there is channel storage
    channelStorageFRAC = self.channelStorage > 0

    #-Test if there is water in the floodplain
    floodplainFRAC = area > channelArea

    #-Determine water height (mm) based on channel dimensions and channel storage
    waterHeight = self.reinfiltration.water_height(self, pcr, area, channelArea, floodplainFRAC) * 1e3

    #-Determine water height above floodplain (mm)
    waterHeightFloodplain = pcr.ifthenelse(floodplainFRAC, waterHeight - self.channelDepth * 1e3, 0)

    #-Assume infiltration capacity to be equal to saturated hydraulic conductivity
    infilCapacity = self.RootKsat

    #-Determine infiltration potential (mm)
    infilPotential = pcr.max(0, pcr.min(infilCapacity, self.RootSat - self.RootWater))

    #-Determine infiltration (mm) in case of floodplain flow
    infilFloodplain = pcr.ifthenelse(floodplainFRAC, self.reinfiltration.infiltration(pcr,infilPotential, waterHeightFloodplain), 0)

    #-Subtract infiltration in floodplain from infiltration potential (mm)
    infilPotential = infilPotential - infilFloodplain

    #-Subtract infiltration in floodplain from water height (mm)
    waterHeight = waterHeight - infilFloodplain

    #-Determine infiltration in the channel section
    infilChannel = self.reinfiltration.infiltration(pcr, infilPotential, waterHeight)

    #-Determine the infiltration volume (m3) based on infiltration in floodplain and channel
    infilVolume = (infilFloodplain * 1e-3 * self.floodplainWidth + infilChannel * 1e-3 * self.channelWidth) * pcr.celllength()

    #-Determine the re-infiltration (mm) per cell based on infiltration volume and cell area
    Infil = infilVolume / pcr.cellarea() * 1e3

    #-Subtract infiltration volume from channel storage (m3)
    self.channelStorage = self.channelStorage - infilVolume

    return Infil