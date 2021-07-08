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

def infiltration(pcr, infilPotential, waterHeight):
    Infil = pcr.ifthenelse(infilPotential > waterHeight, waterHeight, infilPotential)
    return Infil 


# #-init re-infiltration processes
# def init(self, pcr, pcrm, config, np):

# #-initial conditions re-infiltration
# def initial(self, pcr, config):
    

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
    # pcr.report(waterHeight, self.outpath + "h_reinfiltration_" + str(self.counter).zfill(3) + ".map")

    #-Determine water height above floodplain (mm)
    waterHeightFloodplain = pcr.ifthenelse(floodplainFRAC, waterHeight - self.channelDepth * 1e3, 0)
    # pcr.report(waterHeightFloodplain, self.outpath + "waterHeightFloodplain_" + str(self.counter).zfill(3) + ".map")

    #-Assume infiltration capacity to be equal to saturated hydraulic conductivity
    infilCapacity = self.RootKsat

    #-Determine infiltration potential (mm)
    infilPotential = pcr.max(0, pcr.min(infilCapacity, self.RootSat - self.RootWater))
    # pcr.report(infilPotential, self.outpath + "infilPotential_1_" + str(self.counter).zfill(3) + ".map")

    #-Determine infiltration (mm) in case of floodplain flow
    infilFloodplain = pcr.ifthenelse(floodplainFRAC, self.reinfiltration.infiltration(pcr,infilPotential, waterHeightFloodplain), 0)
    # pcr.report(infilFloodplain, self.outpath + "infilFloodplain_" + str(self.counter).zfill(3) + ".map")

    #-Subtract infiltration in floodplain from infiltration potential (mm)
    infilPotential = infilPotential - infilFloodplain
    # pcr.report(infilPotential, self.outpath + "infilPotential_2_" + str(self.counter).zfill(3) + ".map")

    #-Subtract infiltration in floodplain from water height (mm)
    waterHeight = waterHeight - infilFloodplain
    # pcr.report(waterHeight, self.outpath + "waterHeightInfil_" + str(self.counter).zfill(3) + ".map")

    #-Determine infiltration in the channel section
    infilChannel = self.reinfiltration.infiltration(pcr, infilPotential, waterHeight)
    # pcr.report(infilChannel, self.outpath + "infilChannel_" + str(self.counter).zfill(3) + ".map")

    #-Determine the infiltration volume (m3) based on infiltration in floodplain and channel
    infilVolume = (infilFloodplain * 1e-3 * self.floodplainWidth + infilChannel * 1e-3 * self.channelWidth) * pcr.celllength()
    # pcr.report(infilVolume, self.outpath + "infilVolume_" + str(self.counter).zfill(3) + ".map")

    # #-Determine fraction of cell occupied by water
    # waterFRAC = pcr.ifthenelse(floodplainFRAC, self.floodplainWidth, self.channelWidth) / pcr.celllength() * pcr.scalar(channelStorageFRAC)
    # pcr.report(waterFRAC, self.outpath + "waterFRAC_" + str(self.counter).zfill(3) + ".map")

    #-Determine the re-infiltration (mm) per cell based on infiltration volume and cell area
    Infil = infilVolume / pcr.cellarea() * 1e3
    # pcr.report(Infil, self.outpath + "Infil_" + str(self.counter).zfill(3) + ".map")

    #-Subtract infiltration volume from channel storage (m3)
    # pcr.report(self.channelStorage, self.outpath + "channelStorage_1_" + str(self.counter).zfill(3) + ".map")
    self.channelStorage = self.channelStorage - infilVolume
    # pcr.report(self.channelStorage, self.outpath + "channelStorage_2_" + str(self.counter).zfill(3) + ".map")

    return Infil