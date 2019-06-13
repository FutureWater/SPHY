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


print('snow module imported')

#-Function to calculate the potential snow melt
def PotSnowMelt(pcr, temp, ddfs):
    melt = pcr.max(0, temp) * ddfs
    return melt
#-Function to calculate the actual snow melt
def ActSnowMelt(pcr, snowstore, potmelt):
    melt = pcr.min(snowstore, potmelt)
    return melt

#-Function that updates the snow storage
def SnowStoreUpdate(pcr, snowstore, snow, actmelt, temp, snowwatstore):
    snowstore = snowstore + snow - actmelt + pcr.ifthenelse(temp < 0, pcr.scalar(snowwatstore), 0)
    return snowstore

#-Function that determines the maximum amount of water that can be stored in the snowpack
def MaxSnowWatStorage(snowsc, snowstore):
    maxsnowwatstore = snowsc * snowstore
    return maxsnowwatstore

#-Function to calculate the actual snow water storage 
def SnowWatStorage(pcr, temp, maxsnowwatstore, snowwatstore, actmelt, rain):
    snowwatstore = pcr.ifthenelse(temp < 0, 0, pcr.min(maxsnowwatstore, snowwatstore + actmelt + rain))
    return snowwatstore

#-Function to calculate the total snow storage (snowstore + snowwatstore)
def TotSnowStorage(snowstore, snowwatstore, snowfrac, rainfrac):
    totalsnowstore = (snowstore + snowwatstore) * (snowfrac + rainfrac)
    return totalsnowstore

#-Function to calculate runoff from snow
def SnowR(pcr, snowwatstore, maxsnowwatstore, actmelt, rain, oldsnowwatstore, snowfrac):
    snowr = pcr.ifthenelse(snowwatstore == maxsnowwatstore, (((actmelt + rain) - (snowwatstore - oldsnowwatstore)) * snowfrac), 0)
    return snowr

#-init snow processes
def init(self, pcr, config):
    pars = ['Tcrit','SnowSC','DDFS']
    for	i in pars:
        try:
            setattr(self, i, pcr.readmap(self.inpath + config.get('SNOW',i)))
        except:
            setattr(self, i, config.getfloat('SNOW',i))

#-initial snow processes
def initial(self, pcr, config):
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

#-dynamic snow processes
def dynamic(self, pcr, Temp, Precip, Snow_GLAC, ActSnowMelt_GLAC, SnowFrac, RainFrac, SnowR_GLAC):
    #-Snow and rain differentiation
    Snow = pcr.ifthenelse(Temp >= self.Tcrit, 0, Precip)
    Rain = pcr.ifthenelse(Temp < self.Tcrit, 0, Precip)
    #-Report Snow for entire cell (snow+glacier fraction)
    self.reporting.reporting(self, pcr, 'TotSnow', Snow)
    self.reporting.reporting(self, pcr, 'TotSnowF', Snow * (1-self.GlacFrac) + Snow_GLAC)
    #-Snow melt
    PotSnowMelt = self.snow.PotSnowMelt(pcr, Temp, self.DDFS)
    ActSnowMelt = self.snow.ActSnowMelt(pcr, self.SnowStore, PotSnowMelt)
    #-Report snow melt for entire cell (snow+glacier fraction)
    self.reporting.reporting(self, pcr, 'TotSnowMelt', ActSnowMelt)
    self.reporting.reporting(self, pcr, 'TotSnowMeltF', ActSnowMelt * (1-self.GlacFrac) + ActSnowMelt_GLAC)
    #-Update snow store
    self.SnowStore = self.snow.SnowStoreUpdate(pcr, self.SnowStore, Snow, ActSnowMelt, Temp, self.SnowWatStore)
    #-Caclulate the maximum amount of water that can be stored in snowwatstore
    MaxSnowWatStore = self.snow.MaxSnowWatStorage(self.SnowSC, self.SnowStore)
    OldSnowWatStore = self.SnowWatStore
    #-Calculate the actual amount of water stored in snowwatstore
    self.SnowWatStore = self.snow.SnowWatStorage(pcr, Temp, MaxSnowWatStore, self.SnowWatStore, ActSnowMelt, Rain)
    #-Changes in total water storage in snow (SnowStore and SnowWatStore)
    OldTotalSnowStore = self.TotalSnowStore
    self.TotalSnowStore = self.snow.TotSnowStorage(self.SnowStore, self.SnowWatStore, SnowFrac, RainFrac) + self.TotalSnowStore_GLAC  # for entire cell
    #-Report Snow storage
    self.reporting.reporting(self, pcr, 'StorSnow', self.TotalSnowStore)
    #-Determine if cell is covered with snow
    SnowCover = pcr.ifthenelse(self.TotalSnowStore > 0, pcr.scalar(1), pcr.scalar(0))
    self.reporting.reporting(self, pcr, 'SCover', SnowCover)
    #-Snow runoff
    SnowR = self.snow.SnowR(pcr, self.SnowWatStore, MaxSnowWatStore, ActSnowMelt, Rain, OldSnowWatStore, SnowFrac) + SnowR_GLAC  # for entire cell
    SnowR = SnowR * (1-self.openWaterFrac)
    #-Report Snow runoff
    self.reporting.reporting(self, pcr, 'TotSnowRF', SnowR)

    return Rain, SnowR, OldTotalSnowStore