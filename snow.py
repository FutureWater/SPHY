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
__version__ = "2.0"
__email__ = "info@futurewater.nl"
__date__ ='1 January 2017'
############################################################################################

print 'snow module imported'

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
