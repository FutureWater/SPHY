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

#-Function to calculate the potential evapotranspiration
def ETpot(etr, kc):
    etpot = etr * kc
    return etpot

#-Function to calculate the actual evapotranspiration
def ETact(pcr, etpot, rootwater, rootsat, etreddry, rainfrac):
    etredwet = pcr.ifthenelse(rootwater >= rootsat, pcr.scalar(0), 1)
    etact = pcr.ifthenelse(rainfrac > 0, pcr.min(etpot * etreddry * etredwet, rootwater), 0)
    return etact

#-Determine plant water stress for calculation of actual evapotranspiration
def ks(self, pcr, etpot):
    TAW = (self.RootField - self.RootDry)
    p = pcr.max(pcr.min(self.PMap + 0.04 * (5 - etpot), 0.8), 0.1)
    RAW = TAW * p
    RootPWS = self.RootField - RAW
    Ks = pcr.max(pcr.min((self.RootWater - self.RootDry) / (RootPWS - self.RootDry), 1), 0)
    return Ks
