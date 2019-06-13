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

#-Function to calculate surface runoff
def RootRunoff(self, pcr, rainfrac, rain):
    #-Infiltration excess surface runoff
    if self.InfilFLAG == 1:
        #-Infiltration capacity, scaled based on rootwater content and ksat and corrected for paved surface
        Infil_cap = self.K_eff * self.RootKsat / 24 * (1 + ((self.RootSat - self.RootWater) / self.RootSat))**self.Labda_Infil

        #-Infiltration
        Infil_excess = pcr.ifthenelse((self.Alpha * rain) > Infil_cap, rain - ((self.Alpha * rain - Infil_cap)**2) / (self.Alpha**2 * rain), rain)
        Saturated_excess = self.RootSat - self.RootWater
        Infil = pcr.max(0, pcr.min(Infil_excess, Saturated_excess)) * (1-self.pavedFrac)

        #-Surface runoff
        rootrunoff = rain - Infil

    #-Saturation excess surface runoff
    else:
        #-Assume infiltration capacity to be equal to saturated hydraulic conductivity
        Infil_cap = self.RootKsat

        #-Infiltration
        Infil = pcr.max(0, pcr.min(rain, Infil_cap, self.RootSat - self.RootWater))
        #-Runoff
        rootrunoff = pcr.ifthenelse(rainfrac > 0, rain - Infil, 0)

    return rootrunoff, Infil

#-Function to calculate rootzone drainage
def RootDrainage(pcr, rootwater, rootdrain, rootfield, rootsat, drainvel, rootTT):
    rootexcess = pcr.max(rootwater - rootfield, 0)
    rootexcessfrac = rootexcess / (rootsat - rootfield)
    rootlat = rootexcessfrac * drainvel
    rootdrainage = pcr.max(pcr.min(rootexcess, rootlat * (1-pcr.exp(-1/rootTT)) + rootdrain * pcr.exp(-1/rootTT)), 0)
    return rootdrainage

#-Function to calculate rootzone percolation
def RootPercolation(pcr, rootwater, subwater, rootfield, rootTT, subsat):
    rootexcess = pcr.max(rootwater - rootfield, 0)
    rootperc = rootexcess * (1 - pcr.exp(-1 / rootTT))
    rootperc = pcr.ifthenelse(subwater >= subsat, 0, pcr.min(subsat - subwater, rootperc))
    rootperc = pcr.max(pcr.min(rootperc, rootexcess), 0)
    return rootperc

#-Function to calculate the right fraction between the two fluxes
def CalcFrac(pcr, rootwater, rootfield, rootdrain, rootperc):
    rootexcess = pcr.max(rootwater - rootfield, 0)
    frac = ((rootdrain + rootperc) - rootexcess) / (rootdrain + rootperc)
    rootdrain = rootdrain - (rootdrain * frac)
    rootperc = rootperc - (rootperc * frac)
    return rootdrain, rootperc
