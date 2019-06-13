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


# Equations to calculate sediment yield accoding to the soil loss equation (Williams, 1995)
# from rootzone import RootRunoff

print('MUSLE module imported')

#-Modified unviversal soil loss equation to calculate sediment yield (metric tons)
def MUSLE(self, pcr, Q_surf, q_peak):
    sed = 11.8 * (Q_surf * q_peak * self.ha_area)**0.56 * self.K_USLE * self.C_USLE * \
        self.P_USLE * self.LS_USLE * self.CFRG
    sed = (sed / 10000) * pcr.cellarea()  # conversion to ton / cell
    return sed

#-Peak runoff (m3/s)
def q_peak(self, pcr, Runoff):
    q_peak = (self.Alpha_tc * Runoff * pcr.cellarea()) / (3600000 * self.Tc)  # peak runoff in m3/s
    return q_peak

#-LS topographic factor
def LS_ULSE(self, pcr):
    m = 0.6 * (1 - pcr.exp(-35.835 * self.Slope))
    alpha_hill = pcr.atan(self.Slope)
    L_hill = pcr.celllength() / pcr.cos(alpha_hill)
    LS = ((L_hill / 22.1)**m) * (65.41 * pcr.sin(alpha_hill)**2 + 4.56 * pcr.sin(alpha_hill) + 0.065)
    return LS

#-K factor
def K_USLE(self, pcr):
    # ksat_hourly = self.RootKsat / 24 / self.RootKsatFrac
    ksat_hourly = self.RootKsat / 24
    M_textural = (self.RootSiltMap * 100 + 0) * (100 - self.RootClayMap * 100)
    permeability = pcr.scalar(ksat_hourly > 150) * 1
    permeability = permeability + pcr.scalar(pcr.pcrand(ksat_hourly > 50, ksat_hourly < 150)) * 2
    permeability = permeability + pcr.scalar(pcr.pcrand(ksat_hourly > 15, ksat_hourly < 50)) * 3
    permeability = permeability + pcr.scalar(pcr.pcrand(ksat_hourly > 5, ksat_hourly < 15)) * 4
    permeability = permeability + pcr.scalar(pcr.pcrand(ksat_hourly > 1, ksat_hourly < 5)) * 5
    permeability = permeability + pcr.scalar(ksat_hourly < 1) * 6
    s = 2
    K_USLE = ((2.1 * 10**-4 * M_textural**1.14 * (12 - self.RootOMMap) + 3.25 * (s - 2) + 2.5 * (permeability - 3))/100) * 0.1317
    return K_USLE

#-Time of concentration
def Tc(self, pcr):
    slope_adjusted = pcr.ifthenelse(self.Slope < 0.0025, self.Slope + 0.0005, self.Slope)
    alpha_hill = pcr.atan(slope_adjusted)
    L = pcr.celllength() / pcr.cos(alpha_hill)
    #-Kirpich channel flow time of concentration
    T_ch = 0.0195 * L**0.77 * slope_adjusted**-0.385
    #-Kerby overland flow time of concentration
    T_ov = 1.44 * (L * self.N)**0.467 * slope_adjusted**-0.235
    Tc = (T_ch + T_ov) / 60  # conversion to hours
    return Tc

#-init processes musle
def init(self, pcr, config):
    #-read table with MUSLE C-factor and retardance coefficient values per landuse class
    pcr.setglobaloption('matrixtable')
    musle_table = self.inpath + config.get('MUSLE', 'musle_table')
    self.C_USLE = pcr.lookupscalar(musle_table, 1, self.LandUse)
    self.N = pcr.lookupscalar(musle_table, 2, self.LandUse)
    pcr.setglobaloption('columntable')

    #-read P-factor values map or float
    try:
        self.P_USLE = pcr.readmap(self.inpath + config.get('MUSLE', 'P_USLE'))
    except:
        self.P_USLE = config.getfloat('MUSLE', 'P_USLE')

    #-when pedotransfer module is used, calculate the K-factor based on texture maps, else read K-factor values from table
    if self.PedotransferFLAG == 1:
        self.K_USLE = self.musle.K_USLE(self, pcr)
    else:
        self.K_USLE = pcr.readmap(self.inpath + config.get('MUSLE', 'K_USLE'))

    #-calculate other input parameters
    self.LS_USLE = self.musle.LS_ULSE(self, pcr)
    self.CFRG = pcr.exp(-0.053 * (self.RockFrac * 100))
    self.Tc = self.musle.Tc(self, pcr)
    self.Alpha_tc = 1 - pcr.exp(2 * self.Tc * pcr.ln(1 - (self.Alpha/2)))
    self.ha_area = pcr.cellarea() / 10000

#-dynamic processes musle
def dynamic(self, pcr, Runoff):
    #-determine peak runoff
    q_peak = self.musle.q_peak(self, pcr, Runoff) #-peak runoff in m3/s

    #-determine soil erosion
    sed = self.musle.MUSLE(self, pcr, Runoff, q_peak) #-sediment yield in ton

    #-report sediment in transport (ton / cell)
    self.reporting.reporting(self, pcr, 'SedTrans', sed)
