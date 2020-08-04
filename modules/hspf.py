# Soil erosion module using the HSPF soil erosion model
# Copyright (C) 2020 Joris Eekhout / Spanish National Research Council (CEBAS-CSIC)
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


print('HSPF module imported')

#-Detachment of soil particles by raindrop impact (ton/acre)
def DetachmentRaindrop(self, pcr, DELT60, CR, SMPF, KRER, RAIN, JRER):
    DET = DELT60 * pcr.max(0, (1 - CR - self.NoErosion_HSPF)) * SMPF * KRER * (RAIN / DELT60)**JRER
    return DET

#-Detached sediment in storage (ton/acre)
def SedimentStorage(self, pcr, DETS, AFFIX, DET):
    DETS = DETS * (1 - AFFIX) + DET

    return DETS

#-Detachment of soil particles by washoff (ton/acre)
def DetachmentWashoff(self, pcr, STCAP, DETS, SURO, SURS, CR):
    WSSD = pcr.ifthenelse(STCAP > DETS, DETS * SURO / (SURS + SURO), STCAP * SURO / (SURS + SURO)) * (1 - self.NoErosion_HSPF)

    return WSSD

#-Detachment of soil particles from the soil matrix (ton/acre)
def DetachmentSoilScour(self, pcr, SURO, SURS, DELT60, KGER, JGER):
    SCRSD = SURO / (SURS + SURO) * DELT60 * KGER * ((SURS + SURO)/DELT60)**JGER * (1 - self.NoErosion_HSPF)

    #-set values in channels to 0 in case channels should be excluded
    if self.exclChannelsFLAG == 1:
        SCRSD = SCRSD * self.Hillslope

    return SCRSD

#-Transport capacity (ton/acre)
def TransportCapacity(self, pcr, DELT60, KSER, SURO, SURS, JSER):
    STCAP = DELT60 * KSER * ((SURS + SURO)/DELT60)**JSER * (1 - self.NoErosion_HSPF)

    return STCAP

#-K factor (-)
def K_HSPF(self, pcr):
    ksat_hourly = self.RootKsat / 24
    M_textural = (self.RootSiltMap * 100 + 0) * (100 - self.RootClayMap * 100)
    permeability = pcr.scalar(ksat_hourly > 150) * 1
    permeability = permeability + pcr.scalar(pcr.pcrand(ksat_hourly > 50, ksat_hourly < 150)) * 2
    permeability = permeability + pcr.scalar(pcr.pcrand(ksat_hourly > 15, ksat_hourly < 50)) * 3
    permeability = permeability + pcr.scalar(pcr.pcrand(ksat_hourly > 5, ksat_hourly < 15)) * 4
    permeability = permeability + pcr.scalar(pcr.pcrand(ksat_hourly > 1, ksat_hourly < 5)) * 5
    permeability = permeability + pcr.scalar(ksat_hourly < 1) * 6
    s = 2
    K_HSPF = ((2.1 * 10**-4 * M_textural**1.14 * (12 - self.RootOMMap) + 3.25 * (s - 2) + 2.5 * (permeability - 3))/100)
    return K_HSPF

#-init processes hspf
def init(self, pcr, config):
    #-read table with HSPF landuse specific model parameters
    pcr.setglobaloption('matrixtable')
    hspf_table = self.inpath + config.get('HSPF', 'hspf_table')
    self.CR_HSPF = pcr.lookupscalar(hspf_table, 1, self.LandUse)
    try:
        self.KGER_HSPF = config.getfloat('HSPF', 'KGER')
    except:
        self.KGER_HSPF = pcr.lookupscalar(hspf_table, 2, self.LandUse)
    self.NoErosion_HSPF = pcr.lookupscalar(hspf_table, 3, self.LandUse)
    pcr.setglobaloption('columntable')

    #-read other model parameters
    self.JRER_HSPF = config.getfloat('HSPF', 'JRER')
    self.KSER_HSPF = config.getfloat('HSPF', 'KSER')
    self.JSER_HSPF = config.getfloat('HSPF', 'JSER')
    self.JGER_HSPF = config.getfloat('HSPF', 'JGER')
    self.AFFIX_HSPF = config.getfloat('HSPF', 'AFFIX')

    #-initial sediment storage
    self.DETS_HSPF = 0
    self.SURSold = 0
    
    #-define some constants
    self.acre_m2_HSPF = 4046.9
    self.inch_mm_HSPF = 25.4

    #-read P-factor values map or float
    try:
        self.P_HSPF = pcr.readmap(self.inpath + config.get('HSPF', 'P_USLE'))
    except:
        self.P_HSPF = config.getfloat('HSPF', 'P_USLE')

    #-when pedotransfer module is used, calculate the K-factor based on texture maps, else read K-factor values from the config file
    if self.PedotransferFLAG == 1:
        self.K_HSPF = self.hspf.K_HSPF(self, pcr)
    else:
        try:
            self.K_HSPF = pcr.readmap(self.inpath + config.get('HSPF', 'KRER'))
        except:
            self.K_HSPF = config.getfloat('HSPF', 'KRER')


#-dynamic processes hspf
def dynamic(self, pcr, np, Precip, Runoff):
    #-determine daily precipitation in inch
    Precip_inch = Precip / self.inch_mm_HSPF

    #-get acre to m2 transformation value
    acre_m2_HSPF = self.acre_m2_HSPF

    #-determine detachment of soil particles by raindrop impact (ton/acre)
    DET = self.hspf.DetachmentRaindrop(self, pcr, 24, self.CR_HSPF, self.P_HSPF, self.K_HSPF, Precip_inch, self.JRER_HSPF)

    #-determine the sediment storage (ton/acre)
    self.DETS_HSPF = self.hspf.SedimentStorage(self, pcr, self.DETS_HSPF, self.AFFIX_HSPF, DET)

    #-determine the surface water storage
    SURS = 0

    #-determine the surface outflow (inch) (= routed runoff)
    SURO = pcr.max(0.0001, Runoff / self.inch_mm_HSPF)
 
    #-determine transport capacity (ton/acre)
    STCAP = self.hspf.TransportCapacity(self, pcr, 24, self.KSER_HSPF, SURO, SURS, self.JSER_HSPF)

    #-determine detachment of soil particles by washoff (ton/acre)
    WSSD = self.hspf.DetachmentWashoff(self, pcr, STCAP, self.DETS_HSPF, SURO, SURS, self.CR_HSPF)

    #-report detachment of soil particles by raindrop impact (ton / cell)
    self.reporting.reporting(self, pcr, 'DetRn', WSSD * (pcr.cellarea() / acre_m2_HSPF))

    #-update sediment storage
    self.DETS_HSPF = self.DETS_HSPF - WSSD

    #-determine detachment of soil particles from the soil matrix (ton/acre)
    SCRSD = self.hspf.DetachmentSoilScour(self, pcr, SURO, SURS, 24, self.KGER_HSPF, self.JGER_HSPF)

    #-report detachment of soil particles by runoff (ton / cell)
    self.reporting.reporting(self, pcr, 'DetRun', SCRSD * (pcr.cellarea() / acre_m2_HSPF))

    #-determine mass of sediment in transport (ton/acre)
    sed = WSSD + SCRSD

    #-report sediment in transport (ton / cell)
    self.reporting.reporting(self, pcr, 'SedTrans', sed * (pcr.cellarea() / acre_m2_HSPF))

    return sed
