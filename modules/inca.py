# Soil erosion module using the INCA soil erosion model
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


print('INCA module imported')

#-Detachment of soil particles by raindrop impact (kg/km2)
def DetachmentRaindrop(self, pcr, c_X1, p_Sed, E_SP, V, CG):
    p = p_Sed / (24 * 60 * 60) * 1e-3
    S_SP = pcr.max(0, 1 - CG - self.NoErosion_INCA) * c_X1 * p * E_SP**(10 / (10 - V)) * 8.64 * 1e10
    return S_SP

#-Detachment of soil particles by runoff (kg/km2)
def DetachmentRunoff(self, pcr, q_dr, S_SP, S_TC, E_FL, a1, a2, a3, CG, A, L):
    K = a1 * E_FL * 86400 * ((A * q_dr) / L - a2)**a3
    S_FL = pcr.max(0, pcr.ifthenelse(S_TC + K == 0, 0, K * (S_TC - S_SP) / (S_TC + K)) * pcr.max(0, 1 - self.NoErosion_INCA))

    #-set values in channels to 0 in case channels should be excluded
    if self.exclChannelsFLAG == 1:
        S_FL = S_FL * self.Hillslope
    return S_FL, K

#-Transport capacity of the flow (kg/km2)
def TransportCapacity(self, pcr, q_dr, a4, a5, a6, A, L):
    S_TC = a4 * 86400 * ((A * q_dr) / L - a5)**a6
    return S_TC

#-Update sediment store (kg/km2)
def UpdateSedimentStore(self, pcr, S_SP, S_FL, S_TC, K):
    if self.exclChannelsFLAG == 1:
        dStore_dt = pcr.ifthenelse(self.Hillslope == 1, -K * (S_SP - S_TC) / (S_TC + K), S_SP)
    else:
        dStore_dt = -K * (S_SP - S_TC) / (S_TC + K)
    dStore_dt = pcr.ifthenelse(S_TC + K == 0, 0, dStore_dt)
    
    S_store = pcr.max(0, self.S_store + pcr.ifthenelse(self.S_store + S_SP > S_TC, S_SP - S_TC, dStore_dt))
    return S_store

#-Mass of sediment transported (kg/km2)
def sedimentTransported(self, pcr, S_SP, S_FL, S_TC):
    M_Out = pcr.ifthenelse(self.S_store + S_SP > S_TC, S_TC, S_SP + S_FL) * pcr.max(0, 1 - self.NoErosion_INCA)
    return M_Out

#-init processes inca
def init(self, pcr, config):
    #-read table with INCA landuse specific model parameters
    pcr.setglobaloption('matrixtable')
    inca_table = self.inpath + config.get('INCA', 'inca_table')
    self.V_INCA = pcr.lookupscalar(inca_table, 1, self.LandUse)
    self.GC_INCA = pcr.lookupscalar(inca_table, 2, self.LandUse)
    try:
        self.a4_INCA = config.getfloat('INCA', 'a4')
    except:
        self.a4_INCA = pcr.lookupscalar(inca_table, 3, self.LandUse)
    self.NoErosion_INCA = pcr.lookupscalar(inca_table, 4, self.LandUse)
    pcr.setglobaloption('columntable')

    #-read other model parameters
    self.c_x1_INCA = config.getfloat('INCA', 'c_x1')
    self.a1_INCA = config.getfloat('INCA', 'a1')
    self.a2_INCA = config.getfloat('INCA', 'a2')
    self.a3_INCA = config.getfloat('INCA', 'a3')
    self.a5_INCA = config.getfloat('INCA', 'a5')
    self.a6_INCA = config.getfloat('INCA', 'a6')
    self.E_SP_INCA = config.getfloat('INCA', 'E_SP')
    self.E_FL_INCA = config.getfloat('INCA', 'E_FL')

    #-Determine cell size and slope length
    self.A_INCA = pcr.cellarea() * 1e-6
    alpha = pcr.atan(self.Slope)
    self.L_INCA = pcr.celllength() * 1e-3 / pcr.cos(alpha)

    #-Initiate sediment store variable
    self.S_store = self.DEM * 0

#-dynamic processes inca
def dynamic(self, pcr, Precip, Q):
    #-determine canopy cover from LAI
    if self.DynVegFLAG == 1:
        V_INCA = pcr.min(9.999, pcr.min(1, self.LAI) * 10)
    else:
        V_INCA = self.V_INCA

    #-determine detachment of soil particles by raindrop impact (kg/km2)
    S_SP = self.inca.DetachmentRaindrop(self, pcr, self.c_x1_INCA, Precip, self.E_SP_INCA, V_INCA, self.GC_INCA)

    #-report detachment of soil particles by raindrop impact (ton / cell)
    self.reporting.reporting(self, pcr, 'DetRn', S_SP * pcr.cellarea() * 1e-9)

    #-Determine transport capacity of the flow (kg/km2)
    S_TC = self.inca.TransportCapacity(self, pcr, Q, self.a4_INCA, self.a5_INCA, self.a6_INCA, self.A_INCA, self.L_INCA)

    #-determine detachment of soil particles by runoff (kg/km2)
    S_FL, K = self.inca.DetachmentRunoff(self, pcr, Q, S_SP, S_TC, self.E_FL_INCA, self.a1_INCA, self.a2_INCA, self.a3_INCA, self.GC_INCA, self.A_INCA, self.L_INCA)

    #-report detachment of soil particles by runoff (ton / cell)
    self.reporting.reporting(self, pcr, 'DetRun', S_FL * pcr.cellarea() * 1e-9)

    #-determine mass of sediment in transport (kg/km2)
    sed = self.inca.sedimentTransported(self, pcr, S_SP, S_FL, S_TC)

    #-update sediment store (kg/km2)
    self.S_store = self.inca.UpdateSedimentStore(self, pcr, S_SP, S_FL, S_TC, K)

    #-report sediment in transport (ton / cell)
    self.reporting.reporting(self, pcr, 'SedTrans', sed * pcr.cellarea() * 1e-9)

    return sed
