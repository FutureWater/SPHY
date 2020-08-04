# Soil erosion module using the DHSVM soil erosion model
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


print('DHSVM module imported')

#-Detachment of soil particles by raindrop impact (kg/m2/s)
def DetachmentRaindrop(self, pcr, k_r, F_w, C_G, C_C, M_R, M_D):
    D_R = k_r * F_w * pcr.max(0, 1 - C_G - self.NoErosion_DHSVM) * (pcr.max(0, 1 - C_C) * M_R + M_D)
    return D_R

#-Rainfall drop momentum (kg2/s3)
def MomentumRainDrop(self, pcr, np, precSum):
    #-store daily precipitation
    precipDaily = precSum

    #-determine sum over all raster cells 
    precIntTotal = np.sum(pcr.pcr2numpy(precSum, 0))

    #-initiate hour counter and set M_r to 0
    hour = 0
    M_R = self.DEM * 0

    #-while precipitation sum is larger than 0
    while (precIntTotal > 0):
        #-determine hourly rainfall intensity
        I = pcr.max(0, self.Alpha * precipDaily - ((self.Alpha**2 * precipDaily) / 2) * hour)

        #-determine coefficient alpha
        alpha = pcr.scalar(I >= 100) * 11.75 * 1e-8
        alpha = alpha + pcr.scalar(pcr.pcrand(I >= 50, I < 100)) * 6.12 * 1e-8
        alpha = alpha + pcr.scalar(pcr.pcrand(I >= 10, I < 50)) * 3.75 * 1e-8
        alpha = alpha + pcr.scalar(I < 10) * 2.69 * 1e-8

        #-determine coefficient beta
        beta = pcr.scalar(I >= 100) * 1.2821
        beta = beta + pcr.scalar(pcr.pcrand(I >= 50, I < 100)) * 1.4242
        beta = beta + pcr.scalar(pcr.pcrand(I >= 10, I < 50)) * 1.5545
        beta = beta + pcr.scalar(I < 10) * 1.6896

        #-calculate rainfall drop momentum based on hourly precipitation intensity and multiply with 1/24
        M_R_hour = alpha * I**beta * 1/24

        #-add to daily rainfall drop momentum
        M_R = M_R + M_R_hour

        #-determine sum over all raster cells
        precIntTotal = np.sum(pcr.pcr2numpy(I, 0))

        #-increase hour with 1
        hour += 1
    return M_R

#-Leaf drip momentum (kg2 / s3)
def MomentumLeafDrip(self, pcr, D, X, rho, DRAIN, g, C_C):
    #-determine coefficient a
    a = pcr.scalar(pcr.pcrand(D >= 0.033, X >= 7.5)) * 5.14
    a = a + pcr.scalar(pcr.pcrand(D >= 0.033, X < 7.5)) * 1.93
    a = a + pcr.scalar(D < 0.033) * 0

    #-determine coefficient b
    b = pcr.scalar(pcr.pcrand(D >= 0.033, X >= 7.5)) * 660
    b = b + pcr.scalar(pcr.pcrand(D >= 0.033, X < 7.5)) * 1640
    b = b + pcr.scalar(D < 0.033) * 2200

    #-Ratio M/beta (m)
    M_beta = a + b * D
    #-Leaf drip fall speed (m/s)
    V = (M_beta * g * (1 - pcr.exp(-2*X/M_beta)))**0.5
    #-Set the proportion of the drainage that falls as leaf drip equal to the canopy cover
    DRIP = C_C
    #-Determine the momentum squared for leaf drip (kg2 / s3)
    M_D = (((V * rho * self.pi * D**3) / 6)**2 * DRIP * DRAIN) / ((self.pi * D**3) / 6)
    return M_D

#-Detachment of soil particles by runoff (kg/m2/s)
def DetachmentRunoff(self, pcr, beta_de, dy, v_s, TC, C_G):
    D_of = beta_de * dy * v_s * TC * pcr.max(0, 1 - self.NoErosion_DHSVM)
    
    #-set values in channels to 0 in case channels should be excluded
    if self.exclChannelsFLAG == 1:
        D_of = D_of * self.Hillslope
    return D_of

#-Detachment efficiency (-)
def DetachmentEfficiency(self, pcr, C_s):
    beta_de = 0.79 * pcr.exp(-0.6 * C_s)
    return beta_de

#-Settling velocity (m/s)
def SettlingVelocity(self, pcr, g, rho_s, rho, nu, d50):
    #-Initial guess of the settling velocity (m/s)
    v_s = (4/3 * g * (rho_s / rho) - 1)**0.5 * d50
    #-Reynolds number (-)
    Rn = (v_s * d50) / nu
    #-Drag coefficient (-)
    Cd = (24 / Rn) + (3 / (Rn**0.5)) + 0.34
    #-Settling velocity (m/s)
    v_s = ((4/3 * g * ((rho_s/rho) - 1) * d50) / Cd)**0.5
    return v_s

#-Stream power (kg/m/s3)
def StreamPower(self, pcr, rho, g, Q, S):
    SP = rho * g * Q * S
    return SP

#-Transport capacity (m3/m3)
def TransportCapacity(self, pcr, d50, rho, rho_s, S, h, g, SP, SP_cr):
    TC = (0.05 / (d50 * (rho_s / rho - 1)**2)) * ((S * h) / g)**0.5 * pcr.max(0, (SP - SP_cr))
    return TC

#-Water depth (m), flow velocity (m/s)
def Manning(self, pcr, Q, n, WD_ratio, S):
    #-Determine flow area and water depth
    A = ((Q * n * (2 * ((WD_ratio**2 + 1) / WD_ratio)**0.5)**(2/3)) / S**0.5)**(3/4)
    h = (A / WD_ratio)**0.5

    #-Set minimum value for flow area and water depth
    h = pcr.max(self.h_min_DHSVM, h)
    return h

#-init processes dhsvm
def init(self, pcr, config):
    #-read table with DHSVM landuse specific model parameters
    pcr.setglobaloption('matrixtable')
    dhsvm_table = self.inpath + config.get('DHSVM', 'dhsvm_table')
    self.D_DHSVM = pcr.lookupscalar(dhsvm_table, 1, self.LandUse)
    self.X_DHSVM = pcr.lookupscalar(dhsvm_table, 2, self.LandUse)
    self.C_G_DHSVM = pcr.lookupscalar(dhsvm_table, 3, self.LandUse)
    self.C_C_table_DHSVM = pcr.lookupscalar(dhsvm_table, 4, self.LandUse)
    self.n_table_DHSVM = pcr.lookupscalar(dhsvm_table, 5, self.LandUse)
    self.root_cohesion_DHSVM = pcr.lookupscalar(dhsvm_table, 6, self.LandUse)
    self.NoErosion_DHSVM = pcr.lookupscalar(dhsvm_table, 7, self.LandUse)

    #-read table with soil cohesion per soil class
    self.SoilClass_DHSVM = pcr.readmap(self.inpath + config.get('DHSVM', 'SoilClass'))
    dhsvm_cohesion_table = self.inpath + config.get('DHSVM', 'dhsvm_cohesion_table')
    self.soil_cohesion_DHSVM = pcr.lookupscalar(dhsvm_cohesion_table, 1, self.SoilClass_DHSVM)
    pcr.setglobaloption('columntable')

    #-determine overall cohesion from soil and root cohesion
    self.C_s_DHSVM = self.soil_cohesion_DHSVM + self.root_cohesion_DHSVM 

    #-read other model parameters
    self.WD_ratio_DHSVM = config.getfloat('DHSVM', 'WD_ratio')
    self.rho_DHSVM = config.getfloat('DHSVM', 'rho')
    self.rho_s_DHSVM = config.getfloat('DHSVM', 'rho_s')
    self.deltaClay_DHSVM = config.getfloat('DHSVM', 'deltaClay') * 1e-6
    self.deltaSilt_DHSVM = config.getfloat('DHSVM', 'deltaSilt') * 1e-6
    self.deltaSand_DHSVM = config.getfloat('DHSVM', 'deltaSand') * 1e-6
    self.k_r_DHSVM = config.getfloat('DHSVM', 'k_r')
    self.h_min_DHSVM = config.getfloat('DHSVM', 'h_min')
    self.SP_crit_DHSVM = config.getfloat('DHSVM', 'SP_crit')

    #-define constants
    self.F_w_DHSVM = 1
    self.g_DHSVM = 9.81
    self.visc_DHSVM = 1e-06

    #-determine median grain size
    if self.PedotransferFLAG == 1:
        self.D50_DHSVM = pcr.ifthenelse(self.RootClayMap > 0.5, pcr.scalar(self.deltaClay_DHSVM), 0)
        self.D50_DHSVM = pcr.ifthenelse(self.RootClayMap + self.RootSiltMap > 0.5, self.deltaClay_DHSVM + (self.deltaSilt_DHSVM - self.deltaClay_DHSVM) * (0.5 - self.RootClayMap) / self.RootSiltMap, 0) + self.D50_DHSVM
        self.D50_DHSVM = pcr.ifthenelse(self.RootClayMap + self.RootSiltMap < 0.5, self.deltaSilt_DHSVM + (self.deltaSand_DHSVM - self.deltaSilt_DHSVM) * (self.RootSandMap - 0.5) / self.RootSandMap, 0) + self.D50_DHSVM
    else:
        self.D50_DHSVM = config.getfloat('DHSVM', 'D50') * 1e-6

#-dynamic processes dhsvm
def dynamic(self, pcr, np, Precip, Q):
    #-determine canopy cover from LAI
    if self.DynVegFLAG == 1:
        C_C_DHSVM = pcr.min(1, self.LAI)
    else:
        C_C_DHSVM = self.C_C_table_DHSVM

    #-determine rainfall drop momentum (kg2/s3)
    M_R = self.dhsvm.MomentumRainDrop(self, pcr, np, Precip)

    #-determine water drainage rate from canopy (m/s)
    DRAIN = Precip * 1e-3 / (24 * 60 * 60)

    #-determine leaf drop momentum (kg2/s3)
    M_D = self.dhsvm.MomentumLeafDrip(self, pcr, self.D_DHSVM, self.X_DHSVM, self.rho_DHSVM, DRAIN, self.g_DHSVM, C_C_DHSVM)

    #-determine detachment of soil particles by raindrop impact (kg/m2/s)
    D_R = self.dhsvm.DetachmentRaindrop(self, pcr, self.k_r_DHSVM, self.F_w_DHSVM, self.C_G_DHSVM, C_C_DHSVM, M_R, M_D)

    #-report detachment of soil particles by raindrop impact (ton / cell)
    self.reporting.reporting(self, pcr, 'DetRn', D_R * pcr.cellarea() * 1e-3 * (24 * 60 * 60))

    #-determine water depth (m) and width of the flow (m)
    h = self.dhsvm.Manning(self, pcr, Q, self.n_table_DHSVM, self.WD_ratio_DHSVM, self.Slope)

    #-determine detachment efficiency (-)
    beta_de = self.dhsvm.DetachmentEfficiency(self, pcr, self.C_s_DHSVM)

    #-determine settling velocity (m/s)
    v_s = self.dhsvm.SettlingVelocity(self, pcr, self.g_DHSVM, self.rho_s_DHSVM, self.rho_DHSVM, self.visc_DHSVM, self.D50_DHSVM)

    #-determine stream power (kg m/s3)
    SP = self.dhsvm.StreamPower(self, pcr, self.rho_DHSVM, self.g_DHSVM, Q, self.Slope)

    #-determine transport capacity (m3/m3)
    TC = self.dhsvm.TransportCapacity(self, pcr, self.D50_DHSVM, self.rho_DHSVM, self.rho_s_DHSVM, self.Slope, h, self.g_DHSVM, SP, self.SP_crit_DHSVM)

    #-determine detachment of soil particles by runoff (kg/m2/s)
    D_of = self.dhsvm.DetachmentRunoff(self, pcr, beta_de, pcr.celllength(), v_s, TC, self.C_G_DHSVM)

    #-report detachment of soil particles by runoff (ton / cell)
    self.reporting.reporting(self, pcr, 'DetRun', D_of * pcr.cellarea() * 1e-3 * (24 * 60 * 60))

    #-determine mass of sediment in transport (kg/m2/s)
    sed = D_R + D_of

    #-report sediment in transport (ton / cell)
    self.reporting.reporting(self, pcr, 'SedTrans', sed * pcr.cellarea() * 1e-3 * (24 * 60 * 60))

    return sed
