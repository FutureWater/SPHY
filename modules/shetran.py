# Soil erosion module using the SHETRAN soil erosion model
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


print('SHETRAN module imported')

#-Detachment of soil particles by raindrop impact (kg/m2/s)
def DetachmentRaindrop(self, pcr, k_r, F_w, C_g, C_r, M_r, M_d):
    D_r = k_r * F_w * pcr.max(0, 1 - C_g - C_r - self.NoErosion_SHETRAN) * (M_r + M_d)
    return D_r

#-Rainfall drop momentum (kg2/s3)
def MomentumRainDrop(self, pcr, np, C_c, precSum):
    #-store daily precipitation
    precipDaily = precSum

    #-determine sum over all raster cells 
    precIntTotal = np.sum(pcr.pcr2numpy(precSum, 0))

    #-initiate hour counter and set M_r to 0
    hour = 0
    M_r = self.DEM * 0

    #-while precipitation sum is larger than 0
    while (precIntTotal > 0):
        #-determine hourly rainfall intensity
        I = pcr.max(0, self.Alpha * precipDaily - ((self.Alpha**2 * precipDaily) / 2) * hour)

        a_1 = pcr.scalar(I >= 100) * 11.737 * 1e-8
        a_1 = a_1 + pcr.scalar(pcr.pcrand(I >= 50, I < 100)) * 6.1192 * 1e-8
        a_1 = a_1 + pcr.scalar(pcr.pcrand(I >= 10, I < 50)) * 3.7514 * 1e-8
        a_1 = a_1 + pcr.scalar(I < 10) * 2.6893 * 1e-8

        b_1 = pcr.scalar(I >= 100) * 1.2821
        b_1 = b_1 + pcr.scalar(pcr.pcrand(I >= 50, I < 100)) * 1.4242
        b_1 = b_1 + pcr.scalar(pcr.pcrand(I >= 10, I < 50)) * 1.5545
        b_1 = b_1 + pcr.scalar(I < 10) * 1.6896

        #-calculate rainfall drop momentum based on hourly precipitation intensity and multiply with 1/24
        M_r_hour = (1 - C_c) * a_1 * I**b_1 * 1/24

        #-add to daily rainfall drop momentum
        M_r = M_r + M_r_hour

        #-determine sum over all raster cells
        precIntTotal = np.sum(pcr.pcr2numpy(I, 0))

        #-increase hour with 1
        hour += 1
    return M_r

#-Leaf drip momentum (kg2/s3)
def MomentumLeafDrip(self, pcr, d_l, X, rho, DRAINA, g, C_c):
    a_2 = pcr.scalar(pcr.pcrand(d_l >= 0.033, X >= 7.5)) * 5.14
    a_2 = a_2 + pcr.scalar(pcr.pcrand(d_l >= 0.033, X < 7.5)) * 1.93
    a_2 = a_2 + pcr.scalar(d_l < 0.033) * 0

    b_2 = pcr.scalar(pcr.pcrand(d_l >= 0.033, X >= 7.5)) * 660
    b_2 = b_2 + pcr.scalar(pcr.pcrand(d_l >= 0.033, X < 7.5)) * 1640
    b_2 = b_2 + pcr.scalar(d_l < 0.033) * 2200

    M_beta = a_2 + b_2 * d_l
    V_d = (M_beta * g * (1 - pcr.exp(-2*X/M_beta)))**0.5

    L_d = C_c

    M_d = self.pi / 6 * V_d**2 * rho**2 * d_l**3 * L_d * DRAINA
    return M_d

#-Detachment of soil particles by runoff (kg/m2/s)
def DetachmentRunoff(self, pcr, k_f, C_r, tau, tau_cr, C_g):
    D_q = pcr.ifthenelse(tau > tau_cr, k_f * pcr.max(0, 1 - C_g - C_r - self.NoErosion_SHETRAN) * (tau/tau_cr - 1), 0)
    
    #-set values in channels to 0 in case channels should be excluded
    if self.exclChannelsFLAG == 1:
        D_q = D_q * self.Hillslope
    return D_q

#-Water depth (m) and flow width (m)
def Manning(self, pcr, Q, n, WD_ratio, S):
    A = ((Q * n * (2 * ((WD_ratio**2 + 1) / WD_ratio)**0.5)**(2/3)) / S**0.5)**(3/4)
    h = (A / WD_ratio)**0.5
    l = WD_ratio * h
    return h, l

#-Shear stress (N/m2)
def ShearStress(self, pcr, rho, g, h, S):
    tau = rho * g * h * S
    return tau

#-Critical shear stress (N/m2)
def ShearStressCritical(self, pcr, tau, rho_s, rho, g, D_50, nu):
    R_star = pcr.max(0.03, (D_50 * (tau / rho)**0.5) / nu)

    a_3 = pcr.scalar(R_star > 400) * 0.056
    a_3 = a_3 + pcr.scalar(pcr.pcrand(R_star > 135, R_star <= 400)) * 0.03
    a_3 = a_3 + pcr.scalar(pcr.pcrand(R_star > 30, R_star <= 135)) * 0.013
    a_3 = a_3 + pcr.scalar(pcr.pcrand(R_star > 6, R_star <= 30)) * 0.033
    a_3 = a_3 + pcr.scalar(pcr.pcrand(R_star > 1, R_star <= 6)) * 0.1
    a_3 = a_3 + pcr.scalar(pcr.pcrand(R_star >= 0.03, R_star <= 1)) * 0.1

    b_3 = pcr.scalar(R_star > 400) * 0
    b_3 = b_3 + pcr.scalar(pcr.pcrand(R_star > 135, R_star <= 400)) * 0.1
    b_3 = b_3 + pcr.scalar(pcr.pcrand(R_star > 30, R_star <= 135)) * 0.28
    b_3 = b_3 + pcr.scalar(pcr.pcrand(R_star > 6, R_star <= 30)) * 0
    b_3 = b_3 + pcr.scalar(pcr.pcrand(R_star > 1, R_star <= 6)) * -0.62
    b_3 = b_3 + pcr.scalar(pcr.pcrand(R_star >= 0.03, R_star <= 1)) * -0.3

    tau_cr = (rho_s - rho) * g * D_50 * a_3 * R_star**b_3
    return tau_cr


#-Capacity particulate transport rate for overland flow (m3/s)
def Capacity(self, pcr, np, tau, tau_cr, rho, rho_s, g, h, l, Q, D50, S):
    if self.CapacityEquation == 1: # Yalin (1963)
        a = 2.45 * (tau_cr / ((rho_s - rho) * g * D50))**0.5 * (rho_s / rho)**(-0.4)
        delta = pcr.max(tau / tau_cr - 1, 0)

        G_tot = 0.635 * (tau / rho)**0.5 * l * D50 * delta * (1 - 1 / (a * pcr.max(1e-6, delta)) * pcr.ln(1 + a * delta))
    elif self.CapacityEquation == 2: # Engelund & Hansen (1963)
        G_tot = pcr.ifthenelse(h > 0, (0.05 * Q**2 * S**(3/2)) / ((g * h)**0.5 * (rho_s / rho - 1)**2 * D50 * l), 0)
    return G_tot

#-Mass of sediment transported (kg/m2/s)
def sedimentTransported(self, pcr, D_r, D_q, G_tot, rho_s):
    #-determine transport capacity in (kg/m2/s)
    TC = G_tot * rho_s / pcr.cellarea()
    sed = pcr.ifthenelse(D_r + D_q < TC, D_r + D_q, TC)
    return sed

#-init processes shetran
def init(self, pcr, config):
    #-read table with SHETRAN landuse specific model parameters
    pcr.setglobaloption('matrixtable')
    shetran_table = self.inpath + config.get('SHETRAN', 'shetran_table')
    self.d_l_SHETRAN = pcr.lookupscalar(shetran_table, 1, self.LandUse)
    self.X_SHETRAN = pcr.lookupscalar(shetran_table, 2, self.LandUse)
    self.C_g_SHETRAN = pcr.lookupscalar(shetran_table, 3, self.LandUse)
    self.C_c_table_SHETRAN = pcr.lookupscalar(shetran_table, 4, self.LandUse)
    self.n_table_SHETRAN = pcr.lookupscalar(shetran_table, 5, self.LandUse)
    self.NoErosion_SHETRAN = pcr.lookupscalar(shetran_table, 6, self.LandUse)
    pcr.setglobaloption('columntable')

    #-read other model parameters
    self.WD_ratio_SHETRAN = config.getfloat('SHETRAN', 'WD_ratio')
    self.rho_SHETRAN = config.getfloat('SHETRAN', 'rho')
    self.rho_s_SHETRAN = config.getfloat('SHETRAN', 'rho_s')
    self.deltaClay_SHETRAN = config.getfloat('SHETRAN', 'deltaClay') * 1e-6
    self.deltaSilt_SHETRAN = config.getfloat('SHETRAN', 'deltaSilt') * 1e-6
    self.deltaSand_SHETRAN = config.getfloat('SHETRAN', 'deltaSand') * 1e-6
    self.CapacityEquation = config.getint('SHETRAN', 'capacityEquation')
    self.k_r_SHETRAN = config.getfloat('SHETRAN', 'k_r')
    self.k_f_SHETRAN = config.getfloat('SHETRAN', 'k_f')

    #-define some constants
    self.F_w_SHETRAN = 1
    self.g_SHETRAN = 9.81
    self.nu_SHETRAN = 1e-06

    #-determine median grain size
    if self.PedotransferFLAG == 1:
        self.D50_SHETRAN = pcr.ifthenelse(self.RootClayMap > 0.5, pcr.scalar(self.deltaClay_SHETRAN), 0)
        self.D50_SHETRAN = pcr.ifthenelse(self.RootClayMap + self.RootSiltMap > 0.5, self.deltaClay_SHETRAN + (self.deltaSilt_SHETRAN - self.deltaClay_SHETRAN) * (0.5 - self.RootClayMap) / self.RootSiltMap, 0) + self.D50_SHETRAN
        self.D50_SHETRAN = pcr.ifthenelse(self.RootClayMap + self.RootSiltMap < 0.5, self.deltaSilt_SHETRAN + (self.deltaSand_SHETRAN - self.deltaSilt_SHETRAN) * (self.RootSandMap - 0.5) / self.RootSandMap, 0) + self.D50_SHETRAN
    else:
        self.D50_SHETRAN = config.getfloat('SHETRAN', 'D50') * 1e-6


#-dynamic processes shetran
def dynamic(self, pcr, np, Precip, Q):
    #-determine canopy cover from LAI
    if self.DynVegFLAG == 1:
        C_c_SHETRAN = pcr.min(1, self.LAI)
    else:
        C_c_SHETRAN = self.C_c_table_SHETRAN

    #-determine rainfall drop momentum (kg2/s3)
    M_r = self.shetran.MomentumRainDrop(self, pcr, np, C_c_SHETRAN, Precip)

    #-determine water drainage rate from canopy (m/s)
    DRAINA = Precip * 1e-3 / (24 * 60 * 60)

    #-determine leaf drop momentum (kg2/s3)
    M_d = self.shetran.MomentumLeafDrip(self, pcr, self.d_l_SHETRAN, self.X_SHETRAN, self.rho_SHETRAN, DRAINA, self.g_SHETRAN, C_c_SHETRAN)

    #-determine detachment of soil particles by raindrop impact (kg/m2/s)
    D_r = self.shetran.DetachmentRaindrop(self, pcr, self.k_r_SHETRAN, self.F_w_SHETRAN, self.C_g_SHETRAN, self.RockFrac, M_r, M_d)

    #-report detachment of soil particles by raindrop impact (ton / cell)
    self.reporting.reporting(self, pcr, 'DetRn', D_r * pcr.cellarea() * 1e-3 * (24 * 60 * 60))
 
    #-determine water depth (m) and width of the flow (m)
    h, l = self.shetran.Manning(self, pcr, Q, self.n_table_SHETRAN, self.WD_ratio_SHETRAN, self.Slope)

    #-determine shear stress (N/m2)
    tau = self.shetran.ShearStress(self, pcr, self.rho_SHETRAN, self.g_SHETRAN, h, self.Slope)

    #-determine shear stress (N/m2)
    tau_cr = self.shetran.ShearStressCritical(self, pcr, tau, self.rho_s_SHETRAN, self.rho_SHETRAN, self.g_SHETRAN, self.D50_SHETRAN, self.nu_SHETRAN)

    #-determine detachment of soil particles by runoff (kg/m2/s)
    D_q = self.shetran.DetachmentRunoff(self, pcr, self.k_f_SHETRAN, self.RockFrac, tau, tau_cr, self.C_g_SHETRAN)

    #-report detachment of soil particles by runoff (ton / cell)
    self.reporting.reporting(self, pcr, 'DetRun', D_q * pcr.cellarea() * 1e-3 * (24 * 60 * 60))

    #-Determine transport capacity of the flow (m3/s)
    G_tot = self.shetran.Capacity(self, pcr, np, tau, tau_cr, self.rho_SHETRAN, self.rho_s_SHETRAN, self.g_SHETRAN, h, l, Q, self.D50_SHETRAN, self.Slope)

    #-determine mass of sediment in transport (kg/m2/s)
    sed = self.shetran.sedimentTransported(self, pcr, D_r, D_q, G_tot, self.rho_s_SHETRAN)

    #-report sediment in transport (ton / cell)
    self.reporting.reporting(self, pcr, 'SedTrans', sed * pcr.cellarea() * 1e-3 * (24 * 60 * 60))

    return sed
