# Sediment transport module that determines sediment flux, which can be used as 
# input for the morphodynamics module
# Copyright (C) 2015-2019 Joris Eekhout / Spanish National Research Council (CEBAS-CSIC)
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

#- Equations for simulating sediment transport
print('Sediment transport module imported')

#-Sediment transport capacity
def TC(self, pcr, runoff):
    #-runoff discharge per unit width in m2/day
    q = (runoff/1000) * pcr.celllength()
    #-transport capacity calculation ton/day
    TC = self.TC_k * q**self.TC_beta * self.Slope**self.TC_gamma
    #-only apply TC to cells with runoff (>0) and in the channel (upstream area > upstream_km2)
    TC_upstream = pcr.ifthen(pcr.pcrand(self.Channel == 1, runoff > 0), TC)
    #-determine average transport capacity per subcatchment
    TCSubcatchment = pcr.areaaverage(TC_upstream, self.subcatchmentRes)
    #-set all other cells to 0
    TCSubcatchment = pcr.cover(TCSubcatchment, 0)
    return TC

#-init processes
def init(self, pcr, config, csv, np):
    #-init processes when reservoir module is used
    if self.ResFLAG == 1:
        #-nominal map with reservoir IDs
        self.ResSedID = pcr.cover(self.ResID, 0)
    else:
        self.ResSedID = self.Locations

    # #-Read input parameters
    self.SedTransEquation = config.getint('SEDIMENT_TRANS', 'SedTransEquation')
    self.SedTransEquationRills = config.getint('SEDIMENT_TRANS', 'SedTransEquationRills')
    self.Y_cr_Yalin = config.getfloat('SEDIMENT_TRANS', 'Y_cr_Yalin')
    self.tau_c_Wong = config.getfloat('SEDIMENT_TRANS', 'tau_c_Wong')
    self.omega_cr_Govers = config.getfloat('SEDIMENT_TRANS', 'omega_cr_Govers')
    self.Y_cr_Abrahams = config.getfloat('SEDIMENT_TRANS', 'Y_cr_Abrahams')
    self.SedConcMax = config.getfloat('SEDIMENT_TRANS', 'SedConcMax')

    #-set sediment equation to 6 in case travel time algorithm is not used
    if self.travelTimeFLAG == 0:
        self.SedTransEquation = 6

    if self.SedTransEquation == 6:
        #-determine upstream area map
        self.UpstreamArea = pcr.accuflux(self.FlowDir, 1) * pcr.cellarea() / 10**6

        #-determine upstream area smaller than upstream_km2 and define channel cells based on upstream area
        self.Upstream_km2 = config.getfloat('SEDIMENT_TRANS', 'upstream_km2')
        self.Channel = self.UpstreamArea > self.Upstream_km2

        #-determine average slope per stream based on stream order and averaged per subcatchment
        self.Basin = pcr.subcatchment(self.FlowDir, self.ResSedID)
        self.StreamOrder = pcr.streamorder(self.FlowDir)
        self.Streams = pcr.scalar(self.Channel) * pcr.scalar(self.Basin) * 10 + pcr.scalar(self.Channel) * pcr.scalar(self.StreamOrder)
        self.SlopeStreams = pcr.areaaverage(self.Slope, pcr.nominal(self.Streams))
        self.SlopeStreams = pcr.ifthenelse(self.Channel == 1, self.SlopeStreams, self.Slope)

        #-read transport capacity parameters
        self.TC_beta = config.getfloat('SEDIMENT_TRANS', 'TC_beta')
        self.TC_gamma = config.getfloat('SEDIMENT_TRANS', 'TC_gamma')
        self.d_TC = config.getfloat('SEDIMENT_TRANS', 'depthTC')

    #-Define sediment size classes
    self.sedimentClasses = ['Clay', 'Silt', 'Sand', 'Gravel']

    #-define some constants
    self.g = 9.81
    self.nu = 1e-06

    #-determine median grain size
    if self.PedotransferFLAG == 1:
        self.D50 = pcr.ifthenelse(self.RootClayMap > 0.5, pcr.scalar(self.deltaClay), 0)
        self.D50 = pcr.ifthenelse(self.RootClayMap + self.RootSiltMap > 0.5, self.deltaClay + (self.deltaSilt - self.deltaClay) * (0.5 - self.RootClayMap) / self.RootSiltMap, 0) + self.D50
        self.D50 = pcr.ifthenelse(self.RootClayMap + self.RootSiltMap < 0.5, self.deltaSilt + (self.deltaSand - self.deltaSilt) * (self.RootSandMap - 0.5) / self.RootSandMap, 0) + self.D50
    else:
        self.D50 = config.getfloat('SEDIMENT_TRANS', 'D50') * 1e-6

    #-determine roughness factor for mmf only
    if self.SedTransEquation == 6:
        #-Read flag if channels should be excluded from the detachment by runoff calculation
        self.manningChannelsFLAG = config.getint('SEDIMENT_TRANS', 'manningChannelFLAG')

        #-read manning value for channels
        if self.manningChannelsFLAG == 1:
            self.input.input(self, config, pcr, 'manningChannel', 'ROUTING', 'channelManning', 0)

        #-Determine flow velocity for transport capacity calculation
        self.n_veg_TC = self.roughness.manningVegetation(self.d_TC, self.Diameter, self.NoElements)
        self.n_veg_TC = pcr.ifthenelse(self.NoVegetation == 1, 0, self.n_veg_TC)
        self.n_veg_TC = pcr.ifthenelse(self.NoErosion == 1, 0, self.n_veg_TC)
        self.n_veg_TC = pcr.ifthenelse(self.n_table > 0, self.n_table, self.n_veg_TC)
        self.n_TC = (self.n_soil**2 + self.n_veg_TC**2)**0.5

        #-set manning value of channels to predetermined value
        if self.manningChannelsFLAG == 1:
            self.n_TC = pcr.ifthenelse(self.Channel == 1, self.manningChannel, self.n_TC)
        self.v_TC = self.mmf.FlowVelocity(self, pcr, self.n_TC, self.d_TC)

        #-Determine flow velocity after harvest, manning for tilled conditions is used
        if self.harvest_FLAG:
            self.n_veg_TC_harvest = self.roughness.manningVegetation(self.d_field, self.Diameter_harvest, self.NoElements_harvest)
            self.n_veg_TC_harvest = pcr.ifthenelse(self.Tillage_harvest == 1, 0, self.n_field_harvest)
            self.n_TC_harvest = (self.n_soil**2 + self.n_veg_TC_harvest**2)**0.5

            #-set manning value of channels to predetermined value
            if self.manningChannelsFLAG == 1:
                self.n_TC_harvest = pcr.ifthenelse(pcr.pcrand(self.Channel == 1, self.n_TC_harvest > 0), self.manningChannel, self.n_TC_harvest)
            self.v_TC_harvest = self.mmf.FlowVelocity(self, pcr, self.n_TC_harvest, self.d_TC)

        #-Determine flow velocity for bare soil conditions (reference conditions)
        self.v_b = self.mmf.FlowVelocity(self, pcr, self.n_bare, self.d_bare)

        #-Determine roughness factor for transport capacity calculation
        self.roughnessFactor = self.v_TC / self.v_b

        #-read WD ratio for water depth and flow velocity calculation
        self.WD_ratio_SHETRAN = config.getfloat('SHETRAN', 'WD_ratio')

    #-import conservation module
    import modules.conservation
    self.conservation = modules.conservation
    del modules.conservation

    #-read init processes sediment transport
    self.conservation.init(self, pcr, config)

#-Determine transport capacity (g/l)
def Capacity(self, pcr, rho, rho_s, g, h, w, Q, D50, S, SedTransEquation):
    if SedTransEquation == 1: # Govers (1990)
        c = ((D50 * 1e6 + 5) / 0.32)**(-0.6)
        d = ((D50 * 1e6 + 5) / 300)**(0.25)
        omega = self.flowVelocity * 100 * self.slopeChannel
        omega_cr = self.omega_cr_Govers
        TC_f = c * pcr.max(0, omega - omega_cr)**d * self.rho_s
        TC = pcr.max(TC_f / (1 - TC_f/self.rho_s), 0)
    elif SedTransEquation == 2: # Yalin (1963)
        R = self.hydraulicRadius
        s = rho_s / rho
        Y = (R * S) / (D50 * (s-1))
        U_star = pcr.sqrt(g * R * S)
        Y_cr = self.Y_cr_Yalin
        ar = 2.45 / S**0.4 * pcr.sqrt(pcr.max(Y_cr * (Y / Y_cr - 1), 1e-10))
        P = 0.635 * (Y / Y_cr - 1) * (1 - (pcr.ln(1 + ar)) / ar)
        qs = (rho_s - rho) * D50 * U_star * P
        TC = qs * w / pcr.max(Q, 1e-10)
    elif SedTransEquation == 3: # Abrahams et al. (2001)
        R = self.hydraulicRadius
        s = rho_s / rho
        Y = pcr.max((R * S) / (D50 * (s-1)), 1e-10)
        U_star = pcr.max(pcr.sqrt(g * R * S), 1e-10)
        Y_cr = self.Y_cr_Abrahams
        Cr = 0
        Dr = 0.1
        a = 10**(-0.42 * (Cr/Dr)**0.2)
        b = 3.4
        c = 1 + 0.42 * (Cr/Dr)**0.2
        w_i = (g * (s - 1) * D50)**0.5
        d = pcr.ifthenelse((w_i / U_star) > 3, pcr.scalar(0), -0.5)
        qb = a * D50 * U_star * Y * pcr.max(1 - Y_cr / Y, 0)**b * (self.flowVelocity / U_star)**c * (w_i / U_star)**d
        TC = (qb * w * rho_s) / pcr.max(Q, 1e-10)
    elif SedTransEquation == 4: # Wong and Parker (2006)
        R = (rho_s / rho) - 1
        tau_star = (S * h) / (R * D50)
        alpha = 3.97
        beta = 1.5
        tau_c = self.tau_c_Wong
        q_star = alpha * pcr.max(tau_star - tau_c, 0)**beta
        qb = q_star * (R * g * D50)**0.5 * D50
        TC = (qb * w * rho_s) / pcr.max(Q, 1e-10)
    elif SedTransEquation == 5: # Wilcock and Crowe (2003)
        s = (rho_s / rho)
        tau_rm_star = 0.021 + 0.015 * pcr.exp(-20 * self.RootSandMap)
        tau_rm = tau_rm_star * (s-1) * rho * g * self.D50
        b = 0.67 / (1 + pcr.exp(1.5 - D50 / self.D50))
        tau_r = tau_rm * (D50 / self.D50)**b
        tau = rho * g * h * S
        phi = tau / tau_r
        W = pcr.ifthenelse(phi < 1.35, 0.002 * phi**7.5, 14 * (1 - 0.894/(phi**0.5))**4.5)
        u_star = (tau / rho)**0.5
        qb = (W * u_star**3) / ((s-1) * g)
        TC = (qb * w * rho_s) / pcr.max(Q, 1e-10)
    return TC


#-dynamic sediment transport processes
def dynamic(self, pcr, np, Q, Sed):
    if self.SedTransEquation == 6:
        #-change the flow factor for harvested areas to actual and tillage conditions
        if self.harvest_FLAG == 1:
            self.roughnessFactorUpdate = pcr.ifthenelse(self.Harvested == 1, self.v_TC_harvest / self.v_b, self.roughnessFactor)
        else:
            self.roughnessFactorUpdate = self.roughnessFactor

        #-overwrite roughness factor with the value for structural conservation measures
        if self.conservationFLAG == 1:
            self.roughnessFactorUpdate = pcr.ifthenelse(self.conservationMeasures > 0, self.v_TC_conservation / self.v_b, self.roughnessFactorUpdate)
        
        #-determine runoff
        Runoff = (Q * 3600 * 24) / pcr.cellarea() * 1000

        #-determine transport capacity
        self.TC = self.mmf.TransportCapacity(self, pcr, self.roughnessFactorUpdate, self.RootClayMap + self.RootSiltMap + self.RootSandMap, Runoff)

        #-report the transport capacity
        self.reporting.reporting(self, pcr, 'TC', self.TC)

    else:
        # #-in case travel time is not used
        # if self.travelTimeFLAG == 0:
        #     #-import SHETRAN module
        #     import modules.shetran
        #     self.shetran = modules.shetran
        #     del modules.shetran

        #     #-determine water depth (m) and flow depth (m)
        #     h, l = self.shetran.Manning(self, pcr, Q, self.n_TC, self.WD_ratio_SHETRAN, self.Slope)

        #-For loop over the sediment classes
        for sedimentClass in self.sedimentClasses:
            #-Define median grain size for sediment class
            D50 = getattr(self, "delta" + sedimentClass)

            #-Determine transport capacity of the flow (g/l)
            if self.travelTimeFLAG == 1:
                TC = self.sediment_transport.Capacity(self, pcr, self.rho, self.rho_s, self.g, self.waterDepth, self.channelWidth, Q, D50, self.slopeChannel, self.SedTransEquation)
            else:
                TC = self.sediment_transport.Capacity(self, pcr, self.rho, self.rho_s, self.g, h, l, Q, D50, self.Slope, self.SedTransEquation)

            #-Determine transport capacity in the rills when rills are simulated
            if self.RillFLAG == 1:
                #-Determine transport capacity of the flow (g/l)
                if self.travelTimeFLAG == 1:
                    TC_Rills = self.sediment_transport.Capacity(self, pcr, self.rho, self.rho_s, self.g, self.waterDepth, self.channelWidth, Q, D50, self.slopeChannel, self.SedTransEquationRills)
                else:
                    TC_Rills = self.sediment_transport.Capacity(self, pcr, self.rho, self.rho_s, self.g, h, l, Q, D50, self.Slope, self.SedTransEquationRills)

                #-Update transport capacity for the hillslopes
                TC = pcr.ifthenelse(self.channelHillslope == 2, TC_Rills, TC)

            #-apply maximum allowable sediment concentration (g/l = kg/m3)
            TC = pcr.min(TC, self.SedConcMax)

            #-determine transport capacity (ton/day)
            TC = TC * Q * 1e-3 * (24 * 60 * 60)

            #-store transport capacity for sediment class
            setattr(self, "TC" + sedimentClass, TC)

            #-report TC per sediment class
            self.reporting.reporting(self, pcr, 'TC' + sedimentClass, pcr.scalar(getattr(self, "TC" + sedimentClass)))

        #-sum TC for all fractions
        TC = self.TCClay + self.TCSilt + self.TCSand + self.TCGravel

        #-report the transport capacity (ton/day)
        self.reporting.reporting(self, pcr, 'TC', TC)