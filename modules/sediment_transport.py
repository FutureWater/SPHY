# Sediment transport module that determines sediment flux and sediment yield at the stations,
# considering reservoir sedimentation when the reservoir module is used
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
    # return TCSubcatchment

#-Sediment transport
def SedTrans(self, pcr, np, sed, TC):
    #-determine sediment transport without reservoirs
    if self.ResFLAG == 0:
        #-rout sediment based on transport capacity
        sedimentFlux = pcr.accucapacityflux(self.FlowDir, sed, TC)
        sedimentYield = pcr.accucapacitystate(self.FlowDir, sed, TC)
        sedDep = sedimentYield

    #-determine sediment transport with reservoirs
    else:
        #-store sed in sedTrans to be used in routing algorithm
        sedTrans = sed

        #-initiate empty map to be used in for-loop
        sedimentYield = self.DEM * 0
        sedimentFlux = self.DEM * 0
        subFinished = self.DEM * 0
        sedDep = self.DEM * 0

        #-increase the transport capacity in the reservoir cells such that all sediment is transported towards the end of the reservoir
        TC = pcr.ifthenelse(pcr.scalar(self.SedReservoirs) > 0, 1e10, TC)

        #-Determine total soil erosion
        sedTotal = np.sum(pcr.pcr2numpy(sed, 0))

        #-Only rout sediment when soil erosion > 0
        if (sedTotal > 0):
            #-loop through the catchments, rout sediment and determine sedimentation in reservoirs based on trapping efficiency
            for step in range(max(self.subcatchmentSteps)+1): #-for-loop through the steps
                #-determine the subcatchments for this step
                subcatchments = np.where(self.subcatchmentSteps == step)[0]

                #-rout sediment based on transport capacity
                sedTransCapFlux = pcr.accucapacityflux(self.FlowDir, sedTrans, TC)
                sedTransCapState = pcr.accucapacitystate(self.FlowDir, sedTrans, TC)

                #-initiate empty map to be used in for-loop
                stepBool = self.DEM * 0

                #-for-loop through the subcatchments per step
                for subcatchment in subcatchments:
                    #-create boolean map with location of reservoir
                    reservoirBool = pcr.scalar(pcr.ifthenelse(self.ResSedID == int(self.subcatchmentOrder[subcatchment]), pcr.scalar(1), pcr.scalar(0)))

                    #-extract routed sediment value at the reservoir from sedTransCapFlux
                    reservoirFluxTC = pcr.ifthen(reservoirBool == 1, sedTransCapFlux)

                    #-store trapped sediment in sedimentYield (multiply routed sediment value with trapping efficiency to be stored in reservoir cell)
                    sedimentYield = pcr.ifthenelse(reservoirBool == 1, reservoirFluxTC * self.TrappingEff, sedimentYield)

                    #-update subFinished and give subcatchment cells value 1
                    subFinished = pcr.ifthenelse(pcr.scalar(self.subcatchmentRes) == int(self.subcatchmentOrder[subcatchment]), pcr.scalar(1), subFinished)

                    #-update sedTrans, set all subcatchment cells to 0
                    sedTrans = pcr.ifthenelse(pcr.scalar(subFinished) == 1, 0, sedTrans)
                    #-add reservoir outflow to cell downstream of reservoir (multiply routed sediment value with outflow efficiency)
                    sedTrans = sedTrans + pcr.upstream(self.FlowDir, pcr.ifthenelse(reservoirBool == 1, reservoirFluxTC * self.OutflowEff, pcr.scalar(0)))

                    #-create boolean map with location of reservoir
                    stepBool = stepBool + pcr.scalar(pcr.ifthenelse(self.subcatchmentRes == int(self.subcatchmentOrder[subcatchment]), self.subcatchmentRes == int(self.subcatchmentOrder[subcatchment]), pcr.boolean(0)))

                # store sedTransCapFlux in sedimentFlux
                sedimentFlux = sedTransCapFlux * subFinished + sedimentFlux

                #-rout sediment based on transport capacity
                sedDep = sedDep + sedTransCapState * stepBool

            #-rout sediment based on transport capacity
            sedTransCapFlux = pcr.accucapacityflux(self.FlowDir, sedTrans, TC)

            # store sedTransCapFlux in sedimentFlux
            sedimentFlux = sedTransCapFlux * (1 - subFinished) + sedimentFlux
    
    return sedimentYield, sedDep, sedimentFlux


#-init processes
def init(self, pcr, config, csv, np):
    #-init processes when reservoir module is used
    if self.ResFLAG == 1:
        #-nominal map with reservoir IDs
        self.ResSedID = pcr.cover(self.ResID, 0)
    else:
        self.ResSedID = self.Locations

    # #-determine upstream area map
    # self.UpstreamArea = pcr.accuflux(self.FlowDir, 1) * pcr.cellarea() / 10**6

    # #-Read input parameters
    # self.rho = config.getfloat('SEDIMENT_TRANS', 'rho')
    # self.rho_s = config.getfloat('SEDIMENT_TRANS', 'rho_s')
    # self.deltaClay = config.getfloat('SEDIMENT_TRANS', 'deltaClay') * 1e-6
    # self.deltaSilt = config.getfloat('SEDIMENT_TRANS', 'deltaSilt') * 1e-6
    # self.deltaSand = config.getfloat('SEDIMENT_TRANS', 'deltaSand') * 1e-6
    self.SedTransEquation = config.getint('SEDIMENT_TRANS', 'SedTransEquation')
    self.SedTransEquationRills = config.getint('SEDIMENT_TRANS', 'SedTransEquationRills')

    self.Y_cr_Yalin = config.getfloat('SEDIMENT_TRANS', 'Y_cr_Yalin')
    self.tau_c_Wong = config.getfloat('SEDIMENT_TRANS', 'tau_c_Wong')
    self.omega_cr_Govers = config.getfloat('SEDIMENT_TRANS', 'omega_cr_Govers')
    self.Y_cr_Abrahams = config.getfloat('SEDIMENT_TRANS', 'Y_cr_Abrahams')

    self.SedConcMax = config.getfloat('SEDIMENT_TRANS', 'SedConcMax')

    if self.SedTransEquation == 6:
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
    # pcr.report(self.D50, self.outpath + "D50.map")
    
    #-init processes when reservoir module is used
    if self.ResFLAG == 1:
        #-nominal map with reservoir IDs and extent
        if self.ETOpenWaterFLAG == 1:
            self.SedReservoirs = pcr.cover(self.openWaterNominal, 0)
        else:
            self.SedReservoirs = pcr.readmap(self.inpath + config.get('RESERVOIR', 'reservoirs'))
        self.SedReservoirs = pcr.cover(self.SedReservoirs, 0)

        #-read table with the trapping efficiency per reservoir
        self.TrapEffTab = self.inpath + config.get('SEDIMENT_TRANS', 'TrapEffTab')
        self.TrappingEff = pcr.cover(pcr.lookupscalar(self.TrapEffTab, self.ResSedID), 0)

        #-construct map where all cells have 1 and only the reservoir cells have trapping efficiency value obtained from the table
        self.OutflowEff = pcr.cover(1-pcr.lookupscalar(self.TrapEffTab, self.ResSedID), 1)

        #-determine subcatchment map
        self.subcatchmentRes = pcr.subcatchment(self.FlowDir, self.ResSedID)

        #-read reservoir order for sediment transport and add the values to self.subcatchmentOrder and self.subcatchmentSteps
        self.ResOrder = config.get('SEDIMENT_TRANS', 'ResOrder')
        self.subcatchmentOrder = []
        self.subcatchmentSteps = []

        #-loop through the rows of the text file
        with open(self.inpath + self.ResOrder, 'rt') as f:
            next(f) # skip headings
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                self.subcatchmentOrder = np.append(self.subcatchmentOrder, int(row[0]))
                self.subcatchmentOrder = self.subcatchmentOrder.astype(np.int)
                self.subcatchmentSteps = np.append(self.subcatchmentSteps, int(row[1]))
                self.subcatchmentSteps = self.subcatchmentSteps.astype(np.int)
        
        #-loop through the steps and define map with step per subcatchment
        self.subcatchmentStepsMap = self.DEM * 0
        for step in range(max(self.subcatchmentSteps)+1):
            #-determine the subcatchments for this step
            subcatchments = np.where(self.subcatchmentSteps == step)[0]

            #-for-loop through the subcatchments per step
            for subcatchment in subcatchments:
                #-update subFinished and give subcatchment cells value 1
                self.subcatchmentStepsMap = pcr.ifthenelse(pcr.scalar(self.subcatchmentRes) == int(self.subcatchmentOrder[subcatchment]), pcr.scalar(step + 1), self.subcatchmentStepsMap)
        #-set all subcatchments with value 0 to the maximum step 
        self.subcatchmentStepsMap = pcr.ifthenelse(self.subcatchmentStepsMap == 0, int(max(self.subcatchmentSteps)+1), self.subcatchmentStepsMap)
        #-map with step value for stations that are not reservoirs
        self.LocationsNoResSteps = pcr.scalar(self.LocationsNoRes) * self.subcatchmentStepsMap

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

    #-import conservation module
    import modules.conservation
    self.conservation = modules.conservation
    del modules.conservation

    #-read init processes sediment transport
    self.conservation.init(self, pcr, config)

# #-initial conditions sediment transport
# def initial(self, pcr, config):
#     try:
#         self.SYieldR = pcr.readmap(self.inpath + config.get('SEDIMENT_TRANS', 'Sed_init'))
#     except:
#         self.SYieldR = config.getfloat('SEDIMENT_TRANS', 'Sed_init')

# #-dynamic sediment transport processes musle
# def dynamic_musle(self, pcr):
#     #-transport capacity
#     TC = self.sediment_transport.TC(self, pcr, (Q * 3600 * 24) / pcr.cellarea() * 1000)

#     #-report the transport capacity per subcatchment
#     self.reporting.reporting(self, pcr, 'TC', TC)

#     #-determine sediment yield at reservoirs
#     tempvar = self.sediment_transport.SedTrans(self, pcr, np, sed, TC)
#     sedimentYield = tempvar[0]
#     sedDep = tempvar[1]

#     #-report the deposition in channel cells
#     self.reporting.reporting(self, pcr, 'SedDep', sedDep)

#     #-report sediment yield in the reservoirs
#     self.reporting.reporting(self, pcr, 'SYieldRA', sedimentYield)

# #-dynamic sediment transport processes mmf
# def dynamic_mmf(self, pcr, Runoff, np, G):
#     #-change the flow factor for harvested areas to actual and tillage conditions
#     if self.harvest_FLAG == 1:
#         self.roughnessFactorUpdate = pcr.ifthenelse(self.Harvested == 1, self.v_TC_harvest / self.v_b, self.roughnessFactor)
#     else:
#         self.roughnessFactorUpdate = self.roughnessFactor
    
#     #-determine transport capacity
#     TC = self.mmf.TransportCapacity(self, pcr, self.roughnessFactorUpdate, self.RootClayMap + self.RootSiltMap + self.RootSandMap, Runoff)

#     #-report the transport capacity
#     self.reporting.reporting(self, pcr, 'TC', TC)

#     #-determine sediment yield at stations
#     sedYield, sedDep, sedFlux = self.sediment_transport.SedTrans(self, pcr, np, G * pcr.cellarea() / 1000, TC)

#     #-report the sediment deposition by transport capacity (ton/day)
#     self.reporting.reporting(self, pcr, 'SedDep', sedDep)

#     #-report sediment yield in the stations (ton/day)
#     self.reporting.reporting(self, pcr, 'SedYld', sedYield)

#     #-report sediment flux in the stations (ton/day)
#     self.reporting.reporting(self, pcr, 'SedFlux', sedFlux)

# #-Shear stress (N/m2)
# def ShearStress(rho, g, h, S):
#     tau = rho * g * h * S
#     return tau

# #-Critical shear stress (N/m2)
# def ShearStressCritical(pcr, tau, rho_s, rho, g, D_50, nu):
#     R_star = pcr.max(0.03, (D_50 * (tau / rho)**0.5) / nu)

#     a_3 = pcr.scalar(R_star > 400) * 0.056
#     a_3 = a_3 + pcr.scalar(pcr.pcrand(R_star > 135, R_star <= 400)) * 0.03
#     a_3 = a_3 + pcr.scalar(pcr.pcrand(R_star > 30, R_star <= 135)) * 0.013
#     a_3 = a_3 + pcr.scalar(pcr.pcrand(R_star > 6, R_star <= 30)) * 0.033
#     a_3 = a_3 + pcr.scalar(pcr.pcrand(R_star > 1, R_star <= 6)) * 0.1
#     a_3 = a_3 + pcr.scalar(pcr.pcrand(R_star >= 0.03, R_star <= 1)) * 0.1

#     b_3 = pcr.scalar(R_star > 400) * 0
#     b_3 = b_3 + pcr.scalar(pcr.pcrand(R_star > 135, R_star <= 400)) * 0.1
#     b_3 = b_3 + pcr.scalar(pcr.pcrand(R_star > 30, R_star <= 135)) * 0.28
#     b_3 = b_3 + pcr.scalar(pcr.pcrand(R_star > 6, R_star <= 30)) * 0
#     b_3 = b_3 + pcr.scalar(pcr.pcrand(R_star > 1, R_star <= 6)) * -0.62
#     b_3 = b_3 + pcr.scalar(pcr.pcrand(R_star >= 0.03, R_star <= 1)) * -0.3

#     tau_cr = (rho_s - rho) * g * D_50 * a_3 * R_star**b_3
#     return tau_cr

#-Determine transport capacity (g/l)
def Capacity(self, pcr, rho, rho_s, g, h, w, Q, D50, S, SedTransEquation):
    if SedTransEquation == 1: # Govers (1990)
        c = ((D50 * 1e6 + 5) / 0.32)**(-0.6)
        d = ((D50 * 1e6 + 5) / 300)**(0.25)
        omega = self.flowVelocity * 100 * self.slopeChannel
        # omega_cr = 0.4
        omega_cr = self.omega_cr_Govers
        TC_f = c * pcr.max(0, omega - omega_cr)**d * self.rho_s
        TC = pcr.max(TC_f / (1 - TC_f/self.rho_s), 0)
    elif SedTransEquation == 2: # Yalin (1963)
        R = self.hydraulicRadius
        s = rho_s / rho
        Y = (R * S) / (D50 * (s-1))
        U_star = pcr.sqrt(g * R * S)
        # Y_cr = 0.06
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
        # Y_cr = 0.06
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
        # tau_c = 0.0495
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
        TC = self.mmf.TransportCapacity(self, pcr, self.roughnessFactorUpdate, self.RootClayMap + self.RootSiltMap + self.RootSandMap, Runoff)

        #-report the transport capacity
        self.reporting.reporting(self, pcr, 'TC', TC)

        #-determine sediment yield at stations
        sedYield, sedDep, sedFlux = self.sediment_transport.SedTrans(self, pcr, np, Sed, TC)

        #-report the sediment deposition by transport capacity (ton/day)
        self.reporting.reporting(self, pcr, 'SedDep', sedDep)

        #-report sediment yield in the stations (ton/day)
        self.reporting.reporting(self, pcr, 'SedYld', sedYield)

        #-report sediment flux in the stations (ton/day)
        self.reporting.reporting(self, pcr, 'SedFlux', sedFlux)

    else:
        #-determine shear stress (N/m2)
        if self.travelTimeFLAG == 0:
        #     tau = self.sediment_transport.ShearStress(self.rho, self.g, self.waterDepth, self.Slope)
        # else:
            #-determine water depth (m) and flow depth (m)
            h, l = self.shetran.Manning(self, pcr, Q, self.n_table_SHETRAN, self.WD_ratio_SHETRAN, self.Slope)

        #     #-determine shear stress (N/m2)
        #     tau = self.shetran.ShearStress(self, pcr, self.rho, self.g, h, self.Slope)

        # pcr.report(self.waterDepth, self.outpath + "waterDepth_" + str(self.counter).zfill(3) + ".map")
        # pcr.report(tau, self.outpath + "tau_" + str(self.counter).zfill(3) + ".map")

        #-For loop over the sediment classes
        for sedimentClass in self.sedimentClasses:
            #-Define median grain size for sediment class
            D50 = getattr(self, "delta" + sedimentClass)

            #-Determine transport capacity of the flow (g/l)
            if self.travelTimeFLAG == 1:
                TC = self.sediment_transport.Capacity(self, pcr, self.rho, self.rho_s, self.g, self.waterDepth, self.channelWidth, Q, D50, self.slopeChannel, self.SedTransEquation)
            else:
                TC = self.sediment_transport.Capacity(self, pcr, self.rho, self.rho_s, self.g, h, l, Q, D50, self.Slope, self.SedTransEquation)
            # pcr.report(TC, self.outpath + "TC_channel_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")

            #-Determine transport capacity in the rills when rills are simulated
            if self.RillFLAG == 1:
                #-Determine transport capacity of the flow (g/l)
                if self.travelTimeFLAG == 1:
                    TC_Rills = self.sediment_transport.Capacity(self, pcr, self.rho, self.rho_s, self.g, self.waterDepth, self.channelWidth, Q, D50, self.slopeChannel, self.SedTransEquationRills)
                else:
                    TC_Rills = self.sediment_transport.Capacity(self, pcr, self.rho, self.rho_s, self.g, h, l, Q, D50, self.Slope, self.SedTransEquationRills)

                #-Update transport capacity for the hillslopes
                TC = pcr.ifthenelse(self.channelHillslope == 2, TC_Rills, TC)
                # pcr.report(TC_Rills, self.outpath + "TC_hillslope_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")

            #-apply maximum allowable sediment concentration (g/l = kg/m3)
            TC = pcr.min(TC, self.SedConcMax)
            # pcr.report(TC, self.outpath + "TC_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")

            #-determine transport capacity (ton/day)
            # pcr.report(TC, self.outpath + "sedConc_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")
            TC = TC * Q * 1e-3 * (24 * 60 * 60)
            # pcr.report(TC, self.outpath + "TC_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")

            #-store transport capacity for sediment class
            setattr(self, "TC" + sedimentClass, TC)

        #     #-When MMF is used, obtain hillslope erosion per sediment class from attribute
        #     if sedimentClass == 'Gravel':
        #         Sed = self.ones * 0
        #     else:
        #         if self.ErosionModel == 2:
        #             Sed = getattr(self, "Sed" + sedimentClass)
        #         else:
        #             #-Multiply hillslope erosion by sediment class fraction
        #             Sed = Sed * getattr(self, "Root" + sedimentClass + "Map")
        #     pcr.report(Sed, self.outpath + "Sed_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")

        #     #-determine sediment yield at stations
        #     sedYield, sedDep, sedFlux = self.sediment_transport.SedTrans(self, pcr, np, Sed, TC)
        #     pcr.report(sedFlux, self.outpath + "sedFlux_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")
        #     pcr.report(sedDep, self.outpath + "sedDep_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")

        #     #-Assign sediment yield, deposition and flux values to sediment class
        #     setattr(self, "sedYield" + sedimentClass, sedYield)
        #     setattr(self, "sedDep" + sedimentClass, sedDep)
        #     setattr(self, "sedFlux" + sedimentClass, sedFlux)

        #-report the transport capacity (ton/day)
        TC = self.TCClay + self.TCSilt + self.TCSand + self.TCGravel
        self.reporting.reporting(self, pcr, 'TC', TC)
        # pcr.report(TC, self.outpath + "TC_" + str(self.counter).zfill(3) + ".map")

        # #-report sediment yield in the stations (ton/day)
        # sedYield = self.sedYieldClay + self.sedYieldSilt + self.sedYieldSand
        # self.reporting.reporting(self, pcr, 'SedYld', sedYield)

        # #-report sediment flux in the stations (ton/day)
        # sedFlux = self.sedFluxClay + self.sedFluxSilt + self.sedFluxSand
        # self.reporting.reporting(self, pcr, 'SedFlux', sedFlux)
        # # pcr.report(sedFlux, self.outpath + "sedFlux_" + str(self.counter).zfill(3) + ".map")

        # return TC