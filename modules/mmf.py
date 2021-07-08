# Soil erosion module using the Morgan-Morgan-Finney soil erosion model
# Copyright (C) 2017-2019 Joris Eekhout / Spanish National Research Council (CEBAS-CSIC)
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

print('MMF module imported')

#-Effective rainfall (Rf, mm)
def RainEff(self, pcr, precip):
    Rf = precip * pcr.cos(self.Slope)
    return Rf

#-Leaf drainage (LD; mm)
def LeafDrain(self, pcr, Rf):
    LD = Rf * self.CC
    return LD

#-Direct throughfall (DT; mm)
def DirectThroughfall(self, pcr, Rf, LD):
    DT = Rf - LD
    return DT

#-Kinetic energy of the direct throughfall (KE_DT; J/m2)
def KineticEnergyDT(self, pcr, DT, PrecInt):
    # KE_DT = DT * (8.95 + 8.44 * pcr.log10(PrecInt)) # Marshall and Palmer (1948)
    # KE_DT = DT * (0.283 * (1 - 0.52 * pcr.exp(-0.042 * PrecInt))) * 100 # Van Dijk et al (2002)
    # KE_DT = DT * (0.384 * (1 - 0.54 * pcr.exp(-0.029 * PrecInt))) * 100 # Cerro et al (1998)
    # KE_DT = DT * (0.1418 * PrecInt**0.172) * 100 # Meshesha et al. (2016)
    KE_DT = DT * (0.29 * (1 - 0.72 * pcr.exp(-0.05 * PrecInt))) * 100 # Brown and Foster (1987)
    return KE_DT

#-Kinetic energy of the leaf drainage (KE_LD; J/m2)
def KineticEnergyLD(self, pcr, LD, PlantHeight):
    KE_LD = pcr.ifthenelse(PlantHeight < 0.15, 0, LD * (15.8 * PlantHeight**0.5 - 5.87))
    return KE_LD

#-Detachment of soil particles by raindrop impact (F; kg/m2)
def DetachmentRaindrop(self, pcr, K, texture, KE):
    F = K * texture * pcr.max(0, 1 - (self.NoErosion + self.Cover)) * KE * 1e-3
    return F

#-Detachment of soil particles by runoff (H; kg/m2)
def DetachmentRunoff(self, pcr, DR, texture, Q):
    H = DR * texture * Q**1.5 * pcr.max(0, 1 - (self.NoErosion + self.Cover)) * pcr.sin(self.Slope)**(0.3) * 1e-3

    #-set values in channels to 0 in case channels should be excluded
    if self.exclChannelsFLAG == 1:
        H = H * self.Hillslope

    return H

#-Manning for tilled conditions (manningTilled; s/m1/3)
def ManningTillage(self, pcr):
    manningTilled = pcr.exp(-2.1132 + 0.0349 * self.RFR)
    return manningTilled

#-Manning for vegetated conditions (manningVegetated; s/m1/3)
def manningVegetation(self, pcr, waterDepth, diameter, noElements):
    manningVegetated = (waterDepth**(0.67)) / ((2 * 9.81) / (diameter * noElements))**0.5
    return manningVegetated

#-Flow velocity based on manning (v; m/s)
def FlowVelocity(self, pcr, manning, waterDepth):
    v = 1 / manning * waterDepth**(0.67) * pcr.ifthenelse(self.Slope == 0, 1e-5, self.Slope)**0.5
    return v

#-Particle fall number (N_f; -)
def ParticleFallNumber(self, pcr, delta, velocity, waterDepth):
    v_s = (float(1)/18 * (delta ** 2) * (self.rho_s - self.rho) * 9.81) / (self.eta)
    N_f = (pcr.celllength() / pcr.cos(self.Slope) * v_s) / (velocity * waterDepth)
    return N_f

#-Percentage of the detached sediment that is deposited (DEP; %)
def Deposition(self, pcr, N_f):
    DEP = pcr.min(100, 44.1 * N_f ** 0.29)
    return DEP

#-Delivery of detached particles to runoff (G; kg/m2)
def MaterialTransport(self, pcr, F, H, DEP):
    G = (F + H) * (1 - DEP / 100)
    D = (F + H) * DEP / 100
    return G, D

#-Transport capacity of the runoff (TC; ton/ha)
def TransportCapacity(self, pcr, roughnessFactor, texture, Q):
    q = (Q/1000) * pcr.celllength() #-runoff discharge per unit width in m2/day
    TC = roughnessFactor * q**self.TC_beta * self.SlopeStreams**self.TC_gamma #-determine transport capacity
    return TC

#-init processes mmf
def init(self, pcr, config):
    #-if pedotransfer functions are not used read the sand and clay maps
    if self.PedotransferFLAG == 0:
        self.RootSandMap = pcr.readmap(self.inpath + config.get('PEDOTRANSFER','RootSandMap')) / 100
        self.RootClayMap = pcr.readmap(self.inpath + config.get('PEDOTRANSFER','RootClayMap')) / 100
        self.RootSiltMap = 1 - self.RootSandMap - self.RootClayMap

    #-read precipitation intensity map or float
    try:
        self.PrecInt = pcr.readmap(self.inpath + config.get('MMF', 'PrecInt'))
    except:
        self.PrecInt = config.getfloat('MMF', 'PrecInt')

    #-read flag if the canopy cover should be determined based on LAI
    self.CanopyCoverLAIFlag = config.getfloat('MMF', 'CanopyCoverLAIFlag')

    #-read table with MMF input parameters per landuse class
    pcr.setglobaloption('matrixtable')
    MMF_table = self.inpath + config.get('MMF', 'MMF_table')
    self.PlantHeight = pcr.lookupscalar(MMF_table, 1, self.LandUse)
    self.NoElements = pcr.lookupscalar(MMF_table, 2, self.LandUse)
    self.Diameter = pcr.lookupscalar(MMF_table, 3, self.LandUse)
    self.CC_table = pcr.lookupscalar(MMF_table, 4, self.LandUse)
    self.GC_table = pcr.lookupscalar(MMF_table, 5, self.LandUse)
    self.NoErosion = pcr.lookupscalar(MMF_table, 6, self.LandUse)
    self.Tillage = pcr.lookupscalar(MMF_table, 7, self.LandUse)
    self.n_table = pcr.lookupscalar(MMF_table, 8, self.LandUse)
    self.NoVegetation = pcr.lookupscalar(MMF_table, 9, self.LandUse)
    pcr.setglobaloption('columntable')

    #-nominal map with reservoir IDs and extent
    if self.ResFLAG == 1:
        if self.ETOpenWaterFLAG == 1:
            self.Reservoirs = pcr.ifthenelse(pcr.scalar(self.openWaterNominal) > 0, pcr.scalar(1), pcr.scalar(0))
            self.Reservoirs = pcr.cover(self.Reservoirs, 0)
        else:
            self.Reservoirs = pcr.readmap(self.inpath + config.get('RESERVOIR', 'reservoirs'))
        self.NoErosion = pcr.min(self.NoErosion + pcr.scalar(self.Reservoirs), 1)

    #-read table with MMF input parameters per landuse class for the period after harvest
    self.harvest_FLAG = config.getfloat('MMF', 'harvestFLAG')
    pcr.setglobaloption('matrixtable')
    if self.harvest_FLAG:
        MMF_harvest_table = self.inpath + config.get('MMF', 'MMF_harvest')
        self.Sowing = pcr.lookupscalar(MMF_harvest_table, 1, self.LandUse)
        self.Harvest = pcr.lookupscalar(MMF_harvest_table, 2, self.LandUse)
        self.PlantHeight_harvest = pcr.lookupscalar(MMF_harvest_table, 3, self.LandUse)
        self.NoElements_harvest = pcr.lookupscalar(MMF_harvest_table, 4, self.LandUse)
        self.Diameter_harvest = pcr.lookupscalar(MMF_harvest_table, 5, self.LandUse)
        self.CC_harvest = pcr.lookupscalar(MMF_harvest_table, 6, self.LandUse)
        self.GC_harvest = pcr.lookupscalar(MMF_harvest_table, 7, self.LandUse)
        self.Tillage_harvest = pcr.lookupscalar(MMF_harvest_table, 8, self.LandUse)
    pcr.setglobaloption('columntable')

    #-read other model parameters
    self.K_c = config.getfloat('MMF', 'K_c')
    self.K_z = config.getfloat('MMF', 'K_z')
    self.K_s = config.getfloat('MMF', 'K_s')
    self.DR_c = config.getfloat('MMF', 'DR_c')
    self.DR_z = config.getfloat('MMF', 'DR_z')
    self.DR_s = config.getfloat('MMF', 'DR_s')
    self.deltaClay = config.getfloat('MMF', 'deltaClay')
    self.deltaSilt = config.getfloat('MMF', 'deltaSilt')
    self.deltaSand = config.getfloat('MMF', 'deltaSand')
    self.n_bare = config.getfloat('MMF', 'manning')
    self.d_bare = config.getfloat('MMF', 'depthBare')
    self.d_field = config.getfloat('MMF', 'depthInField')
    self.d_TC = config.getfloat('MMF', 'depthTC')
    self.RFR = config.getfloat('MMF', 'RFR')
    self.rho_s = config.getfloat('MMF', 'rho_s')
    self.rho = config.getfloat('MMF', 'rho')
    self.eta = config.getfloat('MMF', 'eta')

    #-Determine manning for soil and tilled conditions
    self.n_tilled = self.mmf.ManningTillage(self, pcr)
    self.n_soil = pcr.ifthenelse(self.Tillage == 1, self.n_tilled, self.n_bare)

    #-Determine flow velocity for in field deposition
    self.n_veg_field = self.mmf.manningVegetation(self, pcr, self.d_field, self.Diameter, self.NoElements)
    self.n_veg_field = pcr.ifthenelse(self.NoVegetation == 1, 0, self.n_veg_field)
    self.n_veg_field = pcr.ifthenelse(self.NoErosion == 1, 0, self.n_veg_field)
    self.n_veg_field = pcr.ifthenelse(self.n_table > 0, self.n_table, self.n_veg_field)
    self.n_field = (self.n_soil**2 + self.n_veg_field**2)**0.5
    self.v_field = self.mmf.FlowVelocity(self, pcr, self.n_field, self.d_field)

    #-Determine flow velocity after harvest
    if self.harvest_FLAG:
        self.n_veg_field_harvest = self.mmf.manningVegetation(self, pcr, self.d_field, self.Diameter_harvest, self.NoElements_harvest)
        self.n_veg_field_harvest = pcr.ifthenelse(self.Tillage_harvest == 1, 0, self.n_veg_field_harvest)
        self.n_field_harvest = (self.n_soil**2 + self.n_veg_field_harvest**2)**0.5
        self.v_field_harvest = self.mmf.FlowVelocity(self, pcr, self.n_field_harvest, self.d_field)


#-dynamic processes
def dynamic(self, pcr, Precip, Runoff):
    #-determine canopy cover from LAI
    if self.CanopyCoverLAIFlag == 1 and self.DynVegFLAG == 1:
        self.CC = pcr.min(1, self.LAI)
    else:
        self.CC = self.CC_table

    #-determine areas that have been harvested
    if self.harvest_FLAG:
        self.Harvested = self.Slope * 0
        self.Harvested = pcr.ifthenelse(self.Harvest < self.Sowing, pcr.ifthenelse(pcr.pcrand(self.Harvest < self.curdate.timetuple().tm_yday, self.Sowing > self.curdate.timetuple().tm_yday), 1, self.Harvested), self.Harvested)
        self.Harvested = pcr.ifthenelse(self.Harvest > self.Sowing, pcr.ifthenelse(pcr.pcror(self.curdate.timetuple().tm_yday > self.Harvest, self.curdate.timetuple().tm_yday < self.Sowing), 1, self.Harvested), self.Harvested)
        self.Harvested = pcr.ifthenelse(self.Harvest == 0, 0, self.Harvested)
    
    #-set canopy cover to value from MMF harvest table for months between harvest and sowing
    if self.CanopyCoverLAIFlag == 0 and self.harvest_FLAG:
        self.CC = pcr.ifthenelse(self.Harvested == 1, self.CC_harvest, self.CC_table)

    #-set ground cover to value from MMF harvest table for months between harvest and sowing
    if self.harvest_FLAG:
        self.GC = pcr.ifthenelse(self.Harvested == 1, self.GC_harvest, self.GC_table)
    else:
        self.GC = self.GC_table
    
    #-define cover as  fraction of soil covered by ground cover and rock
    if self.SnowFLAG == 1:
        SCover = pcr.scalar(self.TotalSnowStore > 0)
        self.Cover = pcr.min(SCover + self.GC + self.RockFrac, 1)
    else:
        self.Cover = pcr.min(self.GC + self.RockFrac, 1)

    #-determine effective rainfall
    Rf = self.mmf.RainEff(self, pcr, Precip)

    #-determine leaf drainage
    LD = self.mmf.LeafDrain(self, pcr, Rf)

    #-determine direct throughfall
    DT = self.mmf.DirectThroughfall(self, pcr, Rf, LD)
    
    #-obtain precipitation intensity from direct throughfall and fraction of the rain in the highest intensity
    if self.InfilFLAG == 1:
        self.PrecInt = DT * self.Alpha
    
    #-determine kinetic energy of the direct throughfall
    KE_DT = self.mmf.KineticEnergyDT(self, pcr, DT, self.PrecInt)
    KE_DT = pcr.ifthenelse(DT == 0, 0, KE_DT)

    #-update plant height for months between harvest and sowing
    if self.harvest_FLAG:
        self.PlantHeightUpdate = pcr.ifthenelse(self.Harvested == 1, self.PlantHeight_harvest, self.PlantHeight)
    else:
        self.PlantHeightUpdate = self.PlantHeight
    
    #-determine kinetic energy of the leaf drainage
    KE_LD = self.mmf.KineticEnergyLD(self, pcr, LD, self.PlantHeightUpdate)

    #-determine total kinetic energy
    KE = KE_DT + KE_LD

    #-determine detachment of soil particles by raindrop impact
    F_c = self.mmf.DetachmentRaindrop(self, pcr, self.K_c, self.RootClayMap, KE)
    F_z = self.mmf.DetachmentRaindrop(self, pcr, self.K_z, self.RootSiltMap, KE)
    F_s = self.mmf.DetachmentRaindrop(self, pcr, self.K_s, self.RootSandMap, KE)
    F = F_c + F_z + F_s

    #-report detachment of soil particles by raindrop impact (ton / cell)
    self.reporting.reporting(self, pcr, 'DetRn', F * pcr.cellarea() / 1000)

    #-determine detachment of soil particles by runoff
    H_c = self.mmf.DetachmentRunoff(self, pcr, self.DR_c, self.RootClayMap, Runoff)
    H_z = self.mmf.DetachmentRunoff(self, pcr, self.DR_z, self.RootSiltMap, Runoff)
    H_s = self.mmf.DetachmentRunoff(self, pcr, self.DR_s, self.RootSandMap, Runoff)
    H = H_c + H_z + H_s

    #-report detachment of soil particles by runoff (ton / cell)
    self.reporting.reporting(self, pcr, 'DetRun', H * pcr.cellarea() / 1000)

    #-replace velocity for vegetated conditions for tilled soil conditions in case of harvested areas
    if self.harvest_FLAG:
        self.v_update = pcr.ifthenelse(self.Harvested == 1, self.v_field_harvest, self.v_field)
    else:
        self.v_update = self.v_field

    #-determine particle fall number
    N_f_c = self.mmf.ParticleFallNumber(self, pcr, self.deltaClay, self.v_update, self.d_field)
    N_f_z = self.mmf.ParticleFallNumber(self, pcr, self.deltaSilt, self.v_update, self.d_field)
    N_f_s = self.mmf.ParticleFallNumber(self, pcr, self.deltaSand, self.v_update, self.d_field)

    #-determine percentage of the detached sediment that is deposited within the cell of origin
    DEP_c = self.mmf.Deposition(self, pcr, N_f_c)
    DEP_z = self.mmf.Deposition(self, pcr, N_f_z)
    DEP_s = self.mmf.Deposition(self, pcr, N_f_s)

    #-determine delivery of detached particles to runoff and sediment that is deposited within the cell of origin
    tempvar = self.mmf.MaterialTransport(self, pcr, F_c, H_c, DEP_c)
    G_c = tempvar[0]
    D_c = tempvar[1]
    tempvar = self.mmf.MaterialTransport(self, pcr, F_z, H_z, DEP_z)
    G_z = tempvar[0]
    D_z = tempvar[1]
    tempvar = self.mmf.MaterialTransport(self, pcr, F_s, H_s, DEP_s)
    G_s = tempvar[0]
    D_s = tempvar[1]
    G = G_c + G_z + G_s
    D = D_c + D_z + D_s

    #-report in field deposition of detached particles (ton / cell)
    self.reporting.reporting(self, pcr, 'SDepFld', D * pcr.cellarea() / 1000)

    #-report sediment in transport (ton / cell)
    self.reporting.reporting(self, pcr, 'SedTrans', G * pcr.cellarea() / 1000)

    return G