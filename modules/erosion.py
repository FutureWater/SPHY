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

print('erosion module imported')


#-Determine the number of rills per meter
def numberOfRills(pcr, Flow, MC, S, RR, Re):
    N = 0.66 + 0.69 * pcr.ln(Flow) + 0.91 * pcr.ln(MC) + 2.04 * pcr.ln(S) - 0.37 * pcr.ln(RR) - 0.37 * pcr.ln(Re)
    return N

#-Determine rill dimensions based on minimum and maximum rill size
def rillDimensions(pcr, self):
    #-Determine maximum accuflux on the hillslopes
    rill = self.channelHillslope == 1
    accufluxMax = pcr.areamaximum(pcr.accuflux(self.FlowDir, 1), rill)
    
    #-Determine fraction of accuflux with respect to maximum accuflux on hillslope
    accufluxFraction = pcr.accuflux(self.FlowDir, 1) / accufluxMax

    #-Determine rill width based on minimum and maximum rill size
    rillWidth = self.minRillWidth + (self.maxRillWidth - self.minRillWidth) * accufluxFraction

    return rillWidth

#-init processes erosion module
def init(self, pcr, config, csv, np):
    #-read soil erosion model selector (1 for MUSLE, 2 for MMF)
    self.ErosionModel = config.getfloat('EROSION', 'ErosionModel')

    #-read rock fraction map
    self.RockFrac = pcr.readmap(self.inpath + config.get('EROSION', 'RockFrac'))

    #-Read flag if channels should be excluded from the detachment by runoff calculation
    self.exclChannelsFLAG = config.getint('EROSION', 'exclChannelsFLAG')
    
    #-determine hillslope map if channels should be excluded
    if self.exclChannelsFLAG == 1:
        #-determine upstream area map
        self.UpstreamArea = pcr.accuflux(self.FlowDir, 1) * pcr.cellarea() / 10**6

        #-determine upstream area larger than upstream_km2 and define hillslope cells based on upstream area
        self.Upstream_km2 = config.getfloat('EROSION', 'upstream_km2')
        self.Hillslope = pcr.scalar(self.UpstreamArea <= self.Upstream_km2)

    #-Read soil erosion parameters
    self.deltaClay = config.getfloat('EROSION', 'deltaClay') * 1e-6
    self.deltaSilt = config.getfloat('EROSION', 'deltaSilt') * 1e-6
    self.deltaSand = config.getfloat('EROSION', 'deltaSand') * 1e-6
    self.deltaGravel = config.getfloat('EROSION', 'deltaGravel') * 1e-6
    self.rho_s = config.getfloat('EROSION', 'rho_s')
    self.rho = config.getfloat('EROSION', 'rho')

    #-read table with soil erosion input parameters per landuse class
    pcr.setglobaloption('matrixtable')
    EROSION_table = self.inpath + config.get('EROSION', 'EROSION_table')
    self.NoElements = pcr.lookupscalar(EROSION_table, 1, self.LandUse)
    self.Diameter = pcr.lookupscalar(EROSION_table, 2, self.LandUse)
    self.NoErosion = pcr.lookupscalar(EROSION_table, 3, self.LandUse)
    self.Tillage = pcr.lookupscalar(EROSION_table, 4, self.LandUse)
    self.n_table = pcr.lookupscalar(EROSION_table, 5, self.LandUse)
    self.NoVegetation = pcr.lookupscalar(EROSION_table, 6, self.LandUse)
    pcr.setglobaloption('columntable')

    #-read table with MMF input parameters per landuse class for the period after harvest
    self.harvest_FLAG = config.getfloat('EROSION', 'harvestFLAG')
    pcr.setglobaloption('matrixtable')
    if self.harvest_FLAG:
        EROSION_harvest_table = self.inpath + config.get('EROSION', 'EROSION_harvest')
        self.Sowing = pcr.lookupscalar(EROSION_harvest_table, 1, self.LandUse)
        self.Harvest = pcr.lookupscalar(EROSION_harvest_table, 2, self.LandUse)
        self.NoElements_harvest = pcr.lookupscalar(EROSION_harvest_table, 3, self.LandUse)
        self.Diameter_harvest = pcr.lookupscalar(EROSION_harvest_table, 4, self.LandUse)
        self.Tillage_harvest = pcr.lookupscalar(EROSION_harvest_table, 5, self.LandUse)
    pcr.setglobaloption('columntable')

    # #-read table with MMF input parameters per landuse class
    # pcr.setglobaloption('matrixtable')
    # EROSION_table = self.inpath + config.get('EROSION', 'EROSION_table')
    # self.PlantHeight = pcr.lookupscalar(EROSION_table, 1, self.LandUse)
    # self.NoElements = pcr.lookupscalar(EROSION_table, 2, self.LandUse)
    # self.Diameter = pcr.lookupscalar(EROSION_table, 3, self.LandUse)
    # self.CC_table = pcr.lookupscalar(EROSION_table, 4, self.LandUse)
    # self.GC_table = pcr.lookupscalar(EROSION_table, 5, self.LandUse)
    # self.NoErosion = pcr.lookupscalar(EROSION_table, 6, self.LandUse)
    # self.Tillage = pcr.lookupscalar(EROSION_table, 7, self.LandUse)
    # self.n_table = pcr.lookupscalar(EROSION_table, 8, self.LandUse)
    # self.NoVegetation = pcr.lookupscalar(EROSION_table, 9, self.LandUse)
    # pcr.setglobaloption('columntable')

    # #-read table with MMF input parameters per landuse class for the period after harvest
    # self.harvest_FLAG = config.getfloat('EROSION', 'harvestFLAG')
    # pcr.setglobaloption('matrixtable')
    # if self.harvest_FLAG:
    #     EROSION_harvest_table = self.inpath + config.get('EROSION', 'EROSION_harvest')
    #     self.Sowing = pcr.lookupscalar(EROSION_harvest_table, 1, self.LandUse)
    #     self.Harvest = pcr.lookupscalar(EROSION_harvest_table, 2, self.LandUse)
    #     self.PlantHeight_harvest = pcr.lookupscalar(EROSION_harvest_table, 3, self.LandUse)
    #     self.NoElements_harvest = pcr.lookupscalar(EROSION_harvest_table, 4, self.LandUse)
    #     self.Diameter_harvest = pcr.lookupscalar(EROSION_harvest_table, 5, self.LandUse)
    #     self.CC_harvest = pcr.lookupscalar(EROSION_harvest_table, 6, self.LandUse)
    #     self.GC_harvest = pcr.lookupscalar(EROSION_harvest_table, 7, self.LandUse)
    #     self.Tillage_harvest = pcr.lookupscalar(EROSION_harvest_table, 8, self.LandUse)
    # pcr.setglobaloption('columntable')

    #-nominal map with reservoir IDs and extent
    if self.ResFLAG == 1:
        if self.ETOpenWaterFLAG == 1:
            self.Reservoirs = pcr.ifthenelse(pcr.scalar(self.openWaterNominal) > 0, pcr.scalar(1), pcr.scalar(0))
            self.Reservoirs = pcr.cover(self.Reservoirs, 0)
        else:
            self.Reservoirs = pcr.readmap(self.inpath + config.get('RESERVOIR', 'reservoirs'))
        self.NoErosion = pcr.min(self.NoErosion + pcr.scalar(self.Reservoirs), 1)

    #-Read rill flag
    self.RillFLAG = config.getint('EROSION', 'RillFLAG')

    #-If rill flag is equal to 1 than calculate random roughness, rill dimensions and adjust floodplain width and channel depth
    if self.RillFLAG == 1 and self.travelTimeFLAG == 1:
        #-Determine random roughness based on RFR (Table 2.2; Morgan, 2005)
        self.randomRoughness = pcr.exp(0.29 + 0.099 * self.Tillage) * 1e-3

        #-Read min and max rill width
        self.minRillWidth = config.getfloat('EROSION', 'minRillWidth')
        self.maxRillWidth = config.getfloat('EROSION', 'maxRillWidth')

        #-Determine rill width based on min and max and set rill depth equal to rill width
        self.rillWidth = self.erosion.rillDimensions(pcr, self)
        self.rillDepth = self.rillWidth

        #-Set floodplain width to cell width for hillslope cells
        self.floodplainWidth = pcr.ifthenelse(self.channelHillslope == 2, pcr.celllength(), self.floodplainWidth)

        #-Set channel depth to rill depth for hillslope cells
        self.channelDepth = pcr.ifthenelse(self.channelHillslope == 2, self.rillDepth, self.channelDepth)

        #-Set channel width to rill width for hillslope cells
        self.channelWidth = pcr.ifthenelse(self.channelHillslope == 2, self.rillWidth, self.channelWidth)

    #-import roughness module
    import modules.roughness
    self.roughness = modules.roughness
    del modules.roughness

    #-read init processes sediment transport
    self.roughness.init(self, pcr, config)

    #-read MUSLE input parameters
    if self.ErosionModel == 1:
        #-import musle module
        import modules.musle
        self.musle = modules.musle
        del modules.musle

        #-read init processes musle
        self.musle.init(self, pcr, config)

    #-read MMF input parameters
    if self.ErosionModel == 2:
        #-import mmf module
        import modules.mmf
        self.mmf = modules.mmf
        del modules.mmf

        #-read init processes mmf
        self.mmf.init(self, pcr, config)

    #-read INCA input parameters
    if self.ErosionModel == 3:
        #-import INCA module
        import modules.inca
        self.inca = modules.inca
        del modules.inca

        #-read init processes INCA
        self.inca.init(self, pcr, config)

    #-read SHETRAN input parameters
    if self.ErosionModel == 4:
        #-import SHETRAN module
        import modules.shetran
        self.shetran = modules.shetran
        del modules.shetran

        #-read init processes SHETRAN
        self.shetran.init(self, pcr, config)

    #-read DHSVM input parameters
    if self.ErosionModel == 5:
        #-import DHSVM module
        import modules.dhsvm
        self.dhsvm = modules.dhsvm
        del modules.dhsvm

        #-read init processes DHSVM
        self.dhsvm.init(self, pcr, config)

    #-read HSPF input parameters
    if self.ErosionModel == 6:
        #-import HSPF module
        import modules.hspf
        self.hspf = modules.hspf
        del modules.hspf

        #-read init processes HSPF
        self.hspf.init(self, pcr, config)

#-dynamic erosion processes
def dynamic(self, pcr, np, Precip, Q_m3, Q_mm):
    #-determine canopy cover from LAI
    if self.DynVegFLAG:
        self.CC = pcr.min(1, self.LAI)
    else:
        self.CC = self.CC_table

    #-determine areas that have been harvested
    if self.harvest_FLAG:
        self.Harvested = self.ones * 0
        self.Harvested = pcr.ifthenelse(self.Harvest < self.Sowing, pcr.ifthenelse(pcr.pcrand(self.Harvest < self.curdate.timetuple().tm_yday, self.Sowing > self.curdate.timetuple().tm_yday), 1, self.Harvested), self.Harvested)
        self.Harvested = pcr.ifthenelse(self.Harvest > self.Sowing, pcr.ifthenelse(pcr.pcror(self.curdate.timetuple().tm_yday > self.Harvest, self.curdate.timetuple().tm_yday < self.Sowing), 1, self.Harvested), self.Harvested)
        self.Harvested = pcr.ifthenelse(self.Harvest == 0, 0, self.Harvested)
    
    #-set canopy cover to value from MMF harvest table for months between harvest and sowing
    if self.DynVegFLAG == 0 and self.harvest_FLAG:
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

    #-update plant height for months between harvest and sowing
    if self.harvest_FLAG:
        self.PlantHeightUpdate = pcr.ifthenelse(self.Harvested == 1, self.PlantHeight_harvest, self.PlantHeight)
    else:
        self.PlantHeightUpdate = self.PlantHeight

    #-MUSLE
    if self.ErosionModel == 1:
        #-read dynamic processes musle
        Sed = self.musle.dynamic(self, pcr, Q_mm)

    #-MMF
    if self.ErosionModel == 2:
        #-determine soil erosion in transport (G)
        Sed = self.mmf.dynamic(self, pcr, Precip, Q_mm)

    #-INCA
    if self.ErosionModel == 3:
        #-determine soil erosion
        Sed = self.inca.dynamic(self, pcr, Precip, Q_m3)

    #-SHETRAN
    if self.ErosionModel == 4:
        #-determine soil erosion
        Sed = self.shetran.dynamic(self, pcr, np, Precip, Q_m3)

    #-DHSVM
    if self.ErosionModel == 5:
        #-determine soil erosion
        Sed = self.dhsvm.dynamic(self, pcr, np, Precip, Q_m3)

    #-HSPF
    if self.ErosionModel == 6:
        #-determine soil erosion
        Sed = self.hspf.dynamic(self, pcr, np, Precip, Q_mm)

    # pcr.report(Sed, self.outpath + "Sed_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(self.flowVelocity, self.outpath + "flowVelocity_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(self.waterDepth, self.outpath + "waterDepth_" + str(self.counter).zfill(3) + ".map")
    
    return Sed