# Conservation module
# Copyright (C) 2017-2025 Joris Eekhout / Spanish National Research Council (CEBAS-CSIC)
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

print('conservation module imported')


#-pedotransfer conservation
def pedotransfer(self, pcr, config):
    #-in of changing organic matter is applied
    if self.changeOMFLAG == 1:
        #-read change in organic matter map and multiply with rootzone OM map
        self.input.input(self, config, pcr, 'changeOM', 'CONSERVATION', 'changeOM', 0)
        self.RootOMMap = self.RootOMMap * (1 + self.changeOM / 100)

    #-in of changing bulk density is applied
    if self.changeBDFLAG == 1:
        #-read change in bulk density map and multiply with rootzone BD map
        self.input.input(self, config, pcr, 'changeBD', 'CONSERVATION', 'changeBD', 0)
        self.RootBulkMap = self.RootBulkMap * (1 + self.changeBD / 100)


#-setup ponds using the reservoir module
def ponds_init(self, pcr, config, np):
    #-read change in organic matter map and multiply with rootzone OM map
    self.input.input(self, config, pcr, 'ponds', 'CONSERVATION', 'ponds', 0)
    self.pondsID = pcr.cover(self.ponds, 0)
    self.ponds = pcr.scalar(self.ponds) > 0
    self.ReInfiltrationFLAG = config.getint('CONSERVATION', 'ReInfiltrationFLAG')
    self.ReInfil_b = config.getfloat('CONSERVATION', 'ReInfil_b')

    #-read table with conservation input parameters per conservation measure class
    pcr.setglobaloption('matrixtable')
    ponds_table = self.inpath + config.get('CONSERVATION', 'ponds_table')
    self.pondDepth = pcr.lookupscalar(ponds_table, 1, self.pondsID)
    self.pondArea = pcr.lookupscalar(ponds_table, 2, self.pondsID)
    self.pondKr = pcr.lookupscalar(ponds_table, 3, self.pondsID)
    self.pondTrappEff = pcr.cover(pcr.lookupscalar(ponds_table, 4, self.pondsID), 0)
    pcr.setglobaloption('columntable')

    if self.ResFLAG == 0:
        #-turn on reservoir module
        self.ResFLAG = 1
        self.PondsAndReservoirs = 0

        #-import reservoirs module
        import modules.reservoirs
        self.reservoirs = modules.reservoirs
        del modules.reservoirs
        
        #-define reservoir parameters for ponds
        self.ResSimple = True
        self.ResAdvanced = False
        self.ResKr = pcr.ifthen(self.ponds, self.pondKr)
        self.ResB = pcr.ifthen(self.ponds, self.ones * 1.5)
        self.ResSmax = pcr.ifthen(self.ponds, self.pondDepth * self.pondArea)
        self.ResID = pcr.nominal(self.pondsID)
        self.ResFunc = pcr.cover(pcr.scalar(self.ponds), 0)
        self.StorRES = self.ResSmax * 0.5
        self.QFRAC = pcr.cover(pcr.ifthen(self.ponds, pcr.scalar(0)), 1)
        self.reservoirTrappEff = self.ones * 0
    else:
        #-define reservoir parameters for ponds
        self.PondsAndReservoirs = 1
        self.ResSimple = True
        self.ResKr = pcr.ifthenelse(self.ponds, self.pondKr, self.ResKr)
        self.ResB = pcr.ifthenelse(self.ponds, self.ones * 1.5, self.ResB)
        self.ResSmax = pcr.ifthenelse(self.ponds, self.pondDepth * self.pondArea, self.ResSmax)
        self.ResID = pcr.nominal(pcr.ifthenelse(pcr.cover(self.ponds, 0), np.max(pcr.pcr2numpy(self.ResID, -9999)) + pcr.scalar(pcr.cover(self.pondsID, 0)), pcr.scalar(self.ResID)))
        self.ResFunc = pcr.ifthenelse(self.ponds, 1, self.ResFunc)

#-setup ponds using the reservoir module
def ponds_initial(self, pcr):
    if self.ResFLAG == 1:
        self.StorRES = pcr.ifthenelse(self.ponds, self.ResSmax * 0.5, self.StorRES)
        self.QFRAC = pcr.ifthenelse(self.ponds, pcr.scalar(0), self.QFRAC)

#-Determine reinfiltration in ponds
def ponds_reinfiltration(self, pcr):
    #-Determine relative saturation
    relSat = pcr.min(pcr.max(self.RootWater / self.RootSat, 0), 1)
    
    #-Determine unsaturated hydraulic conductivity (mm/day)
    RootKUnSat = pcr.max(self.RootKsat * (relSat**self.ReInfil_b), 0)
    
    #-Determine infiltration amount (m3/day), which cannot exceed pond volume
    infilPond = pcr.cover(pcr.min(RootKUnSat * 1e-3 * self.pondArea, self.StorRES), 0)

    #-Update reservoir storage and rootwater for ponds
    self.StorRES = pcr.ifthenelse(self.ponds, self.StorRES - infilPond, self.StorRES)
    self.RootWater += pcr.upstream(self.FlowDir, infilPond / pcr.cellarea() * 1e3)

    return infilPond / pcr.cellarea() * 1e3


#-input parameters for vegetation cover
def vegetation_cover_init(self, pcr, config):
    #-in case vegetation cover are applied
    if self.vegetationCoverFLAG == 1:
        #-read vegetation cover map
        self.VegetationCover = pcr.readmap(self.inpath + config.get('CONSERVATION', 'vegetationCover'))
    
        #-read table with vegetation cover input parameters per vegetation cover class
        pcr.setglobaloption('matrixtable')
        vegetation_cover_table = self.inpath + config.get('CONSERVATION', 'vegetationCover_table')
        self.Sowing_VC = pcr.lookupscalar(vegetation_cover_table, 1, self.VegetationCover)
        self.Harvest_VC = pcr.lookupscalar(vegetation_cover_table, 2, self.VegetationCover)
        self.PlantHeight_VC = pcr.lookupscalar(vegetation_cover_table, 3, self.VegetationCover)
        self.NoElements_VC = pcr.lookupscalar(vegetation_cover_table, 4, self.VegetationCover)
        self.Diameter_VC = pcr.lookupscalar(vegetation_cover_table, 5, self.VegetationCover)
        self.GC_VC = pcr.lookupscalar(vegetation_cover_table, 6, self.VegetationCover)
        self.strip_VC = pcr.cover(pcr.lookupscalar(vegetation_cover_table, 7, self.VegetationCover), 0)
        pcr.setglobaloption('columntable')

        #-Determine manning for in field deposition
        manningHillslopeVegetation = self.mmf.manningVegetation(self.d_field, self.Diameter_VC, self.NoElements_VC)
        self.n_field_VC = (self.n_soil**2 + manningHillslopeVegetation**2)**0.5

        #-Determine flow velocity for in field deposition
        self.v_field_VC = self.mmf.FlowVelocity(self, pcr, self.n_field_VC, self.d_field)

#-dynamic processes for vegetation cover
def vegetation_cover_dynamic_harvested(self, pcr):
    #-determine areas that have been harvested
    self.Harvested_VC = self.ones * 0
    self.Harvested_VC = pcr.ifthenelse(self.Harvest_VC < self.Sowing_VC, pcr.ifthenelse(pcr.pcrand(self.Harvest_VC < self.curdate.timetuple().tm_yday, self.Sowing_VC > self.curdate.timetuple().tm_yday), 1, self.Harvested_VC), self.Harvested_VC)
    self.Harvested_VC = pcr.ifthenelse(self.Harvest_VC > self.Sowing_VC, pcr.ifthenelse(pcr.pcror(self.curdate.timetuple().tm_yday > self.Harvest_VC, self.curdate.timetuple().tm_yday < self.Sowing_VC), 1, self.Harvested_VC), self.Harvested_VC)
    self.Harvested_VC = pcr.ifthenelse(self.Harvest_VC == 0, 0, self.Harvested_VC)
    self.Harvested_VC = pcr.cover(self.Harvested_VC, 0)


#-dynamic processes for vegetation cover
def vegetation_cover_dynamic_mmf(self, pcr):
    #-set ground cover to vegetation cover value for months between sowing and harvest of vegetation cover
    self.GC = pcr.ifthenelse(pcr.pcrand(pcr.pcrand(self.VegetationCover > 0, self.Harvested_VC == 0), self.strip_VC == 0), self.GC_VC, self.GC)

    #-set plant height to vegetation cover value for months between sowing and harvest of vegetation cover
    self.PlantHeightUpdate = pcr.ifthenelse(pcr.pcrand(pcr.pcrand(self.VegetationCover > 0, self.Harvested_VC == 0), self.strip_VC == 0), self.PlantHeight_VC, self.PlantHeightUpdate)

    #-set flow velocity to vegetation cover value for months between sowing and harvest of vegetation cover
    self.v_update = pcr.ifthenelse(pcr.pcrand(pcr.pcrand(self.VegetationCover > 0, self.Harvested_VC == 0), self.strip_VC == 0), self.v_field_VC, self.v_update)

#-input parameters for land use change
def land_use_change_init(self, pcr, np, config):
    #-read land use change map
    self.landUseChange = pcr.readmap(self.inpath + config.get('CONSERVATION', 'LandUseChange'))

    #-store original land use map
    self.LandUseOriginal = self.LandUse

    #-apply land use change to land use map
    self.LandUse = pcr.ifthenelse(self.landUseChange > 0, self.landUseChange, self.LandUseOriginal)

    #-determine the land use change classes
    self.landUseChangeClasses = np.unique(pcr.pcr2numpy(self.landUseChange, np.nan))

    #-remove nan and 0
    self.landUseChangeClasses = self.landUseChangeClasses[~np.isnan(self.landUseChangeClasses)]
    self.landUseChangeClasses = self.landUseChangeClasses[self.landUseChangeClasses != 0] 

    
#-dynamic processes for land use change
def land_use_change_dynamic(self, pcr, map, landuse, landusechange, landusechangeclasses):
    #-determine average map per land use class
    map_avg = pcr.areaaverage(map, pcr.nominal(landuse * pcr.scalar(self.clone)))
    
    #-determine the anomaly from average map per land use class 
    map_anomaly = map - map_avg
    
    #-for-loop over the land use change classes
    for cls in landusechangeclasses:
        #-determine the average map for the class
        map_avg_class = pcr.mapmaximum(pcr.ifthen(self.LandUse == int(cls), map_avg))
        
        #-apply average to changed land use
        map = pcr.ifthenelse(landusechange == int(cls), map_avg_class, map)

        #-add map anomaly to change land use
        map = pcr.ifthenelse(landusechange == int(cls), map + map_anomaly, map)

    return map