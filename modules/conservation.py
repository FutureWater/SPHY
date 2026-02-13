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

import numpy as np

print('conservation module imported')


#-pedotransfer conservation
def pedotransfer(self, pcr, config):
    #-read change in organic matter map and multiply with rootzone OM map
    self.input.input(self, config, pcr, 'changeOM', 'CONSERVATION', 'changeOM', 0)
    self.RootOMMap = self.RootOMMap * (1 + self.changeOM / 100)

    #-read change in bulk density map and multiply with rootzone BD map
    self.input.input(self, config, pcr, 'changeBD', 'CONSERVATION', 'changeBD', 0)
    self.RootBulkMap = self.RootBulkMap * (1 + self.changeBD / 100)


# #-structural measures
# def structural(self, pcr, config):
#     #-read change in organic matter map and multiply with rootzone OM map
#     self.input.input(self, config, pcr, 'structuralMap', 'CONSERVATION', 'structural', 0)
#     self.structuralMap = pcr.cover(self.structuralMap, 0)
#     self.ReInfiltrationFLAG = config.getint('CONSERVATION', 'ReInfiltrationFLAG')
#     self.ReInfil_b = config.getfloat('CONSERVATION', 'ReInfil_b')

#     #-read table with conservation input parameters per conservation measure class
#     pcr.setglobaloption('matrixtable')
#     structural_table = self.inpath + config.get('CONSERVATION', 'structural_table')
#     self.structuralType = pcr.lookupscalar(structural_table, 1, self.structuralMap)
#     self.pondDepth = pcr.lookupscalar(structural_table, 2, self.structuralMap)
#     self.pondArea = pcr.lookupscalar(structural_table, 3, self.structuralMap)
#     self.NoElements_conservation = pcr.lookupscalar(structural_table, 4, self.structuralMap)
#     self.Diameter_conservation = pcr.lookupscalar(structural_table, 5, self.structuralMap)
#     self.n_table_conservation = pcr.lookupscalar(structural_table, 6, self.structuralMap)
#     pcr.setglobaloption('columntable')

#     self.PondsFLAG = np.any(np.unique(pcr.pcr2numpy(self.structuralMap, -9999)) == 1)
#     self.BufferFLAG = np.any(np.unique(pcr.pcr2numpy(self.structuralMap, -9999)) == 2)

#-setup ponds using the reservoir module
def ponds_init(self, pcr, config):
    #-read change in organic matter map and multiply with rootzone OM map
    # self.input.input(self, config, pcr, 'ponds', 'CONSERVATION', 'ponds', 0)
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
    pcr.setglobaloption('columntable')

    if self.ResFLAG == 0:
        #-turn on reservoir module
        self.ResFLAG = 1

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
    else:
        #-define reservoir parameters for ponds
        self.ResSimple = True
        self.ResKr = pcr.ifthenelse(self.ponds, self.pondKr, self.ResKr)
        self.ResB = pcr.ifthenelse(self.ponds, self.ones * 1.5, self.ResB)
        self.ResSmax = pcr.ifthenelse(self.ponds, self.pondDepth * self.pondArea, self.ResSmax)
        self.ResID = pcr.nominal(pcr.ifthenelse(pcr.cover(self.ponds, 0), np.max(pcr.pcr2numpy(self.ResID, -9999)) + pcr.scalar(pcr.cover(self.pondsID, 0)), pcr.scalar(self.ResID)))
        self.ResFunc = pcr.ifthenelse(self.ponds, 1, self.ResFunc)
        self.StorRES = pcr.ifthenelse(self.ponds, self.ResSmax * 0.5, self.StorRES)
        self.QFRAC = pcr.ifthenelse(self.ponds, pcr.scalar(0), self.QFRAC)

#-Determine reinfiltration in ponds
def ponds_reinfiltration(self, pcr, config):
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


#-input parameters for cover crops
def cover_crops_init(self, pcr, config):
    #-in case cover crops are applied
    if self.coverCropsFLAG == 1:
        #-read cover crops map
        self.CoverCrops = pcr.readmap(self.inpath + config.get('CONSERVATION', 'coverCrops'))
    
        #-read table with cover crops input parameters per cover crop class
        pcr.setglobaloption('matrixtable')
        cover_crops_table = self.inpath + config.get('CONSERVATION', 'coverCrops_table')
        self.Sowing_CC = pcr.lookupscalar(cover_crops_table, 1, self.CoverCrops)
        self.Harvest_CC = pcr.lookupscalar(cover_crops_table, 2, self.CoverCrops)
        self.PlantHeight_CC = pcr.lookupscalar(cover_crops_table, 3, self.CoverCrops)
        self.NoElements_CC = pcr.lookupscalar(cover_crops_table, 4, self.CoverCrops)
        self.Diameter_CC = pcr.lookupscalar(cover_crops_table, 5, self.CoverCrops)
        self.GC_CC = pcr.lookupscalar(cover_crops_table, 6, self.CoverCrops)
        pcr.setglobaloption('columntable')

        #-Determine manning for in field deposition
        manningHillslopeVegetation = self.mmf.manningVegetation(self.d_field, self.Diameter_CC, self.NoElements_CC)
        self.n_field_CC = (self.n_soil**2 + manningHillslopeVegetation**2)**0.5

        #-Determine flow velocity for in field deposition
        self.v_field_CC = self.mmf.FlowVelocity(self, pcr, self.n_field_CC, self.d_field)


#-dynamic processes for cover crops
def cover_crops_dynamic(self, pcr):
    #-determine areas that have been harvested
    self.Harvested_CC = self.ones * 0
    self.Harvested_CC = pcr.ifthenelse(self.Harvest_CC < self.Sowing_CC, pcr.ifthenelse(pcr.pcrand(self.Harvest_CC < self.curdate.timetuple().tm_yday, self.Sowing_CC > self.curdate.timetuple().tm_yday), 1, self.Harvested_CC), self.Harvested_CC)
    self.Harvested_CC = pcr.ifthenelse(self.Harvest_CC > self.Sowing_CC, pcr.ifthenelse(pcr.pcror(self.curdate.timetuple().tm_yday > self.Harvest_CC, self.curdate.timetuple().tm_yday < self.Sowing_CC), 1, self.Harvested_CC), self.Harvested_CC)
    self.Harvested_CC = pcr.ifthenelse(self.Harvest_CC == 0, 0, self.Harvested_CC)
    self.Harvested_CC = pcr.cover(self.Harvested_CC, 0)
    
    #-set ground cover to cover crop value for months between sowing and harvest of cover crops
    self.GC = pcr.ifthenelse(pcr.pcrand(self.CoverCrops > 0, self.Harvested_CC == 0), self.GC_CC, self.GC)

    #-set plant height to cover crop value for months between sowing and harvest of cover crops
    self.PlantHeightUpdate = pcr.ifthenelse(pcr.pcrand(self.CoverCrops > 0, self.Harvested_CC == 0), self.PlantHeight_CC, self.PlantHeightUpdate)

    #-set flow velocity to cover crop value for months between sowing and harvest of cover crops
    self.v_update = pcr.ifthenelse(pcr.pcrand(self.CoverCrops > 0, self.Harvested_CC == 0), self.v_field_CC, self.v_update)


# #-sediment transport conservation
# def sediment_transport(self, pcr, config):
#     #-read conservation measures map
#     self.input.input(self, config, pcr, 'conservationMeasures', 'CONSERVATION', 'conservationMeasures', 0)

#     #-read table with conservation input parameters per conservation measure class
#     pcr.setglobaloption('matrixtable')
#     CONSERVATION_table = self.inpath + config.get('CONSERVATION', 'CONSERVATION_table')
#     self.NoElements_conservation = pcr.lookupscalar(CONSERVATION_table, 1, self.conservationMeasures)
#     self.Diameter_conservation = pcr.lookupscalar(CONSERVATION_table, 2, self.conservationMeasures)
#     self.n_table_conservation = pcr.lookupscalar(CONSERVATION_table, 3, self.conservationMeasures)
#     pcr.setglobaloption('columntable')

#     #-Determine flow velocity for conservation measures
#     self.n_veg_TC_conservation = self.roughness.manningVegetation(self.d_field, self.Diameter_conservation, self.NoElements_conservation)
#     self.n_veg_TC_conservation = pcr.ifthenelse(self.n_table_conservation > 0, self.n_table_conservation, self.n_veg_TC_conservation)
#     self.n_TC_conservation = (self.n_soil**2 + self.n_veg_TC_conservation**2)**0.5
#     self.v_TC_conservation = self.mmf.FlowVelocity(self, pcr, self.n_TC_conservation, self.d_TC)
