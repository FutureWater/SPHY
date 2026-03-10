# Morphodynamics module that determines morphological change in rills and channels and
# reservoir sedimentation when the reservoir module is used.
# Copyright (C) 2021-2023 Joris Eekhout / Spanish National Research Council (CEBAS-CSIC)
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

print('Morphodynamics module imported')


#-init morphodynamics processes
def init(self, pcr, pcrm, config, csv, np):
    if self.travelTimeFLAG == 1:
        #-read input parameters
        self.input.input(self, config, pcr, 'rillThickness', 'MORPHODYNAMICS', 'rillThickness', 0)
        self.input.input(self, config, pcr, 'channelThickness', 'MORPHODYNAMICS', 'channelThickness', 0)
        self.bedThickness = self.channelDepth * pcr.ifthenelse(self.channelHillslope == 1, self.channelThickness, self.rillThickness)

        #-set gravel fraction to 0 for hillslope
        self.RootGravelMap = 0

        #-For loop over the sediment classes
        for sedimentClass in self.sedimentClasses:
            #-read channel material percentage
            self.input.input(self, config, pcr, 'channel' + sedimentClass, 'MORPHODYNAMICS', 'channel' + sedimentClass, 0)
            
            #-assign channel material to channels and use soil texture fraction for hillslopes
            setattr(self, "channel" + sedimentClass, pcr.ifthenelse(self.channelHillslope == 1, getattr(self, "channel" + sedimentClass) / 100, getattr(self, "Root" + sedimentClass + "Map")))

            #-if reservoir is used
            if self.ResFLAG == 1:
                #-set channel depth in reservoirs to 0
                setattr(self, "channel" + sedimentClass, pcr.ifthenelse(pcr.scalar(self.Reservoirs) > 0, 0, getattr(self, "channel" + sedimentClass)))

        #-determine max channel depth
        self.channelDepthMax = self.channelDepth + self.bedThickness

    #-init processes when reservoir module is used
    if self.pondsFLAG == 1:
        #-nominal map with reservoir IDs and extent
        if self.ResFLAG == 1:
            self.sedResId = self.ResID
        else:
            self.sedResId = pcr.readmap(self.inpath + config.get('MORPHODYNAMICS', 'sedRes'))
        self.sedResId = pcr.cover(self.sedResId, 0)

        #-define map where reservoirs are located (=1)
        self.sedRes = pcr.ifthenelse(pcr.scalar(self.sedResId) > 0, pcr.scalar(1), pcr.scalar(0))

        #-read table with the trapping efficiency per reservoir
        self.TrapEffTab = self.inpath + config.get('MORPHODYNAMICS', 'TrapEffTab')
        self.TrappingEff = pcr.cover(pcr.lookupscalar(self.TrapEffTab, self.sedResId), 0)

        #-construct map where all cells have 1 and only the reservoir cells have trapping efficiency value obtained from the table
        self.OutflowEff = pcr.cover(1-pcr.lookupscalar(self.TrapEffTab, self.sedResId), 1)

        #-determine subcatchment map
        self.subcatchmentRes = pcr.subcatchment(self.FlowDir, self.sedResId)

        #-determine steps per reservoir
        self.reservoirStep = pcr.ifthen(self.sedRes == 1, pcr.accuflux(self.FlowDir, self.sedRes) * self.sedRes)

        #-determine unique steps
        self.reservoirStepsArray = np.unique(pcr.pcr2numpy(self.reservoirStep, 1))


#-initial morphodynamics processes
def initial(self, pcr):
    if self.travelTimeFLAG == 1:
        #-For loop over the sediment classes
        for sedimentClass in self.sedimentClasses:
            #-assign bed thickness to all sediment classes
            setattr(self, "sedimentStoreInitial" + sedimentClass, getattr(self, "channel" + sedimentClass) * self.bedThickness * self.channelWidth * pcr.celllength() * self.rho_s * 1e-3)
            setattr(self, "sedimentStore" + sedimentClass, getattr(self,  "sedimentStoreInitial" + sedimentClass))
        
        #-determine overall sediment storage
        self.sedimentStore = self.sedimentStoreClay + self.sedimentStoreSilt + self.sedimentStoreSand + self.sedimentStoreGravel
