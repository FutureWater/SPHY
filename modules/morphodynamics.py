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


#-Sediment transport
def SedTrans(self, pcr, np, sed, TC):
    #-determine sediment transport without reservoirs
    if self.SedReservoirsFLAG == 0:
        #-rout sediment based on transport capacity
        sedimentFlux = pcr.accucapacityflux(self.FlowDir, sed, TC)
        sedDep = pcr.accucapacitystate(self.FlowDir, sed, TC)
        sedimentYield = self.ones * 0

    #-determine sediment transport with reservoirs
    else:
        #-store sed in sedTrans to be used in routing algorithm
        sedTrans = sed

        #-initiate empty map to be used in for-loop
        sedimentYield = self.ones * 0
        sedimentFlux = self.ones * 0
        subFinished = self.ones * 0
        sedDep = self.ones * 0

        #-increase the transport capacity in the reservoir cells such that all sediment is transported towards the end of the reservoir
        if self.ResFLAG == 1:
            TC = pcr.ifthenelse(pcr.scalar(self.sedResId) > 0, 1e10, TC)

        #-Determine total soil erosion
        sedTotal = np.sum(pcr.pcr2numpy(sed, 0))

        #-Only rout sediment when soil erosion > 0
        if (sedTotal > 0):
            #-loop through the catchments, rout sediment and determine sedimentation in reservoirs based on trapping efficiency
            for step in self.reservoirStepsArray: #-for-loop through the steps
                #-determine the reservoirs for this step
                reservoirs = np.unique(pcr.pcr2numpy(pcr.ifthen(self.reservoirStep == int(step), self.sedResId), -9999))[1:]

                #-set TC in finished subcatchments to 0
                TC = pcr.ifthenelse(subFinished == 1, 0, TC)

                #-rout sediment based on transport capacity
                sedTransCapFlux = pcr.accucapacityflux(self.FlowDir, sedTrans, TC)
                sedTransCapState = pcr.accucapacitystate(self.FlowDir, sedTrans, TC)

                #-initiate empty map to be used in for-loop
                stepBool = self.ones * 0

                #-for-loop through the reservoirs per step
                for reservoir in reservoirs:
                    #-create boolean map with location of reservoir
                    reservoirBool = pcr.scalar(pcr.ifthenelse(self.sedResId == int(reservoir), pcr.scalar(1), pcr.scalar(0)))

                    #-extract routed sediment value at the reservoir from sedTransCapFlux
                    reservoirFluxTC = pcr.ifthen(reservoirBool == 1, sedTransCapFlux)

                    #-store trapped sediment in sedimentYield (multiply routed sediment value with trapping efficiency to be stored in reservoir cell)
                    sedimentYield = pcr.ifthenelse(reservoirBool == 1, reservoirFluxTC * self.TrappingEff, sedimentYield)

                    #-update subFinished and give subcatchment cells value 1
                    subFinished = pcr.ifthenelse(pcr.scalar(self.subcatchmentRes) == int(reservoir), pcr.scalar(1), subFinished)

                    #-update sedTrans, set all subcatchment cells to 0
                    sedTrans = pcr.ifthenelse(pcr.scalar(subFinished) == 1, 0, sedTrans)

                    #-add reservoir outflow to cell downstream of reservoir (multiply routed sediment value with outflow efficiency)
                    sedTrans = sedTrans + pcr.upstream(self.FlowDir, pcr.ifthenelse(reservoirBool == 1, reservoirFluxTC * self.OutflowEff, pcr.scalar(0)))

                    #-create boolean map with location of reservoir
                    stepBool = stepBool + pcr.scalar(pcr.ifthenelse(self.subcatchmentRes == int(reservoir), self.subcatchmentRes == int(reservoir), pcr.boolean(0)))

                #-store sedTransCapFlux in sedimentFlux
                sedimentFlux = sedTransCapFlux * subFinished + sedimentFlux

                #-store sedTransCapState in sedDep
                sedDep = sedDep + sedTransCapState * stepBool

            #-rout sediment based on transport capacity
            sedTransCapFlux = pcr.accucapacityflux(self.FlowDir, sedTrans, TC)

            # store sedTransCapFlux in sedimentFlux
            sedimentFlux = sedTransCapFlux * (1 - subFinished) + sedimentFlux
    
    return sedimentYield, sedDep, sedimentFlux


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

    if config.get('MORPHODYNAMICS', 'sedRes') != "" or self.ResFLAG == 1:
        self.SedReservoirsFLAG = 1
    else:
        self.SedReservoirsFLAG = 0

    #-init processes when reservoir module is used
    if self.SedReservoirsFLAG == 1:
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


#-dynamic sediment transport processes
def dynamic(self, pcr, pcrm, np, Q, Sed):
    if self.travelTimeFLAG == 0:
        #-determine sediment yield at stations
        sedYield, sedDep, sedFlux = self.morphodynamics.SedTrans(self, pcr, np, Sed, self.TC)

        #-report the sediment deposition by transport capacity (ton/day)
        self.reporting.reporting(self, pcr, 'SedDep', sedDep)

        #-report sediment yield in the stations (ton/day)
        self.reporting.reporting(self, pcr, 'SedYld', sedYield)

        #-report sediment flux in the stations (ton/day)
        self.reporting.reporting(self, pcr, 'SedFlux', sedFlux)
    else:
        #-For loop over the sediment classes
        for sedimentClass in self.sedimentClasses:

            #-When MMF is used, obtain hillslope erosion per sediment class from attribute
            if sedimentClass == 'Gravel':
                Sed = self.ones * 0
            else:
                if self.ErosionModel == 2:
                    #-Get hillslope erosion by sediment class fraction from MMF
                    Sed = getattr(self, "Sed" + sedimentClass)
                else:
                    #-Multiply hillslope erosion by sediment class fraction
                    Sed = Sed * getattr(self, "Root" + sedimentClass + "Map")

            #-get transport capacity for sediment class
            TC = getattr(self, "TC" + sedimentClass)

            #-determine initial sediment store
            sedimentStoreInitial = getattr(self, "sedimentStoreInitial" + sedimentClass) #-reset channel storage to initial value

            #-determine total amount available for transport
            material = sedimentStoreInitial + Sed #* 0

            #-route sediment through the catchment
            sedYield, sedDep, sedFlux = self.morphodynamics.SedTrans(self, pcr, np, material, TC)

            #-determine channel change
            channelChange = sedDep - sedimentStoreInitial

            #-Assign sediment yield, deposition and flux values to sediment class
            setattr(self, "sedYield" + sedimentClass, sedYield)
            setattr(self, "sedDep" + sedimentClass, sedDep)
            setattr(self, "channelChange" + sedimentClass, channelChange)
            setattr(self, "sedFlux" + sedimentClass, sedFlux)

        #-determine total channel change
        channelChange = self.channelChangeClay + self.channelChangeSilt + self.channelChangeSand + self.channelChangeGravel
        self.reporting.reporting(self, pcr, 'ChChng', channelChange)

        #-report the sediment deposition by transport capacity (ton/day)
        sedDep = self.sedDepClay + self.sedDepSilt + self.sedDepSand + self.sedDepGravel
        self.reporting.reporting(self, pcr, 'SedDep', sedDep)

        #-report sediment yield in the stations (ton/day)
        sedYield = self.sedYieldClay + self.sedYieldSilt + self.sedYieldSand + self.sedYieldGravel
        self.reporting.reporting(self, pcr, 'SedYld', sedYield)

        #-report sediment flux in the stations (ton/day)
        sedFlux = self.sedFluxClay + self.sedFluxSilt + self.sedFluxSand + self.sedFluxGravel
        self.reporting.reporting(self, pcr, 'SedFlux', sedFlux)

        #-report sediment concentration in the stations (g/l)
        sedConcentration = pcr.cover(sedFlux / (Q * 3600 * 24 * 1e-3), 0)
        self.reporting.reporting(self, pcr, 'SedConc', sedConcentration)