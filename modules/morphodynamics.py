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

print('Morphodynamics module imported')


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


# #-Function to rout water using the accuflux travel time algorithm
# def morphodynamics_iteration(self, pcr, np, material, TC):
#     #-Set maximum absolute and relative differences to a high dummy value 
#     diffAbsMax = 1e6
#     diffPercMax = 1e6

#     #-Initiate counter
#     counter = 0

#     #######################################################################################################################
#     #-While loop to obtain flow velocity
#     while diffAbsMax > self.thresholdAbs and diffPercMax > self.thresholdPerc and counter < self.maxIterations:
#         #-Increase counter
#         counter += 1

#         #-determine sediment yield at stations
#         sedYield, sedDep, sedFlux = self.morphodynamics.SedTrans(self, pcr, np, material, TC)

#         #-Determine absolute difference (%) between current and previous discharge
#         diffAbs = pcr.abs(sedFlux - TC)
#         # pcr.report(diffAbs, self.outpath + "diffAbs_loop_" + str(counter).zfill(3) + ".map")
#         diffPerc = diffAbs / TC * 100
#         # pcr.report(diffPerc, self.outpath + "diffPerc_loop_" + str(counter).zfill(3) + ".map")

#         # #-Store discharge to compare with the discharge of the next iteration
#         # fluxOld = sedFlux

#         material = materialOld - TC

#         #-Determine map maximum for absolute and relative maximum
#         diffAbsMax = pcr.cellvalue(pcr.mapmaximum(diffAbs), 1)[0]
#         diffPercMax = pcr.cellvalue(pcr.mapmaximum(diffPerc), 1)[0]


#     self.counterMax = max([self.counterMax, counter])
#     if (self.counterMax == counter):
#         self.counterOfMax = self.counter
#     print(counter)
#     print(str(self.counterMax) + " time step " + str(self.counterOfMax))

#     return material, sedYield, sedDep, sedFlux

#-init morphodynamics processes
def init(self, pcr, pcrm, config):
    #-read input parameters
    self.input.input(self, config, pcr, 'bedThickness', 'MORPHODYNAMICS', 'bedThickness', 0)

    # #-update slope based on initial channel depth
    # self.slopeChannel = pcrm.slopeToDownstreamNeighbourNotFlat(self.DEM - self.channelDepth, self.FlowDir, 0.0001)

    #-set gravel fraction to 0 for hillslope
    self.RootGravelMap = 0

    #-For loop over the sediment classes
    for sedimentClass in self.sedimentClasses:
        #-read channel material percentage
        self.input.input(self, config, pcr, 'channel' + sedimentClass, 'MORPHODYNAMICS', 'channel' + sedimentClass, 0)
        
        #-assign channel material to channels and use soil texture fraction for hillslopes
        setattr(self, "channel" + sedimentClass, pcr.ifthenelse(self.channelHillslope == 1, getattr(self, "channel" + sedimentClass) / 100, getattr(self, "Root" + sedimentClass + "Map")))


    #-determine max channel depth
    self.channelDepthMax = self.channelDepth + self.bedThickness


#-initial morphodynamics processes
def initial(self, pcr):
    #-For loop over the sediment classes
    for sedimentClass in self.sedimentClasses:
        #-assign bed thickness to all sediment classes
        setattr(self, "sedimentStoreInitial" + sedimentClass, getattr(self, "channel" + sedimentClass) * self.bedThickness * self.channelWidth * pcr.celllength() * self.rho_s * 1e-3)
        setattr(self, "sedimentStore" + sedimentClass, getattr(self,  "sedimentStoreInitial" + sedimentClass))
    
    #-determine overall sediment storage
    self.sedimentStore = self.sedimentStoreClay + self.sedimentStoreSilt + self.sedimentStoreSand + self.sedimentStoreGravel


#-dynamic sediment transport processes
def dynamic(self, pcr, pcrm, np, Q, Sed):
    #-For loop over the sediment classes
    for sedimentClass in self.sedimentClasses:

        #-When MMF is used, obtain hillslope erosion per sediment class from attribute
        if sedimentClass == 'Gravel':
            Sed = self.ones * 0
        else:
            if self.ErosionModel == 2:
                Sed = getattr(self, "Sed" + sedimentClass)
            else:
                #-Multiply hillslope erosion by sediment class fraction
                Sed = Sed * getattr(self, "Root" + sedimentClass + "Map")

        #-get transport capacity for sediment class
        TC = getattr(self, "TC" + sedimentClass)

        #-get bed thickness for sediment class
        sedimentStore = getattr(self, "sedimentStore" + sedimentClass)

        #-set sediment storage for hillslopes to minimum defined by initial channel bed depth
        # sedimentStoreInitial = pcr.ifthenelse(self.channelHillslope == 2, getattr(self, "sedimentStoreInitial" + sedimentClass), sedimentStore)
        sedimentStoreInitial = pcr.ifthenelse(self.channelHillslope == 2, 0, sedimentStore) #-no sediment store on hillslopes

        #-determine the amount of sediment available from the channel bed
        # pcr.report(sedimentStore, self.outpath + "sedimentStore_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")

        #-determine total amount available for transport
        material = sedimentStoreInitial + Sed #* 0

        #-route sediment through the catchment
        sedYield, sedDep, sedFlux = self.morphodynamics.SedTrans(self, pcr, np, material, TC)
        # pcr.report(Sed, self.outpath + "Sed_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")
        # pcr.report(TC, self.outpath + "TC_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")
        # pcr.report(sedFlux, self.outpath + "sedFlux_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")
        # pcr.report(sedDep, self.outpath + "sedDep_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")

        #-store sediment depostion as sediment storage
        sedimentStore = sedDep
        # pcr.report(sedimentStore, self.outpath + "sedimentStore_" + sedimentClass + "_" + str(self.counter).zfill(3) + ".map")

        setattr(self, "sedimentStore" + sedimentClass, sedimentStore)

        channelChange = sedimentStore - sedimentStoreInitial

        #-Assign sediment yield, deposition and flux values to sediment class
        setattr(self, "sedYield" + sedimentClass, sedYield)
        setattr(self, "sedDep" + sedimentClass, sedDep)
        setattr(self, "channelChange" + sedimentClass, channelChange)
        setattr(self, "sedFlux" + sedimentClass, sedFlux)

    # #-store old sediment storage
    # sedimentStoreOld = self.sedimentStore

    #-determine total sediment store of all grain size classes
    self.sedimentStore = self.sedimentStoreClay + self.sedimentStoreSilt + self.sedimentStoreSand + self.sedimentStoreGravel
    self.reporting.reporting(self, pcr, 'SdStor', self.sedimentStore)
    # pcr.report(self.sedimentStore, self.outpath + "sedimentStore_" + str(self.counter).zfill(3) + ".map")

    #-determine total channel erosion
    # channelChange = self.sedimentStore - sedimentStoreOld
    channelChange = self.channelChangeClay + self.channelChangeSilt + self.channelChangeSand + self.channelChangeGravel
    self.reporting.reporting(self, pcr, 'ChChng', channelChange)

    #-determine bed thickness based on sediment store
    bedThicknessOld = self.bedThickness
    self.bedThickness = (self.sedimentStore) / (self.channelWidth * pcr.celllength() * self.rho_s * 1e-3)
    # pcr.report(self.bedThickness, self.outpath + "bedThickness_" + str(self.counter).zfill(3) + ".map")

    #-determine channel bed change based on new bed thickness
    channelBedChange = pcr.ifthenelse(self.channelHillslope == 1, self.bedThickness - bedThicknessOld, 0)
    # pcr.report(channelBedChange, self.outpath + "channelBedChange_" + str(self.counter).zfill(3) + ".map")

    #-determine new channel depth
    self.channelDepth = pcr.min(pcr.max(self.channelDepth - channelBedChange, 0.01), self.channelDepthMax)
    self.reporting.reporting(self, pcr, 'ChDpth', self.channelDepth)
    # pcr.report(self.channelDepth, self.outpath + "channelDepth_" + str(self.counter).zfill(3) + ".map")

    # #-update slope based on new channel depth
    # self.slopeChannel = pcrm.slopeToDownstreamNeighbourNotFlat(self.DEMSmooth - self.channelDepth, self.FlowDir, 0.0001)
    # self.reporting.reporting(self, pcr, 'ChSlp', self.slopeChannel)

    #-report the sediment deposition by transport capacity (ton/day)
    sedDep = self.sedDepClay + self.sedDepSilt + self.sedDepSand + self.sedDepGravel
    self.reporting.reporting(self, pcr, 'SedDep', sedDep)
    # pcr.report(sedFlux, self.outpath + "sedDep_" + str(self.counter).zfill(3) + ".map")

    #-report sediment yield in the stations (ton/day)
    sedYield = self.sedYieldClay + self.sedYieldSilt + self.sedYieldSand + self.sedYieldGravel
    self.reporting.reporting(self, pcr, 'SedYld', sedYield)

    #-report sediment flux in the stations (ton/day)
    sedFlux = self.sedFluxClay + self.sedFluxSilt + self.sedFluxSand + self.sedFluxGravel
    self.reporting.reporting(self, pcr, 'SedFlux', sedFlux)
    # pcr.report(sedFlux, self.outpath + "sedFlux_" + str(self.counter).zfill(3) + ".map")

    #-report sediment concentration in the stations (g/l)
    sedConcentration = pcr.cover(sedFlux / (Q * 3600 * 24 * 1e-3), 0)
    self.reporting.reporting(self, pcr, 'SedConc', sedConcentration)
    # pcr.report(sedConcentration, self.outpath + "sedConc_" + str(self.counter).zfill(3) + ".map")
