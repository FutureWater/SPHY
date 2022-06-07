# The Spatial Processes in HYdrology (SPHY) model:
# A spatially distributed hydrological model 
# Copyright (C) 2013-2019  FutureWater
# Email: sphy@futurewater.nl
#
# Authors (alphabetical order):
# P. Droogers, J. Eekhout, W. Immerzeel, S. Khanal, A. Lutz, G. Simons, W. Terink
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

print('travel time routing module imported')

#-Function to determine flow area (m2)
def area(pcr, waterDepth, channelDepth, channelWidth, floodplainWidth):
    area = pcr.ifthenelse(waterDepth <= channelDepth, waterDepth * channelWidth, channelDepth * channelWidth + (waterDepth - channelDepth) * floodplainWidth)
    area = pcr.ifthenelse(waterDepth < 0, 0, area)
    return area

#-Function to determine the wetted perimeter (m)
def wettedPerimeter(pcr, waterDepth, channelDepth, channelWidth, floodplainWidth, numberOfChannels):
    wettedPerimeter = pcr.ifthenelse(waterDepth <= channelDepth, 2 * waterDepth * numberOfChannels + channelWidth, 2 * waterDepth * numberOfChannels + floodplainWidth)
    wettedPerimeter = pcr.ifthenelse(waterDepth < 0, 0, wettedPerimeter)
    return wettedPerimeter

#-Function to determine the hydraulic radius (m)
def hydraulicRadius(pcr, area, wettedPerimeter):
    hydraulicRadius = pcr.cover(area / wettedPerimeter, 0)
    return hydraulicRadius

#-Function to determine the water depth (m)
def waterDepth(pcr, discharge, velocity, channelDepth, channelWidth, floodplainWidth):
    area = discharge / velocity
    waterDepth = pcr.ifthenelse(area <= channelDepth * channelWidth, area / channelWidth, (area - channelDepth * channelWidth) / floodplainWidth + channelDepth)
    return waterDepth

#-Function to determine flow velocity from water depth using the Manning's equation
def flowVelocity(self, pcr, waterDepth):
    #-Determine area, wetted perimeter and hydraulic radius for compound cross-section
    area = self.travel_time_routing.area(pcr, waterDepth, self.channelDepth, self.channelWidth, self.floodplainWidth)
    wettedPerimeter = self.travel_time_routing.wettedPerimeter(pcr, waterDepth, self.channelDepth, self.channelWidth, self.floodplainWidth, self.Ncell)
    hydraulicRadius = self.travel_time_routing.hydraulicRadius(pcr, area, wettedPerimeter)
    # pcr.report(waterDepth, self.outpath + "waterDepth_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(area, self.outpath + "area_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(wettedPerimeter, self.outpath + "wettedPerimeter_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(hydraulicRadius, self.outpath + "hydraulicRadius_" + str(self.counter).zfill(3) + ".map")

    #-Determine area, wetted perimeter and hydraulic radius for the channel
    areaChannel = self.travel_time_routing.area(pcr, waterDepth, waterDepth, self.channelWidth, 0)
    wettedPerimeterChannel = self.travel_time_routing.wettedPerimeter(pcr, waterDepth, waterDepth, self.channelWidth, 0, self.Ncell)
    hydraulicRadiusChannel = self.travel_time_routing.hydraulicRadius(pcr, areaChannel, wettedPerimeterChannel)
    # pcr.report(areaChannel, self.outpath + "areaChannel_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(wettedPerimeterChannel, self.outpath + "wettedPerimeterChannel_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(hydraulicRadiusChannel, self.outpath + "hydraulicRadiusChannel_" + str(self.counter).zfill(3) + ".map")

    #-Determine area, wetted perimeter and hydraulic radius for the floodplain
    areaFP = self.travel_time_routing.area(pcr, waterDepth - self.channelDepth, waterDepth - self.channelDepth, (self.floodplainWidth - self.channelWidth) / 2, 0)
    wettedPerimeterFP = self.travel_time_routing.wettedPerimeter(pcr, waterDepth - self.channelDepth, waterDepth - self.channelDepth, (self.floodplainWidth - self.channelWidth) / 2, 0, 1)
    hydraulicRadiusFP = self.travel_time_routing.hydraulicRadius(pcr, areaFP, wettedPerimeterFP)
    # pcr.report(areaFP, self.outpath + "areaFP_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(wettedPerimeterFP, self.outpath + "wettedPerimeterFP_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(hydraulicRadiusFP, self.outpath + "hydraulicRadiusFP_" + str(self.counter).zfill(3) + ".map")

    #-Determine compound Manning's coefficient for compound channels (Lotter method)
    manningCompound = (wettedPerimeter * hydraulicRadius**(5./3.)) / (((wettedPerimeterFP * hydraulicRadiusFP**(5./3.)) / self.manningFP) * 2 + (wettedPerimeterChannel * hydraulicRadiusChannel**(5./3.)) / self.manningChannel)

    #-Assign compound Manning's n to cells where water depth exceeds channel depth
    manning = pcr.max(self.manningChannel, pcr.ifthenelse(waterDepth > self.channelDepth, manningCompound, self.manningChannel))

    #-Determine characteristic distance
    flowVelocity = (hydraulicRadius)**(2./3.) * ((self.slopeChannel)**(0.5))/ manning * self.dT # meter/day

    return flowVelocity, hydraulicRadius

#-Function to rout water using the accuflux travel time algorithm
def flow_velocity_iteration(self, pcr, qOld):
    #-Set maximum absolute and relative differences to a high dummy value 
    diffAbsMax = 1e6
    diffPercMax = 1e6

    #-Initiate counter
    counter = 0

    #######################################################################################################################
    #-While loop to obtain flow velocity
    while diffAbsMax > self.thresholdAbs and diffPercMax > self.thresholdPerc and counter < self.maxIterations:
        #-Increase counter
        counter += 1

        # if self.ErosionFLAG and self.RillFLAG:
        #     #-Determine number of rills per meter
        #     N = pcr.max(self.erosion.numberOfRills(pcr, qOld * 1000 * 60 / pcr.celllength(), self.RootWater / self.RootDepthFlat, \
        #         self.slopeChannel, self.randomRoughness, self.GC), 0)
            
        #     #-Determine number of rills per cell, round values and set minimum value to 1
        #     self.Ncell = pcr.max(pcr.roundoff(pcr.cover(N * pcr.celllength(), 1)), 1)

        #     #-Set number of rills in channel cells to 1
        #     self.Ncell = pcr.ifthenelse(self.channelHillslope == 1, 1, self.Ncell)

        #     #-Determine channel width based on number of rills and rill dimensions
        #     self.channelWidth = pcr.max(pcr.ifthenelse(self.channelHillslope == 2, self.Ncell * self.rillWidth, self.channelWidth), 1e-3)
        # pcr.report(N, self.outpath + "N_loop_" + str(counter).zfill(3) + ".map")
        # pcr.report(qOld * 1000 / 60 / pcr.celllength(), self.outpath + "Flow_loop_" + str(counter).zfill(3) + ".map")
        # pcr.report(self.RootWater / self.RootDepthFlat, self.outpath + "MC_loop_" + str(counter).zfill(3) + ".map")

        #-Determine flow velocity (m/day)
        flowVelocity, hydraulicRadius = self.travel_time_routing.flowVelocity(self, pcr, self.waterDepth)
        # pcr.report(flowVelocity, self.outpath + "flowVelocity_loop_" + str(counter).zfill(3) + ".map")

        #-Determine discharge (m3/s)
        Q = pcr.accutraveltimefractionflux(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity), self.QFRAC) / self.dT
        # pcr.report(qOld, self.outpath + "qOld_loop_" + str(counter).zfill(3) + ".map")
        # pcr.report(Q, self.outpath + "Q_loop_" + str(counter).zfill(3) + ".map")

        #-Determine absolute difference (%) between current and previous discharge
        diffAbs = pcr.abs(Q - qOld)
        # pcr.report(diffAbs, self.outpath + "diffAbs_loop_" + str(counter).zfill(3) + ".map")
        diffPerc = diffAbs / qOld * 100
        # pcr.report(diffPerc, self.outpath + "diffPerc_loop_" + str(counter).zfill(3) + ".map")

        #-Determine flow velocity in m/s
        u = pcr.max(flowVelocity / self.dT, 1e-5)
        # pcr.report(u, self.outpath + "u_loop_" + str(counter).zfill(3) + ".map")

        #-Store discharge to compare with the discharge of the next iteration
        qOld = Q

        #-Determine map maximum for absolute and relative maximum
        diffAbsMax = pcr.cellvalue(pcr.mapmaximum(diffAbs), 1)[0]
        diffPercMax = pcr.cellvalue(pcr.mapmaximum(diffPerc), 1)[0]

        #-Determine water depth (m)
        self.waterDepth = self.travel_time_routing.waterDepth(pcr, Q, u, self.channelDepth, self.channelWidth, self.floodplainWidth)
        # pcr.report(self.waterDepth, self.outpath + "h_loop_" + str(self.counter).zfill(3) + "_" + str(counter).zfill(3) + ".map")

        #-Determine hillslope roughness and assign to hillslopes (both channel and floodplain)
        if self.ErosionFLAG and self.RillFLAG:
            self.roughness.dynamic(self, pcr)
            # self.manningHillslope = self.manningHillslopeTest #-rill roughness equal to channel roughness
            self.manningChannel = pcr.ifthenelse(self.channelHillslope == 2, self.manningHillslopeTest, self.manningChannel)
            if self.conservationFLAG == 1:
                self.manningChannel = pcr.ifthenelse(self.conservationMeasures != 0, self.n_TC_conservation, self.manningChannel)
            self.manningFP = pcr.ifthenelse(self.channelHillslope == 2, self.manningHillslope, self.manningFP)

        #-Determine flow velocity (m/day)
        flowVelocity, hydraulicRadius = self.travel_time_routing.flowVelocity(self, pcr, self.waterDepth)
        # pcr.report(flowVelocity, self.outpath + "flowVelocity_loop_" + str(counter).zfill(3) + ".map")


    self.counterMax = max([self.counterMax, counter])
    if (self.counterMax == counter):
        self.counterOfMax = self.counter
    self.counterTotal += counter
    print(counter)
    print("max " + str(self.counterMax) + " steps in time step " + str(self.counterOfMax))
    print("average " + str(round(self.counterTotal / self.counter, 1)) + " steps")
    
    # pcr.report(self.Ncell, self.outpath + "Ncell_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(u, self.outpath + "flowVelocity_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(self.waterDepth, self.outpath + "waterDepth_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(self.channelDepth, self.outpath + "channelDepth_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(self.channelWidth, self.outpath + "channelWidth_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(self.manningChannel, self.outpath + "manningChannel_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(self.manningFP, self.outpath + "manningFP_" + str(self.counter).zfill(3) + ".map")

    # pcr.report(qOld * 1000 / 60 / pcr.celllength(), self.outpath + "Flow_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(self.RootWater / self.RootDepthFlat, self.outpath + "MC_" + str(self.counter).zfill(3) + ".map")

    return Q, u, hydraulicRadius

#-init routing processes
def init(self, pcr, pcrm, config, np):
    #-Read map input
    self.input.input(self, config, pcr, 'manningChannel', 'ROUTING', 'channelManning', 0)
    self.input.input(self, config, pcr, 'manningHillslopeTest', 'ROUTING', 'hillslopeManning', 0)
    self.input.input(self, config, pcr, 'manningFP', 'ROUTING', 'floodplainManning', 0)
    self.input.input(self, config, pcr, 'channelWidth', 'ROUTING', 'channelWidth', 0)
    self.input.input(self, config, pcr, 'channelDepth', 'ROUTING', 'channelDepth', 0)
    self.input.input(self, config, pcr, 'floodplainWidth', 'ROUTING', 'floodplainWidth', pcr.celllength())
    self.input.input(self, config, pcr, 'channelHillslope', 'ROUTING', 'channelHillslope', 1)

    if pcr.pcr2numpy(pcr.mapmaximum(self.channelHillslope), 1)[0, 0] == 2:
        #-Set floodplain width to cell width for hillslope cells
        self.floodplainWidth = pcr.ifthenelse(self.channelHillslope == 2, pcr.celllength(), self.floodplainWidth)

        #-Set channel depth to small value for hillslope cells
        self.channelDepth = pcr.ifthenelse(self.channelHillslope == 2, 1e-3, self.channelDepth)

    #-Read stability thresholds
    self.thresholdPerc = config.getfloat('ROUTING', 'thresholdPercentage')
    self.thresholdAbs = config.getfloat('ROUTING', 'thresholdAbsolute')
    self.maxIterations = config.getfloat('ROUTING', 'maxIterations')

    #-Calculate timestep used in the model (sec)
    self.dT = (3600 * 24)

    #-Use windowaverage to smooth elevation of channel cells
    self.channelSmooth = config.getint('ROUTING', 'channelSmooth')
    # self.DEMChannelSmooth = pcr.ifthen(self.channelHillslope == 1, self.DEM)
    # pcr.setglobaloption("unitcell")
    # self.DEMChannelSmooth = pcr.windowaverage(self.DEMChannelSmooth, self.channelSmooth)
    # pcr.setglobaloption("unittrue")
    # self.DEMChannelSmooth = pcr.ifthenelse(self.channelHillslope == 1, self.DEMChannelSmooth, self.DEM)
    # self.slopeChannel = pcrm.slopeToDownstreamNeighbourNotFlat(self.DEMChannelSmooth - self.channelDepth, self.FlowDir, 0.0001)

    #-determine average slope per stream based on stream order and averaged per subcatchment
    # self.Basin = pcr.subcatchment(self.FlowDir, self.ResSedID)
    
    #-determine stream order of all channels
    self.streamOrder = pcr.ifthen(self.channelHillslope == 1, pcr.streamorder(self.FlowDir))
    # pcr.report(pcr.streamorder(self.FlowDir), self.outpath + "streamOrder.map")
    
    #-clump stream order map
    self.clumpedStreamOrder = pcr.clump(self.streamOrder)
    # pcr.report(self.clumpedStreamOrder, self.outpath + "clumpedStreamOrder.map")
    # exit()


    #-determine accuflux map
    self.channel = pcr.ifthen(self.channelHillslope == 1, self.ones)
    self.accuflux = pcr.accuflux(self.FlowDir, 1)
    self.accufluxChannel = pcr.accuflux(pcr.ifthen(self.channelHillslope == 1, self.FlowDir), 1)
    # pcr.report(self.accufluxChannel, self.outpath + "accufluxChannel.map")

    #-determine number of channel segments
    self.noChannelSegments = int(pcr.pcr2numpy(pcr.mapmaximum(pcr.scalar(self.clumpedStreamOrder)), -9999)[0, 0])

    # #-initiate smoothed channel map
    # self.channelDEMSmooth = self.ones * 0

    #-initiate channel slope map
    self.slopeChannel = self.ones * 0

    #-for-loop over each invidual channel segment
    for segment in range(1, self.noChannelSegments + 1):
        #-extract channel segment
        self.channelSegment = pcr.ifthen(self.clumpedStreamOrder == segment, self.ones)

        #-while-loop over channel segments
        while(pcr.pcr2numpy(pcr.maptotal(pcr.scalar(self.channelSegment)), -9999)[0, 0] > 0):
            #-find cell with lowest accuflux value
            self.upstreamChannel = pcr.ifthenelse((self.accuflux * self.channelSegment) == pcr.mapminimum(pcr.scalar(self.channelSegment) * self.accuflux), self.ones, self.ones * 0)
            # pcr.report(self.upstreamChannel, self.outpath + "upstreamChannel_" + str(segment) + ".map")

            #-determine path from upstream cell to outlet
            self.pathChannel = pcr.path(self.FlowDir, pcr.boolean(self.upstreamChannel))
            # pcr.report(self.pathChannel, self.outpath + "pathChannel_" + str(segment) + ".map")

            # #-extract DEM for channel path
            # self.pathChannelDEM = pcr.ifthen(pcr.boolean(self.pathChannel), self.DEM)
            # # pcr.report(self.pathChannelDEM, self.outpath + "pathChannelDEM_" + str(segment) + ".map")
            # # pcr.report(pcr.scalar(self.pathChannel) * self.DEM, self.outpath + "pathChannelDEM2_" + str(segment) + ".map")

            # #-determine number of cells in channel segment
            # self.channelSegmentTotal = pcr.pcr2numpy(pcr.maptotal(pcr.scalar(self.channelSegment)), -9999)[0, 0]
            # # print(self.channelSegmentTotal)
            # # exit()

            # #-if the number of cells of the channel segment is larger than the channel smooth parameter, then interpolate, else use the DEM values
            # if (self.channelSegmentTotal > self.channelSmooth):
            #     #-smooth channel segment
            #     pcr.setglobaloption("unitcell")
            #     self.pathChannelDEMSmooth = pcr.ifthenelse(pcr.boolean(self.pathChannel), pcr.windowaverage(self.pathChannelDEM, self.channelSmooth), 0)
            #     pcr.setglobaloption("unittrue")
            #     # pcr.report(self.pathChannelDEMSmooth, self.outpath + "pathChannelDEMSmooth_" + str(segment) + ".map")
            #     # exit()
            # else:
            #     self.pathChannelDEMSmooth = self.pathChannelDEM

            # # pcr.report(self.pathChannelDEMSmooth, self.outpath + "pathChannelDEMSmooth_" + str(segment) + ".map")
            
            # #-prevent double counting at river junctions
            # self.pathChannelDEMSmooth = pcr.ifthenelse(self.channelDEMSmooth != 0, 0, self.pathChannelDEMSmooth)

            # #-add smoothed channel segment to map with all smoothed channels
            # self.channelDEMSmooth = self.channelDEMSmooth + pcr.cover(self.pathChannelDEMSmooth, 0) 
            # # pcr.report(self.channelDEMSmooth, self.outpath + "channelDEMSmooth_" + str(segment) + ".map")

            # #-extract slope for channel path
            # self.pathChannelSlope = pcr.ifthen(pcr.boolean(self.pathChannel), self.Slope)
            # pcr.report(self.pathChannelSlope, self.outpath + "pathChannelSlope.map")

            # pcr.report(pcr.nominal(self.pathChannel), self.outpath + "pathChannel.map")

            
            self.accufluxPathOutlet = pcr.accuflux(pcr.ifthen(pcr.path(self.FlowDir, pcr.cover(pcr.boolean(self.upstreamChannel), 0)) == 1, self.FlowDir), 1)
            # pcr.report(self.accufluxPathOutlet, self.outpath + "accufluxPathOutlet_" + str(segment) + ".map")
            
            self.accufluxPath = pcr.ifthen(self.pathChannel == 1, self.accufluxPathOutlet)
            # self.accufluxPath = pcr.ifthen(self.pathChannel == 1, self.accufluxChannel)
            # self.accufluxPath = self.accufluxPath - pcr.mapminimum(self.accufluxPath) + 1
            # pcr.report(self.accufluxPath, self.outpath + "accufluxPath_" + str(segment) + ".map")

            #-determine number of cells in path
            self.cellsPathChannel = pcr.pcr2numpy(pcr.mapmaximum(self.accufluxPath), -9999)[0, 0]

            # if (self.cellsPathChannel > self.channelSmooth):
            #     self.cellsPathChannel / self.channelSmooth

            # self.accufluxPath = pcr.accuflux(pcr.ifthen(self.pathChannel == 1, self.FlowDir), 1)
            # pcr.report(self.accufluxChannel, self.outpath + "accufluxPath.map")

            # print(self.cellsPathChannel)

            if (self.cellsPathChannel < self.channelSmooth):
                self.subPathChannel = self.pathChannel
            else:
                self.subPathChannel = pcr.rounddown(self.accufluxPath / ((self.cellsPathChannel + 1) / np.around(self.cellsPathChannel / self.channelSmooth)))
                # pcr.report(self.subPathChannel, self.outpath + "subPathChannel_" + str(segment) + ".map")
                # exit()
                

            #-determine average slope in segment
            self.pathChannelSlope = pcr.areaaverage(self.Slope, pcr.nominal(self.subPathChannel))
            # pcr.report(self.pathChannelSlope, self.outpath + "pathChannelSlope_" + str(segment) + ".map")
            # if (self.cellsPathChannel > (self.channelSmooth * 4)):
            #     exit()

            #-add smoothed channel segment to map with all smoothed channels
            self.slopeChannel = self.slopeChannel + pcr.cover(self.pathChannelSlope, 0) 

            #-subtract path from channel segment map
            self.channelSegment = pcr.ifthen(pcr.boolean(pcr.scalar(self.channelSegment) - pcr.scalar(self.pathChannel)) == 1, self.ones)
            # pcr.report(self.channelSegment, self.outpath + "channelSegment.map")
            # exit()
            
            # print(pcr.pcr2numpy(pcr.maptotal(pcr.scalar(self.channelSegment)), -9999)[0, 0])

    # self.Streams = pcr.scalar(self.Channel) * pcr.scalar(self.Basin) * 10 + pcr.scalar(self.Channel) * pcr.scalar(self.StreamOrder)
    # self.SlopeStreams = pcr.areaaverage(self.Slope, pcr.nominal(self.Streams))
    # self.SlopeStreams = pcr.ifthenelse(self.Channel == 1, self.SlopeStreams, self.Slope)

    # pcr.report(self.channelDEMSmooth, self.outpath + "channelDEMSmooth.map")
    # exit()
    
    # #-add smoothed channel DEM map to DEM
    # self.DEMSmooth = pcr.ifthenelse(self.channelHillslope == 1, self.channelDEMSmooth, self.DEM)
    # # pcr.report(self.DEMSmooth, self.outpath + "DEMSmooth.map")

    # self.slopeChannel = pcrm.slopeToDownstreamNeighbourNotFlat(self.DEMSmooth - self.channelDepth, self.FlowDir, 0.0001)

    # pcr.report(self.slopeChannel, self.outpath + "slopeChannel.map")

    self.slopeChannel = pcr.ifthenelse(self.channelHillslope == 1, self.slopeChannel, self.Slope)
    # pcr.report(self.slopeChannel, self.outpath + "slopeChannel_final.map")
    # exit()

    #-Calculate slope to downstream neighbour (m/m) and channel length through cell (m)
    distanceToDownstreamCell = pcrm.distancetodownstreamcell(self.FlowDir)
    self.dX = distanceToDownstreamCell / pcr.cos(pcr.atan(self.slopeChannel))

    #-Set number of channels per cell to 1
    self.Ncell = 1

    if self.ResFLAG == 0 and self.LakeFLAG == 0:
        self.QFRAC = 1

    #-initiate counter total
    self.counterTotal = 0

  
#-initial conditions routing
def initial(self, pcr, config):
    #-Read initial routed total runoff (m3/s)
    self.input.input(self, config, pcr, 'QRAold', 'ROUT_INIT', 'QRA_init', 0)

    #-Read initial channel water depth (m)
    self.input.input(self, config, pcr, 'waterDepth', 'ROUT_INIT', 'H_init', 0)

    #-Determine initial channel storage (m3)
    self.channelStorage = self.waterDepth * self.channelWidth * self.dX

    self.counterMax = 0

#-dynamic routing processes
def dynamic(self, pcr, TotR):
    # #-Determine hillslope roughness
    # if self.ErosionFLAG and self.RillFLAG:
    #     #-Determine hillslope roughness and assign to hillslopes (both channel and floodplain)
    #     self.roughness.dynamic(self, pcr)
    #     self.manningChannel = pcr.ifthenelse(self.channelHillslope == 2, self.manningHillslope, self.manningChannel)
    #     self.manningFP = pcr.ifthenelse(self.channelHillslope == 2, self.manningHillslope, self.manningFP)

    #-Update channelStorage by adding the total runoff (m3)
    self.channelStorage += TotR * 0.001 * pcr.cellarea()
    # pcr.report(self.channelStorage, self.outpath + "channelStorage_" + str(self.counter).zfill(3) + ".map")

    #-Determine flow velocity (m/day)
    flowVelocity, hydraulicRadius = self.travel_time_routing.flowVelocity(self, pcr, self.waterDepth)

    # estimating channel discharge (m3/s)
    Q = pcr.accutraveltimefractionflux(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity), self.QFRAC) / self.dT

    #-Determine flow velocity for first iteration (m/s)
    u = pcr.max(flowVelocity, 1e-0) / self.dT
    # pcr.report(u, self.outpath + "u1_" + str(self.counter).zfill(3) + ".map")
    
    #-Determine water depth for first iteration (m)
    self.waterDepth = self.travel_time_routing.waterDepth(pcr, Q, u, self.channelDepth, self.channelWidth, self.floodplainWidth)

    #-Determine flow velocity (m/s), water depth (m) and resulting discharge (m3/s) through iteration
    Q, u, hydraulicRadius = self.travel_time_routing.flow_velocity_iteration(self, pcr, Q)
    self.reporting.reporting(self, pcr, 'FlowV', u)

    # pcr.report(Q, self.outpath + "Q_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(u, self.outpath + "flowVelocity_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(self.waterDepth, self.outpath + "h_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(hydraulicRadius, self.outpath + "R_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(self.manningChannel, self.outpath + "manningChannel_" + str(self.counter).zfill(3) + ".map")

    #-Determine flow velocity in m/day
    flowVelocity = u * self.dT

    #-Update channelStorage after routing (m3)
    self.channelStorage = pcr.accutraveltimefractionstate(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity), self.QFRAC) 
    # pcr.report(self.channelStorage, self.outpath + "channelStorage_" + str(self.counter).zfill(3) + ".map")

    #-Report discharge
    self.reporting.reporting(self, pcr, 'QallRAtot', Q)
    if self.mm_rep_FLAG == 1 and self.QTOT_mm_FLAG == 1:
        self.QTOTSubBasinTSS.sample(((Q * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) * 1000)
                 
    return Q, u, hydraulicRadius