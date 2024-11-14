# The travel time routing algorithm
# Copyright (C) 2020-2023 Joris Eekhout / Spanish National Research Council (CEBAS-CSIC)
# Email: jeekhout@cebas.csic.es
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

    #-Determine area, wetted perimeter and hydraulic radius for the channel
    areaChannel = self.travel_time_routing.area(pcr, waterDepth, waterDepth, self.channelWidth, 0)
    wettedPerimeterChannel = self.travel_time_routing.wettedPerimeter(pcr, waterDepth, waterDepth, self.channelWidth, 0, self.Ncell)
    hydraulicRadiusChannel = self.travel_time_routing.hydraulicRadius(pcr, areaChannel, wettedPerimeterChannel)

    #-Determine area, wetted perimeter and hydraulic radius for the floodplain
    areaFP = self.travel_time_routing.area(pcr, waterDepth - self.channelDepth, waterDepth - self.channelDepth, (self.floodplainWidth - self.channelWidth) / 2, 0)
    wettedPerimeterFP = self.travel_time_routing.wettedPerimeter(pcr, waterDepth - self.channelDepth, waterDepth - self.channelDepth, (self.floodplainWidth - self.channelWidth) / 2, 0, 1)
    hydraulicRadiusFP = self.travel_time_routing.hydraulicRadius(pcr, areaFP, wettedPerimeterFP)

    #-Determine compound Manning's coefficient for compound channels (Lotter method)
    manningCompound = (wettedPerimeter * hydraulicRadius**(5./3.)) / (((wettedPerimeterFP * hydraulicRadiusFP**(5./3.)) / self.manningFP) * 2 + (wettedPerimeterChannel * hydraulicRadiusChannel**(5./3.)) / self.manningChannel)

    #-Assign compound Manning's n to cells where water depth exceeds channel depth
    manning = pcr.max(self.manningChannel, pcr.ifthenelse(waterDepth > self.channelDepth, manningCompound, self.manningChannel))

    #-Determine characteristic distance
    flowVelocity = (hydraulicRadius)**(2./3.) * ((self.slopeChannel)**(0.5))/ manning * self.dT # meter/day

    return flowVelocity, hydraulicRadius

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

#-Function to rout water using the accuflux travel time algorithm
def flow_velocity_iteration(self, pcr, qOld):
    #-Set maximum absolute and relative differences to a high dummy value 
    diffAbsMax = 1e6

    #-Initiate counter
    counter = 0

    #######################################################################################################################
    #-While loop to obtain flow velocity
    while diffAbsMax > self.thresholdAbs and counter < self.maxIterations:
        #-Increase counter
        counter += 1

        #-Determine flow velocity (m/day)
        flowVelocity, hydraulicRadius = self.travel_time_routing.flowVelocity(self, pcr, self.waterDepth)

        #-Determine discharge (m3/s)
        Q = pcr.accutraveltimefractionflux(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity), self.QFRAC) / self.dT

        #-Determine absolute difference (%) between current and previous discharge
        diffAbs = pcr.abs(Q - qOld)

        #-Determine flow velocity in m/s
        u = pcr.max(flowVelocity / self.dT, 1e-5)

        #-Store discharge to compare with the discharge of the next iteration
        qOld = Q

        #-Determine map maximum for absolute and relative maximum
        diffAbsMax = pcr.cellvalue(pcr.mapmaximum(diffAbs), 1)[0]

        #-Determine water depth (m)
        self.waterDepth = self.travel_time_routing.waterDepth(pcr, Q, u, self.channelDepth, self.channelWidth, self.floodplainWidth)

        #-Determine hillslope roughness and assign to hillslopes (both channel and floodplain)
        if self.ErosionFLAG:
            self.roughness.dynamic(self, pcr)
        else:
            self.manningHillslope = self.manningRill
        
        #-Update channel and floodplain manning
        self.manningChannel = pcr.ifthenelse(self.channelHillslope == 2, self.manningRill, self.manningChannel)
        self.manningFP = pcr.ifthenelse(self.channelHillslope == 2, self.manningHillslope, self.manningFP)

        #-Update channel manning when conservation module is used
        if self.SedTransFLAG:
            if self.conservationFLAG == 1:
                self.manningChannel = pcr.ifthenelse(self.conservationMeasures != 0, self.n_TC_conservation, self.manningChannel)

        #-Determine flow velocity (m/day)
        flowVelocity, hydraulicRadius = self.travel_time_routing.flowVelocity(self, pcr, self.waterDepth)

    return Q, u, hydraulicRadius

#-init routing processes
def init(self, pcr, pcrm, config, np):
    #-Read map input
    self.input.input(self, config, pcr, 'manningChannel', 'ROUTING', 'channelManning', 0)
    self.input.input(self, config, pcr, 'manningRill', 'ROUTING', 'rillManning', 0)
    self.input.input(self, config, pcr, 'manningFP', 'ROUTING', 'floodplainManning', 0)
    self.input.input(self, config, pcr, 'channelWidth', 'ROUTING', 'channelWidth', 0)
    self.input.input(self, config, pcr, 'channelDepth', 'ROUTING', 'channelDepth', 0)
    self.input.input(self, config, pcr, 'floodplainWidth', 'ROUTING', 'floodplainWidth', pcr.celllength())
    self.input.input(self, config, pcr, 'channelHillslope', 'ROUTING', 'channelHillslope', 1)
    self.input.input(self, config, pcr, 'inundationThreshold', 'ROUTING', 'inundationThreshold', 0)

    #-in case of a differentiation between channels and hillslopes
    if pcr.pcr2numpy(pcr.mapmaximum(self.channelHillslope), 1)[0, 0] == 2:
        #-Set floodplain width to cell width for hillslope cells
        self.floodplainWidth = pcr.ifthenelse(self.channelHillslope == 2, pcr.celllength(), self.floodplainWidth)

        #-Set channel depth to small value for hillslope cells
        self.channelDepth = pcr.ifthenelse(self.channelHillslope == 2, 1e-3, self.channelDepth)

        #-Read min and max rill width
        self.input.input(self, config, pcr, 'minRillWidth', 'ROUTING', 'minRillWidth', 0.05)
        self.input.input(self, config, pcr, 'maxRillWidth', 'ROUTING', 'maxRillWidth', 0.3)

        #-Determine rill width based on min and max and set rill depth equal to rill width
        self.rillWidth = self.travel_time_routing.rillDimensions(pcr, self)
        self.rillDepth = self.rillWidth

        #-Set floodplain width to cell width for hillslope cells
        self.floodplainWidth = pcr.ifthenelse(self.channelHillslope == 2, pcr.celllength(), self.floodplainWidth)

        #-Set channel depth to rill depth for hillslope cells
        self.channelDepth = pcr.ifthenelse(self.channelHillslope == 2, self.rillDepth, self.channelDepth)

        #-Set channel width to rill width for hillslope cells
        self.channelWidth = pcr.ifthenelse(self.channelHillslope == 2, self.rillWidth, self.channelWidth)

    #-Read stability thresholds
    self.thresholdAbs = config.getfloat('ROUTING', 'thresholdIteration')
    self.maxIterations = config.getfloat('ROUTING', 'maxIterations')

    #-Calculate timestep used in the model (sec)
    self.dT = (3600 * 24)

    #-Use windowaverage to smooth elevation of channel cells
    self.channelSmooth = config.getint('ROUTING', 'channelSmooth')
   
    #-determine stream order of all channels
    self.streamOrder = pcr.ifthen(self.channelHillslope == 1, pcr.streamorder(self.FlowDir))
    
    #-clump stream order map
    self.clumpedStreamOrder = pcr.clump(self.streamOrder)

    #-determine accuflux map
    self.channel = pcr.ifthen(self.channelHillslope == 1, self.ones)
    self.accuflux = pcr.accuflux(self.FlowDir, 1)

    #-determine number of channel segments
    self.noChannelSegments = int(pcr.pcr2numpy(pcr.mapmaximum(pcr.scalar(self.clumpedStreamOrder)), -9999)[0, 0])

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

            #-determine path from upstream cell to outlet
            self.pathChannel = pcr.path(self.FlowDir, pcr.boolean(self.upstreamChannel))
            
            #-determine accumulated cells to outlet
            self.accufluxPathOutlet = pcr.accuflux(pcr.ifthen(pcr.path(self.FlowDir, pcr.cover(pcr.boolean(self.upstreamChannel), 0)) == 1, self.FlowDir), 1)
            
            #-only get accuflux values for the path
            self.accufluxPath = pcr.ifthen(self.pathChannel == 1, self.accufluxPathOutlet)

            #-determine number of cells in path
            self.cellsPathChannel = pcr.pcr2numpy(pcr.mapmaximum(self.accufluxPath), -9999)[0, 0]

            #-determine if the length of the path is smaller than channel smooth criterium
            if (self.cellsPathChannel < self.channelSmooth):
                self.subPathChannel = self.pathChannel
            else:
                self.subPathChannel = pcr.rounddown(self.accufluxPath / ((self.cellsPathChannel + 1) / np.around(self.cellsPathChannel / self.channelSmooth)))

            #-determine average slope in segment
            self.pathChannelSlope = pcr.areaaverage(self.Slope, pcr.nominal(self.subPathChannel))

            #-add smoothed channel segment to map with all smoothed channels
            self.slopeChannel = self.slopeChannel + pcr.cover(self.pathChannelSlope, 0) 

            #-subtract path from channel segment map
            self.channelSegment = pcr.ifthen(pcr.boolean(pcr.scalar(self.channelSegment) - pcr.scalar(self.pathChannel)) == 1, self.ones)

    #-update channel slope
    self.slopeChannel = pcr.ifthenelse(self.channelHillslope == 1, self.slopeChannel, self.Slope)

    #-Calculate slope to downstream neighbour (m/m) and channel length through cell (m)
    distanceToDownstreamCell = pcrm.distancetodownstreamcell(self.FlowDir)
    self.dX = distanceToDownstreamCell / pcr.cos(pcr.atan(self.slopeChannel))

    #-Set number of channels per cell to 1
    self.Ncell = 1

    if self.ResFLAG == 0 and self.LakeFLAG == 0:
        self.QFRAC = 1

#-initial conditions routing
def initial(self, pcr, config):
    #-Read initial routed total runoff (m3/s)
    self.input.input(self, config, pcr, 'QRAold', 'ROUT_INIT', 'QRA_init', 0)

    #-Read initial channel water depth (m)
    self.input.input(self, config, pcr, 'waterDepth', 'ROUT_INIT', 'H_init', 0)

    #-Determine initial channel storage (m3)
    self.channelStorage = self.waterDepth * self.channelWidth * self.dX

#-dynamic routing processes
def dynamic(self, pcr, TotR):
    #-Update channelStorage by adding the total runoff (m3)
    self.channelStorage += TotR * 0.001 * pcr.cellarea()

    #-Determine flow velocity (m/day)
    flowVelocity, hydraulicRadius = self.travel_time_routing.flowVelocity(self, pcr, self.waterDepth)

    # estimating channel discharge (m3/s)
    Q = pcr.accutraveltimefractionflux(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity), self.QFRAC) / self.dT

    #-Determine flow velocity for first iteration (m/s)
    u = pcr.max(flowVelocity, 1e-0) / self.dT
    
    #-Determine water depth for first iteration (m)
    self.waterDepth = self.travel_time_routing.waterDepth(pcr, Q, u, self.channelDepth, self.channelWidth, self.floodplainWidth)

    #-Determine flow velocity (m/s), water depth (m) and resulting discharge (m3/s) through iteration
    Q, u, hydraulicRadius = self.travel_time_routing.flow_velocity_iteration(self, pcr, Q)
    self.reporting.reporting(self, pcr, 'FlowV', u)
    self.reporting.reporting(self, pcr, 'HydraulicRadius', hydraulicRadius)
    self.reporting.reporting(self, pcr, 'WaterDepth', self.waterDepth)

    #-Update inundation frequency
    self.inundationFrequency = pcr.ifthenelse(self.waterDepth > (self.channelDepth + self.inundationThreshold), pcr.scalar(1), pcr.scalar(0))
    self.reporting.reporting(self, pcr, 'InundationFreq', self.inundationFrequency)

    #-Determine flow velocity in m/day
    flowVelocity = u * self.dT

    #-Update channelStorage after routing (m3)
    self.channelStorage = pcr.accutraveltimefractionstate(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity), self.QFRAC) 

    #-Report discharge
    self.reporting.reporting(self, pcr, 'QallRAtot', Q)
    if self.mm_rep_FLAG == 1 and self.QTOT_mm_FLAG == 1:
        self.QTOTSubBasinTSS.sample(((Q * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) * 1000)
                 
    return Q, u, hydraulicRadius