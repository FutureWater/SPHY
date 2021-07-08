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
    return area

#-Function to determine the wetted perimeter (m)
def wettedPerimeter(pcr, waterDepth, channelDepth, channelWidth, floodplainWidth):
    wettedPerimeter = pcr.ifthenelse(waterDepth <= channelDepth, 2 * waterDepth + channelWidth, 2 * waterDepth + floodplainWidth)
    return wettedPerimeter

#-Function to determine the water depth (m)
def waterDepth(pcr, discharge, velocity, channelDepth, channelWidth, floodplainWidth):
    area = discharge / velocity
    waterDepth = pcr.ifthenelse(area <= channelDepth * channelWidth, area / channelWidth, (area - channelDepth * channelWidth) / floodplainWidth + channelDepth)
    return waterDepth

#-Function to determine flow velocity from water depth using the Manning's equation
def flowVelocity(self, pcr, waterDepth):
    #-Determine area, wetted perimeter and hydraulic radius for compound cross-section
    area = self.travel_time_routing.area(pcr, waterDepth, self.channelDepth, self.channelWidth, self.floodplainWidth)
    wettedPerimeter = self.travel_time_routing.wettedPerimeter(pcr, waterDepth, self.channelDepth, self.channelWidth, self.floodplainWidth)
    hydraulicRadius = area / wettedPerimeter

    #-Determine area, wetted perimeter and hydraulic radius for one floodplain side
    areaChannel = self.travel_time_routing.area(pcr, waterDepth, waterDepth, self.channelWidth, 0)
    wettedPerimeterChannel = self.travel_time_routing.wettedPerimeter(pcr, waterDepth, waterDepth, self.channelWidth, 0)
    hydraulicRadiusChannel = areaChannel / wettedPerimeterChannel

    #-Determine area, wetted perimeter and hydraulic radius for one floodplain side
    areaFP = self.travel_time_routing.area(pcr, waterDepth - self.channelDepth, waterDepth - self.channelDepth, (self.floodplainWidth - self.channelWidth) / 2, 0)
    wettedPerimeterFP = self.travel_time_routing.wettedPerimeter(pcr, waterDepth - self.channelDepth, waterDepth - self.channelDepth, (self.floodplainWidth - self.channelWidth) / 2, 0)
    hydraulicRadiusFP = areaFP / wettedPerimeterFP

    #-Determine compound Manning's coefficient for compound channels (Lotter method)
    manningCompound = (wettedPerimeter * hydraulicRadius**(5./3.)) / (((wettedPerimeterFP * hydraulicRadiusFP**(5./3.)) / self.manningFP) * 2 + (wettedPerimeterChannel * hydraulicRadiusChannel**(5./3.)) / self.manningChannel)

    #-Assign compound Manning's n to cells where water depth exceeds channel depth
    manning = pcr.max(self.manningChannel, pcr.ifthenelse(waterDepth > self.channelDepth, manningCompound, self.manningChannel))

    #-Determine characteristic distance
    flowVelocity = (hydraulicRadius)**(2./3.) * ((self.slopeToDownstreamNeighbour)**(0.5))/ manning * self.dT # meter/day

    return flowVelocity

#-Determine the number of rills per meter
def numberOfRills(pcr, Flow, MC, S, RR, Re):
    N = 0.66 + 0.69 * pcr.ln(Flow) + 0.91 * pcr.ln(MC) + 2.04 * pcr.ln(S) - 0.37 * pcr.ln(RR) - 0.37 * pcr.ln(Re)
    return N

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

        if self.RillFLAG == 1:
            #-Determine number of rills per meter
            N = pcr.max(self.travel_time_routing.numberOfRills(pcr, qOld * 1000 * 60 / pcr.celllength(), self.RootWater / self.RootDepthFlat, \
                self.slopeToDownstreamNeighbour, self.randomRoughness, self.groundCover), 0)
            
            #-Determine number of rills per cell
            Ncell = pcr.cover(N * pcr.celllength(), 0)

            #-Determine channel width based on number of rills and rill dimensions
            self.channelWidth = pcr.max(pcr.ifthenelse(self.channelRill == 2, Ncell * self.rillWidth, self.channelWidth), 1e-3)

        # pcr.report(N, self.outpath + "N_loop_" + str(counter).zfill(3) + ".map")
        # pcr.report(qOld * 1000 / 60 / pcr.celllength(), self.outpath + "Flow_loop_" + str(counter).zfill(3) + ".map")
        # pcr.report(self.RootWater / self.RootDepthFlat, self.outpath + "MC_loop_" + str(counter).zfill(3) + ".map")

        #-Determine flow velocity (m/day)
        flowVelocity = self.travel_time_routing.flowVelocity(self, pcr, self.waterDepth)

        #-Determine discharge (m3/s)
        Q = pcr.accutraveltimefractionflux(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity), self.QFRAC) / self.dT
        # pcr.report(Q, self.outpath + "Q_loop_" + str(counter).zfill(3) + ".map")

        #-Determine absolute difference (%) between current and previous discharge
        diffAbs = pcr.abs(Q - qOld)
        # pcr.report(diffAbs, self.outpath + "diffAbs_loop_" + str(counter).zfill(3) + ".map")
        diffPerc = diffAbs / qOld * 100
        # pcr.report(diffPerc, self.outpath + "diffPerc_loop_" + str(counter).zfill(3) + ".map")

        #-Determine flow velocity to m/s
        u = pcr.max(flowVelocity, 1e-0) / self.dT
        # pcr.report(u, self.outpath + "u_loop_" + str(counter).zfill(3) + ".map")

        #-Store discharge to compare with the discharge of the next iteration
        qOld = Q

        #-Determine map maximum for absolute and relative maximum
        diffAbsMax = pcr.cellvalue(pcr.mapmaximum(diffAbs), 1)[0]
        diffPercMax = pcr.cellvalue(pcr.mapmaximum(diffPerc), 1)[0]

        #-Determine water depth (m)
        self.waterDepth = self.travel_time_routing.waterDepth(pcr, Q, u, self.channelDepth, self.channelWidth, self.floodplainWidth)
        # pcr.report(self.waterDepth, self.outpath + "h_loop_" + str(counter).zfill(3) + ".map")

    self.counterMax = max([self.counterMax, counter])
    if (self.counterMax == counter):
        self.counterOfMax = self.counter
    print(counter)
    print(str(self.counterMax) + " time step " + str(self.counterOfMax))
    # pcr.report(Ncell, self.outpath + "Ncell_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(qOld * 1000 / 60 / pcr.celllength(), self.outpath + "Flow_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(self.RootWater / self.RootDepthFlat, self.outpath + "MC_" + str(self.counter).zfill(3) + ".map")

    return Q, u

#-init routing processes
def init(self, pcr, pcrm, config, np):
    #-Read map input
    self.input.input(self, config, pcr, 'manningChannel', 'ROUTING', 'channelManning', 0)
    self.input.input(self, config, pcr, 'manningFP', 'ROUTING', 'floodplainManning', 0)
    self.input.input(self, config, pcr, 'channelWidth', 'ROUTING', 'channelWidth', 0)
    self.input.input(self, config, pcr, 'channelDepth', 'ROUTING', 'channelDepth', 0)
    self.input.input(self, config, pcr, 'floodplainWidth', 'ROUTING', 'floodplainWidth', pcr.celllength())

    self.RillFLAG = config.getint('ROUTING', 'RillFLAG')
    if self.RillFLAG == 1:
        self.input.input(self, config, pcr, 'channelRill', 'ROUTING', 'channelRill', 1)
        self.input.input(self, config, pcr, 'groundCover', 'ROUTING', 'groundCover', 0)
        self.input.input(self, config, pcr, 'randomRoughness', 'ROUTING', 'randomRoughness', 0)
        self.rillWidth = config.getfloat('ROUTING', 'rillWidth')
        self.rillDepth = config.getfloat('ROUTING', 'rillDepth')

        #-Set floodplain width to cell width for rill cells
        self.floodplainWidth = pcr.ifthenelse(self.channelRill == 2, pcr.celllength(), self.floodplainWidth)

        #-Set channel depth to rill depth for rill cells
        self.channelDepth = pcr.ifthenelse(self.channelRill == 2, self.rillDepth, self.channelDepth)

    #-Read stability thresholds
    self.thresholdPerc = config.getfloat('ROUTING', 'thresholdPercentage')
    self.thresholdAbs = config.getfloat('ROUTING', 'thresholdAbsolute')
    self.maxIterations = config.getfloat('ROUTING', 'maxIterations')

    #-Calculate timestep used in the model (sec)
    self.dT = (3600 * 24)

    #-Calculate slope to downstream neighbour (m/m) and channel length through cell (m)
    distanceToDownstreamCell = pcrm.distancetodownstreamcell(self.FlowDir)
    self.slopeToDownstreamNeighbour = pcrm.slopeToDownstreamNeighbourNotFlat(self.DEM, self.FlowDir, 0.0001)
    self.dX = distanceToDownstreamCell / pcr.cos(pcr.atan(self.slopeToDownstreamNeighbour))

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

    self.counterMax = 0

#-dynamic routing processes
def dynamic(self, pcr, TotR):
    #-Update channelStorage by adding the total runoff (m3)
    self.channelStorage += TotR * 0.001 * pcr.cellarea()
    # pcr.report(self.channelStorage, self.outpath + "channelStorage_" + str(self.counter).zfill(3) + ".map")

    #-Determine flow velocity (m/day)
    flowVelocity = self.travel_time_routing.flowVelocity(self, pcr, self.waterDepth)

    # estimating channel discharge (m3/s)
    Q = pcr.accutraveltimefractionflux(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity), self.QFRAC) / self.dT

    #-Determine flow velocity for first iteration (m/s)
    u = pcr.max(flowVelocity, 1e-0) / self.dT
    
    #-Determine water depth for first iteration (m)
    self.waterDepth = self.travel_time_routing.waterDepth(pcr, Q, u, self.channelDepth, self.channelWidth, self.floodplainWidth)

    #-Determine flow velocity (m/s), water depth (m) and resulting discharge (m3/s) through iteration
    Q, u = self.travel_time_routing.flow_velocity_iteration(self, pcr, Q)

    pcr.report(Q, self.outpath + "Q_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(u, self.outpath + "u_" + str(self.counter).zfill(3) + ".map")
    # pcr.report(self.waterDepth, self.outpath + "h_" + str(self.counter).zfill(3) + ".map")

    #-Determine flow velocity in m/day
    flowVelocity = u * self.dT

    #-Update channelStorage after routing (m3)
    self.channelStorage = pcr.accutraveltimefractionstate(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity), self.QFRAC) 
    # pcr.report(self.channelStorage, self.outpath + "channelStorage_" + str(self.counter).zfill(3) + ".map")

    #-Report discharge
    self.reporting.reporting(self, pcr, 'QallRAtot', Q)
    if self.mm_rep_FLAG == 1 and self.QTOT_mm_FLAG == 1:
        self.QTOTSubBasinTSS.sample(((Q * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) * 1000)
                 
    return Q, u