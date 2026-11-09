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

#-Function to determine the number of rills per cell based on discharge
def rillHydraulics(pcr, self, Q, manningRill):
    #-determine maximum flow area, maximum hydraulic radius and maximum discharge per cell for a single rill
    self.rillMaxArea = self.rillWidth * self.rillDepth
    self.rillMaxR = self.rillMaxArea / (self.rillWidth + 2.0 * self.rillDepth)
    self.rillMaxQ = (1.0 / manningRill) * self.rillMaxArea * (self.rillMaxR ** (2.0/3.0)) * pcr.max(pcr.sqrt(self.slopeChannel), 0.01)

    #-determine number of cells based on maximum discharge per rill and current discharge
    N_rills = pcr.max(1.0, Q / pcr.max(self.rillMaxQ, 1e-10))

    #-determine discharge per rill (can be lower than maximum discharge per rill in case of only 1 rill, which is the absolute minimum of rills)
    self.Q_rill = pcr.max(Q, 1e-10) / N_rills
    
    #-determine flow area, water depth, wetted perimeter, hydraulic radius and flow velocity per rill
    A_rill = (self.Q_rill * manningRill / pcr.max(pcr.sqrt(self.slopeChannel), 0.01)) ** 0.6
    h_rill = A_rill / self.rillWidth
    P_rill = self.rillWidth + (2.0 * h_rill)
    R_rill = A_rill / P_rill
    u_rill = self.Q_rill / pcr.max(A_rill, 1e-10)

    #-determine the effective channel width, baed on the rill width and number of rills
    self.channelWidth = pcr.ifthenelse(self.channelHillslope == 2, self.rillWidth * N_rills, self.channelWidth)

    return u_rill, R_rill


#-Function to rout water using the accuflux travel time algorithm
def flow_velocity_iteration(self, pcr, QLoopOld):
    #-Set maximum absolute and relative differences to a high dummy value 
    diffAbsMax = 1e6

    #-Initiate counter
    counter = 0

    #######################################################################################################################
    #-While loop to determine hydraulic variables through iteration
    while diffAbsMax > self.thresholdAbs and counter < self.maxIterations:
        #-Add counter
        counter += 1

        #-Determine hillslope roughness and assign to hillslopes (both channel and floodplain)
        if self.ErosionFLAG and self.ErosionModel == 2:
            #-Determine Manning's roughness for in-field deposition
            self.n_field, self.n_field_harvest = self.mmf.manningField(self, pcr, self.waterDepth)

            #-replace roughness for vegetated conditions for tilled soil conditions in case of harvested areas
            if self.harvest_FLAG:
                self.manningHillslope = pcr.ifthenelse(self.Harvested == 1, self.n_field_harvest, self.n_field)
            else:
                self.manningHillslope = self.n_field

            #-in case of conservation with vegetation cover
            if self.vegetationCoverFLAG == 1:
                #-set hillslope mannning to vegetation cover value for months between sowing and harvest of vegetation cover
                self.manningHillslope = pcr.ifthenelse(pcr.pcrand(self.VegetationCover > 0, self.Harvested_VC == 0), self.n_field_VC, self.n_field)

        else:
            self.manningHillslope = self.manningRill
        
        #-in case of conservation with vegetation cover
        if self.vegetationCoverFLAG == 1:
            #-set hillslope mannning to vegetation cover value for months between sowing and harvest of vegetation cover
            self.manningRillUpdate = pcr.ifthenelse(pcr.pcrand(self.VegetationCover > 0, self.Harvested_VC == 0), self.n_field_VC, self.manningRill)
        else:
            self.manningRillUpdate = self.manningRill
        
        #-Update channel and floodplain manning
        self.manningChannel = pcr.ifthenelse(self.channelHillslope == 2, self.manningRillUpdate, self.manningChannel)
        self.manningFP = pcr.ifthenelse(self.channelHillslope == 2, self.manningHillslope, self.manningFP)

        #-Determine flow velocity (m/day) and hydraulic radius (m)
        flowVelocity, hydraulicRadius = self.travel_time_routing.flowVelocity(self, pcr, self.waterDepth)

        #-Determine discharge with travel time function
        Q = pcr.accutraveltimefractionflux(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity), self.QFRAC) / self.dT

        #-Check convergence
        diffAbs = pcr.abs(Q - QLoopOld)
        diffAbsMax = pcr.cellvalue(pcr.mapmaximum(diffAbs), 1)[0]

        #-Replace discharge estimate for comparison in next iteration
        QLoopOld = Q

        #-Determine velocity (m/s) and update water depth (m) for the next iteration
        u = pcr.max(flowVelocity / self.dT, 1e-5)
        self.waterDepth = self.travel_time_routing.waterDepth(pcr, Q, u, self.channelDepth, self.channelWidth, self.floodplainWidth)

        #-Determine flow velocity (m/day) and hydraulic radius (m)
        flowVelocity, hydraulicRadius = self.travel_time_routing.flowVelocity(self, pcr, self.waterDepth)

        #-Determine flow velocity (m/s) and hydraulic radius (m) for rills 
        self.u_rill, self.R_rill = self.travel_time_routing.rillHydraulics(pcr, self, Q, self.manningRillUpdate)

        #-update existing hydraulic radius and flow velocity variables with values for rills (hillslopes)
        hydraulicRadius = pcr.ifthenelse(self.channelHillslope == 2, self.R_rill, hydraulicRadius)
        flowVelocity = pcr.ifthenelse(self.channelHillslope == 2, self.u_rill * self.dT, flowVelocity)
        u = pcr.ifthenelse(self.channelHillslope == 2, self.u_rill, u)

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
    self.input.input(self, config, pcr, 'kx_lat', 'ROUTING', 'kx_lat', 1) 
    self.input.input(self, config, pcr, 'kx_base', 'ROUTING', 'kx_base', 1) 
    
    #-in case of a differentiation between channels and hillslopes
    if pcr.pcr2numpy(pcr.mapmaximum(self.channelHillslope), 1)[0, 0] == 2:
        #-Read rill width and depth
        self.input.input(self, config, pcr, 'rillWidth', 'ROUTING', 'rillWidth', 0.3)
        self.input.input(self, config, pcr, 'rillDepth', 'ROUTING', 'rillDepth', 0.15)
        
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

    #-determine channel map
    self.channel = pcr.ifthen(self.channelHillslope == 1, self.ones)

    #### Smooth channel slope ####
    #-determine LDD map for channels only
    channel_LDD = pcr.lddmask(self.FlowDir, self.channelHillslope == 1)
    
    #-retrieve slope for channels
    smoothed_slope = pcr.ifthen(self.channelHillslope == 1, self.Slope)
    
    #-create map with one's in channel
    ones_map = pcr.ifthen(self.channelHillslope == 1, pcr.scalar(1.0))

    #-for-loop to determine the smoothed slope using moving average using 50 iterations
    for i in range(50):
        
        #-determine the downstream slope per cell and account for the outlet cell
        slope_down = pcr.downstream(channel_LDD, smoothed_slope)
        slope_down = pcr.cover(slope_down, smoothed_slope)
        
        #-determine the upstream slope per cell
        sum_up = pcr.upstream(channel_LDD, smoothed_slope)

        #-count upstream channel cells
        count_up = pcr.upstream(channel_LDD, ones_map)
        
        #-if tributaries exist, average them. If headwater, keep own value
        slope_up = pcr.ifthenelse(count_up > 0, sum_up / pcr.max(count_up, 1.0), smoothed_slope)
        
        #-determine smoothed slope using moving average
        smoothed_slope = (slope_up + smoothed_slope + slope_down) / 3.0

    #-apply a minimum value to prevent negative values
    smoothed_slope = pcr.max(smoothed_slope, 0.0001)

    #-store in slopeChannel
    self.slopeChannel = pcr.ifthenelse(self.channelHillslope == 1, smoothed_slope, self.Slope)

    #### Smooth hillslope slope ####
    #-apply smooth hillslope cells
    slopeHillslope = pcr.ifthen(self.channelHillslope == 2, self.slopeChannel)

    #-determine the slope in the two downstream located cells
    slope_down1 = pcr.cover(pcr.downstream(self.FlowDir, slopeHillslope), slopeHillslope)
    slope_down2 = pcr.cover(pcr.downstream(self.FlowDir, slope_down1), slope_down1)

    #-count how many cells flow directly into the current cell
    up_count = pcr.upstream(self.FlowDir, pcr.scalar(1.0))
    
    #-determine the slope in the two upstream located cells
    slope_up1 = pcr.ifthenelse(up_count > 0, 
                                pcr.upstream(self.FlowDir, slopeHillslope) / up_count, 
                                slopeHillslope)
    slope_up2 = pcr.ifthenelse(up_count > 0, 
                                pcr.upstream(self.FlowDir, slope_up1) / up_count, 
                                slopeHillslope)

    #-apply 5-cell moving average and store in slopeChannel
    self.slopeChannel = pcr.ifthenelse(self.channelHillslope == 2, (slopeHillslope + slope_down1 + slope_down2 + slope_up1 + slope_up2) / 5.0, self.slopeChannel)


    #-Calculate slope to downstream neighbour (m/m) and channel length through cell (m)
    distanceToDownstreamCell = pcrm.distancetodownstreamcell(self.FlowDir)
    self.dX = distanceToDownstreamCell / pcr.cos(pcr.atan(self.slopeChannel))

    #-Set number of channels per cell to 1
    self.Ncell = 1

    if self.ResFLAG == 0 and self.LakeFLAG == 0:
        self.QFRAC = 1

    #-Give each channel cell a unique value
    self.ChannelIDs = pcr.uniqueid(pcr.boolean(self.channel))

    #-Determine subcatchment for each channel cell
    self.AquiferSubcatchments = pcr.subcatchment(self.FlowDir, pcr.nominal(pcr.cover(self.ChannelIDs, 0)))

#-initial conditions routing
def initial(self, pcr, config):
    #-Read initial routed total runoff (m3/s)
    self.input.input(self, config, pcr, 'QRAold', 'ROUT_INIT', 'QRA_init', 0)

    #-Read initial channel water depth (m)
    self.input.input(self, config, pcr, 'waterDepth', 'ROUT_INIT', 'H_init', 0)

    #-Determine initial channel storage (m3)
    self.channelStorage = self.waterDepth * self.channelWidth * self.dX

    #-Set initial lateral storage to 0
    self.LateralStorage = self.ones * 0

#-dynamic routing processes
def dynamic(self, pcr):
    #-Add root drainage to lateral flow storage
    self.LateralStorage += self.RootDR

    #-Calculate the amount that is being released to the total runoff
    lateralFlow = self.LateralStorage * (1.0 - self.kx_lat)

    #-Update lateral flow storage
    self.LateralStorage -= lateralFlow


    #-Convert generated baseflow from mm to volume (m3)
    BaseR_m3 = self.BaseR * 0.001 * pcr.cellarea()
    
    #-Determine subcatchment sum of base flow and concentrate them in the channel
    BaseRSubcatchments_m3 = pcr.areatotal(BaseR_m3, self.AquiferSubcatchments) * pcr.cover(self.channel, 0)

    #-Convert the released volume back to mm
    baseFlow_mm = (BaseRSubcatchments_m3 / pcr.cellarea()) * 1000.0


    #-Determine routed runoff
    routedRunoff = self.RootRR + lateralFlow + baseFlow_mm + self.SnowR + self.GlacR

    #-In case MMF is used
    if self.ErosionFLAG and self.ErosionModel == 2:
        #-determine areas that have been harvested
        if self.harvest_FLAG:
            self.Harvested = self.ones * 0
            self.Harvested = pcr.ifthenelse(self.Harvest < self.Sowing, pcr.ifthenelse(pcr.pcrand(self.Harvest < self.curdate.timetuple().tm_yday, self.Sowing > self.curdate.timetuple().tm_yday), 1, self.Harvested), self.Harvested)
            self.Harvested = pcr.ifthenelse(self.Harvest > self.Sowing, pcr.ifthenelse(pcr.pcror(self.curdate.timetuple().tm_yday > self.Harvest, self.curdate.timetuple().tm_yday < self.Sowing), 1, self.Harvested), self.Harvested)
            self.Harvested = pcr.ifthenelse(self.Harvest == 0, 0, self.Harvested)

    #-Update channel storage by adding the total runoff (m3)
    self.channelStorage += routedRunoff * 0.001 * pcr.cellarea()

    #-Determine flow velocity (m/day)
    flowVelocity_init, hydraulicRadius_init = self.travel_time_routing.flowVelocity(self, pcr, self.waterDepth)

    #-Estimate channel discharge (m3/s)
    Q_init = pcr.accutraveltimefractionflux(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity_init), self.QFRAC) / self.dT

    #-Determine flow velocity for first iteration (m/s)
    u_init = pcr.max(flowVelocity_init, 1e-0) / self.dT

    #-Determine flow velocity (m/s) and hydraulic radius (m) for rills 
    u_rill, R_rill = self.travel_time_routing.rillHydraulics(pcr, self, Q_init, self.manningRill)
    hydraulicRadius_init = pcr.ifthenelse(self.channelHillslope == 2, R_rill, hydraulicRadius_init)
    flowVelocity_init = pcr.ifthenelse(self.channelHillslope == 2, u_rill * self.dT, flowVelocity_init)
    u_init = pcr.ifthenelse(self.channelHillslope == 2, u_rill, u_init)

    #-Determine water depth for first iteration (m)
    self.waterDepth = self.travel_time_routing.waterDepth(pcr, Q_init, u_init, self.channelDepth, self.channelWidth, self.floodplainWidth)

    #-Determine flow velocity (m/s), water depth (m) and resulting discharge (m3/s) through iteration
    Q, u, hydraulicRadius = self.travel_time_routing.flow_velocity_iteration(self, pcr, Q_init)

    #-Determine flow velocity in m/day
    flowVelocity = u * self.dT

    #-Update channel storage
    self.channelStorage = pcr.accutraveltimefractionstate(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity), self.QFRAC) 

    #-Determine inundation frequency
    self.inundationFrequency = pcr.ifthenelse(self.waterDepth > (self.channelDepth + self.inundationThreshold), pcr.scalar(1), pcr.scalar(0))

    #-Reporting
    self.reporting.reporting(self, pcr, 'FlowV', u)
    self.reporting.reporting(self, pcr, 'HydraulicRadius', hydraulicRadius)
    self.reporting.reporting(self, pcr, 'WaterDepth', self.waterDepth)
    self.reporting.reporting(self, pcr, 'InundationFreq', self.inundationFrequency)
    self.reporting.reporting(self, pcr, 'QallRAtot', Q)

    if self.mm_rep_FLAG == 1 and self.QTOT_mm_FLAG == 1:
        self.QTOTSubBasinTSS.sample(((Q * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) * 1000)
                 
    return Q, u, hydraulicRadius