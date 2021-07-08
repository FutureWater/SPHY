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

print('kinematic wave routing module imported')
import pcraster as pcr

# def waterDepth(self):
#     waterDepth = (self.alpha * (self.QRAold**self.beta)) / self.channelWidth
#     return waterDepth

# def wettedPerimeter(waterDepth, channelWidth):
#     wettedPerimeter = 2 * waterDepth + channelWidth
#     return wettedPerimeter

# def alpha(self, pcr, wettedPerimeter):
#     alpTerm = (self.n / (pcr.sqrt(self.slopeToDownstreamNeighbour)))**self.beta
#     alpPow = (pcr.scalar(2) / pcr.scalar(3)) * self.beta
#     alpha = alpTerm * (wettedPerimeter**alpPow)
#     # alpha = ((self.n * wettedPerimeter**(2/3)) / (1.49 * self.slopeToDownstreamNeighbour**(1/2)))**self.beta
#     return alpha
  
def channelStore(self, pcr):
    channelStore = (self.waterDepth * self.channelArea) / pcr.cellarea()
    return channelStore

def getValDivZero(x,y,y_lim=1E-39,z_def= 0.):
    #-returns the result of a division that possibly involves a zero
    # denominator; in which case, a default value is substituted:
    # x/y= z in case y > y_lim,
    # x/y= z_def in case y <= y_lim, where y_lim -> 0.
    # z_def is set to zero if not otherwise specified
    return pcr.ifthenelse(y > y_lim,x/pcr.max(y_lim,y),z_def)


def calculate_alpha_and_initial_discharge_for_kinematic_wave(self, pcr, channelStorage, water_height, innundatedFraction, floodDepth): 

    # calculate alpha (dimensionless), which is the roughness coefficient 
    # - for kinewatic wave (see: http://pcraster.geo.uu.nl/pcraster/4.0.0/doc/manual/op_kinematic.html)
    # - based on wetted area (m2) and wetted perimeter (m), as well as self.beta (dimensionless)
    # - assuming rectangular channel
    # - flood innundated areas with 

    # channel wetted area (m2)
    # - the minimum wetted area is: water height x channel width (Edwin introduce this) 
    # - channel wetted area is mainly based on channelStorage and channelLength (Rens's approach)
    channel_wetted_area = water_height * self.channelWidth
    channel_wetted_area = pcr.max(channel_wetted_area,\
                                    channelStorage / self.channelLength)

    # wetted perimeter
    flood_only_wetted_perimeter = floodDepth * (2.0) + \
                                    pcr.max(0.0, innundatedFraction*self.cellArea/self.channelLength - self.channelWidth)
    channel_only_wetted_perimeter = \
                pcr.min(self.channelDepth, getValDivZero(channelStorage, self.channelLength*self.channelWidth, 0.0)) * 2.0 + \
                self.channelWidth
    # total channel wetted perimeter (unit: m)
    channel_wetted_perimeter = channel_only_wetted_perimeter + \
                                    flood_only_wetted_perimeter
    # minimum channel wetted perimeter = 10 cm
    channel_wetted_perimeter = pcr.max(0.1, channel_wetted_perimeter)                             
        
    # corrected Manning's coefficient: 
    if self.floodPlain:
        usedManningsN = ((channel_only_wetted_perimeter/channel_wetted_perimeter) *      self.manningsN**(1.5) + \
                            (  flood_only_wetted_perimeter/channel_wetted_perimeter) * self.floodplainManN**(1.5))**(2./3.)
    else:
        usedManningsN = self.manningsN
    
    # alpha (dimensionless) and initial estimate of channel discharge (m3/s)
    #
    alpha = (usedManningsN*channel_wetted_perimeter**(2./3.)*self.slopeToDownstreamNeighbour**(-0.5))**self.beta  # dimensionless
    dischargeInitial = pcr.ifthenelse(alpha > 0.0,\
                                        (channel_wetted_area / alpha)**(1.0/self.beta), 0.0)       # unit: m3
    
    return (alpha, dischargeInitial)    

def estimate_length_of_sub_time_step(self, pcr): 

    # estimate the length of sub-time step (unit: s):
    # - the shorter is the better
    # - estimated based on the initial or latest sub-time step discharge (unit: m3/s)
    # 
    # length_of_sub_time_step = pcr.ifthenelse(self.subDischarge > 0.0, self.water_height * self.cellArea / self.subDischarge, self.dT)
    # TODO: Check this logic with Rens!

    u = self.Q / (self.channelDepth * self.channelWidth) 
    pcr.report(u, self.outpath + "u_" + str(self.counter).zfill(3) + ".map")
        
    length_of_sub_time_step = pcr.ifthenelse(self.subDischarge > 0.0, 
                                self.channelLength / u, self.dT)
    # pcr.report(length_of_sub_time_step, self.outpath + "length_of_sub_time_step.map")

    # determine the number of sub time steps (based on Rens van Beek's method)
    #
    critical_condition = (length_of_sub_time_step < self.dT)  & \
                            (self.water_height > self.critical_water_height) & \
                            (self.FlowDir != pcr.ldd(5))
    # pcr.report(critical_condition, self.outpath + "critical_condition.map")
    #
    number_of_sub_time_steps = self.dT /\
                                pcr.cover(
                                pcr.areaminimum(\
                                pcr.ifthen(critical_condition, \
                                            length_of_sub_time_step),self.clone),\
                                            self.dT/self.limit_num_of_sub_time_steps)   
    
    number_of_sub_time_steps = 1.25 * number_of_sub_time_steps + 1
    number_of_sub_time_steps = pcr.roundup(number_of_sub_time_steps)
    # pcr.report(number_of_sub_time_steps, self.outpath + "number_of_sub_time_steps.map")
    #
    number_of_loops = int(max(1.0, pcr.cellvalue(pcr.mapmaximum(number_of_sub_time_steps), 1)[0]))     # minimum number of sub_time_steps = 1 
    number_of_loops = int(min(self.limit_num_of_sub_time_steps, number_of_loops))
    # print(self.limit_num_of_sub_time_steps)
    # exit()
    
    # actual length of sub-time step (s)
    length_of_sub_time_step = self.dT / number_of_loops

    return (length_of_sub_time_step, number_of_loops)                               

#-init routing processes
def init(self, pcr, pcrm, config, np):
    #-Read kinematic wave parameters
    self.beta = config.getfloat('ROUTING', 'beta')
    self.nrTimeSlices = config.getint('ROUTING', 'NrTimeSlices')

    #-Read Manning's n
    # self.n = config.getfloat('ROUTING', 'channelManning')
    self.input.input(self, config, pcr, 'manningsN', 'ROUTING', 'channelManning', 0)
    self.input.input(self, config, pcr, 'channelWidth', 'ROUTING', 'channelWidth', 0)
    self.input.input(self, config, pcr, 'channelDepth', 'ROUTING', 'channelDepth', 0)

    #-Calculate timestep used in the model (sec)
    self.dT = (3600 * 24)

    #-Calculate channel length through cell (m)
    distanceToDownstreamCell = pcrm.distancetodownstreamcell(self.FlowDir)
    self.slopeToDownstreamNeighbour = pcrm.slopeToDownstreamNeighbourNotFlat(self.DEM, self.FlowDir, 0.0001)
    self.channelLength = distanceToDownstreamCell / pcr.cos(pcr.atan(self.slopeToDownstreamNeighbour))
    
    #-to make the script run
    self.innundatedFraction = 0
    self.floodDepth = 0
    self.floodPlain = False

    # critical water height (m) used to select stable length of sub time step in kinematic wave methods/approaches
    self.critical_water_height = 0.05;  # used in Van Beek et al. (2011)

    # courantNumber criteria for numerical stability in kinematic wave methods/approaches
    self.courantNumber = 0.50

    # empirical values for minimum number of sub-time steps:
    design_flood_speed = 1.00 # m/s
    design_length_of_sub_time_step   = pcr.cellvalue(
                                        pcr.mapminimum(
                                        self.courantNumber * self.channelLength / design_flood_speed), 1)[0]
    # self.limit_num_of_sub_time_steps = np.ceil(self.dT / design_length_of_sub_time_step)
    self.limit_num_of_sub_time_steps = 1000


#-initial conditions routing
def initial(self, pcr, config):
    #-initial routed total runoff
    self.input.input(self, config, pcr, 'QRAold', 'ROUT_INIT', 'QRA_init', 0)

    #-initial channel water depth
    self.input.input(self, config, pcr, 'water_height', 'ROUT_INIT', 'H_init', 0)

    self.newstate = 0

    # wettedPerimeter = self.kinematic.wettedPerimeter(self.waterDepth, self.channelWidth)
    # self.alpha = self.kinematic.alpha(self, pcr, wettedPerimeter)
    self.channelStorage = self.water_height * self.channelWidth * self.channelLength
    self.Q = 0.001
    self.subDischarge = self.Q


#-dynamic routing processes
def dynamic(self, pcr, TotR):
    #-Transform inflow (TotR) from mm to m2/s
    # q = (TotR / 1000 * pcr.cellarea()) / (3600 * 24) / self.dX
    self.runoff = TotR * 0.001

    # channelStorage that will be given to the ROUTING operation:
    channelStorageForRouting = pcr.max(0.0, self.channelStorage)                              # unit: m3

    # estimate of water height (m)
    # - needed to estimate the length of sub-time step and 
    #     also to estimate the channel wetted area (for the calculation of alpha and dischargeInitial)
    self.water_height = channelStorageForRouting / (self.channelWidth * self.channelLength)

    # estimate the length of sub-time step (unit: s):
    length_of_sub_time_step, number_of_loops = estimate_length_of_sub_time_step(self, pcr)

    print(number_of_loops)

    acc_discharge_volume = pcr.scalar(0.0)   # unit: m3

    #######################################################################################################################
    for i_loop in range(number_of_loops):
        
        nonIrrReturnFlow = 0

        # update channelStorageForRouting after runoff and return flow from non irrigation demand
        channelStorageForRouting += (self.runoff + nonIrrReturnFlow) * self.channelWidth * self.channelLength * length_of_sub_time_step/self.dT  # unit: m3
        # pcr.report(channelStorageForRouting, self.outpath + "channelStorageForRouting_" + str(i_loop).zfill(3) + ".map")

        # alpha parameter and initial/estimate discharge variable - needed for kinematic wave calculation 
        alpha, dischargeInitial = \
                calculate_alpha_and_initial_discharge_for_kinematic_wave(self, pcr, channelStorageForRouting, \
                                                                                self.water_height, \
                                                                                self.innundatedFraction, self.floodDepth) 
        # pcr.report(dischargeInitial, self.outpath + "dischargeInitial_" + str(i_loop).zfill(3) + ".map")

        # discharge (m3/s) based on the KINEMATIC WAVE approximation
        #~ logger.debug('start pcr.kinematic')
        self.subDischarge = pcr.kinematic(self.FlowDir, dischargeInitial, 0.0, 
                                            alpha, self.beta, \
                                            1, length_of_sub_time_step, self.channelLength)
        self.subDischarge = pcr.max(0.0, pcr.cover(self.subDischarge, 0.0))

        # update channelStorage (m3) after lateral flows in channels
        storage_change_in_volume  = pcr.upstream(self.FlowDir, self.subDischarge * length_of_sub_time_step) - self.subDischarge * length_of_sub_time_step 
        channelStorageForRouting += storage_change_in_volume 

        # total discharge_volume (m3) until this present i_loop
        acc_discharge_volume += self.subDischarge * length_of_sub_time_step


    # channel discharge (m3/s) = self.Q
    Q = acc_discharge_volume / self.dT
    self.Q = Q
    # Q = self.subDischarge

    # updating channelStorage (after routing)
    self.channelStorage = channelStorageForRouting
    # exit()

    pcr.report(Q, self.outpath + "Q_" + str(self.counter).zfill(3) + ".map")

    #-Report discharge
    self.reporting.reporting(self, pcr, 'QallRAtot', Q)
    
    if self.mm_rep_FLAG == 1 and self.QTOT_mm_FLAG == 1:
        self.QTOTSubBasinTSS.sample(((Q * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) * 1000)

                  
    return Q