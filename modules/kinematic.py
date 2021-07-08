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

def waterDepth(self):
    waterDepth = (self.alpha * (self.QRAold**self.beta)) / self.channelWidth
    return waterDepth

def wettedPerimeter(waterDepth, channelWidth):
    wettedPerimeter = 2 * waterDepth + channelWidth
    return wettedPerimeter

def alpha(self, pcr, wettedPerimeter):
    alpTerm = (self.n / (pcr.sqrt(self.slopeToDownstreamNeighbour)))**self.beta
    alpPow = (pcr.scalar(2) / pcr.scalar(3)) * self.beta
    alpha = alpTerm * (wettedPerimeter**alpPow)
    # alpha = ((self.n * wettedPerimeter**(2/3)) / (1.49 * self.slopeToDownstreamNeighbour**(1/2)))**self.beta
    return alpha
  
def channelStore(self, pcr):
    channelStore = (self.waterDepth * self.channelArea) / pcr.cellarea()
    return channelStore


#-init routing processes
def init(self, pcr, pcrm, config, np):
    #-Read kinematic wave parameters
    self.beta = config.getfloat('ROUTING', 'beta')
    self.nrTimeSlices = config.getint('ROUTING', 'NrTimeSlices')

    #-Read Manning's n
    # self.n = config.getfloat('ROUTING', 'channelManning')
    self.input.input(self, config, pcr, 'n', 'ROUTING', 'channelManning', 0)
    pcr.report(self.n, self.outpath + "manning.map")
    self.channelWidth = config.getfloat('ROUTING', 'channelWidth')

    #-Calculate timestep used in the model (sec)
    self.dT = (3600 * 24)

    #-Calculate channel length through cell (m)
    distanceToDownstreamCell = pcrm.distancetodownstreamcell(self.FlowDir)
    self.slopeToDownstreamNeighbour = pcrm.slopeToDownstreamNeighbourNotFlat(self.DEM, self.FlowDir, 0.0001)
    self.dX = distanceToDownstreamCell / pcr.cos(pcr.atan(self.slopeToDownstreamNeighbour))
  
#-initial conditions routing
def initial(self, pcr, config):
    #-initial routed total runoff
    self.input.input(self, config, pcr, 'QRAold', 'ROUT_INIT', 'QRA_init', 0)

    #-initial channel water depth
    self.input.input(self, config, pcr, 'waterDepth', 'ROUT_INIT', 'H_init', 0)

    self.newstate = 0

    wettedPerimeter = self.kinematic.wettedPerimeter(self.waterDepth, self.channelWidth)
    self.alpha = self.kinematic.alpha(self, pcr, wettedPerimeter)



#-dynamic routing processes
def dynamic(self, pcr, TotR):
    #-Transform inflow (TotR) from mm to m2/s
    q = (TotR / 1000 * pcr.cellarea()) / (3600 * 24) / self.dX
    # pcr.report(TotR, self.outpath + "TotR.map")

    self.waterDepth = self.kinematic.waterDepth(self)
    wettedPerimeter = self.kinematic.wettedPerimeter(self.waterDepth, self.channelWidth)
    self.alpha = self.kinematic.alpha(self, pcr, wettedPerimeter)
    pcr.report(self.alpha, self.outpath + "alpha_" + str(self.counter).zfill(3) + ".map")
    
    #-Rout total discharge
    Q = pcr.kinematic(self.FlowDir, self.QRAold, q, self.alpha, self.beta, self.nrTimeSlices, self.dT, self.dX)
    # Q = pcr.kinwaveflux(self.FlowDir, self.newstate, q, self.alpha, self.beta, self.nrTimeSlices, self.dT, self.dX)
    # self.newstate = pcr.kinwavestate(self.FlowDir, self.newstate, q, self.alpha, self.beta, self.nrTimeSlices, self.dT, self.dX)
    # self.newstate, Q = pcr.kinwavestate, pcr.kinwaveflux(self.FlowDir, self.newstate, q, self.alpha, self.beta, self.nrTimeSlices, self.dT, self.dX)
    pcr.report(Q, self.outpath + "Q_" + str(self.counter).zfill(3) + ".map")
    # exit()

    #-Store discharge for next time step
    self.QRAold = Q

    H = self.kinematic.waterDepth(self)
    A = H * self.channelWidth
    u = Q / A
    pcr.report(H, self.outpath + "H_" + str(self.counter).zfill(3) + ".map")
    pcr.report(A, self.outpath + "A_" + str(self.counter).zfill(3) + ".map")
    pcr.report(u, self.outpath + "u_" + str(self.counter).zfill(3) + ".map")

    # self.waterDepth = self.kinematic.waterDepth(self)
    # wettedPerimeter = self.kinematic.wettedPerimeter(self.waterDepth, self.channelWidth)
    # self.alpha = self.kinematic.alpha(self, pcr, wettedPerimeter)
    # # pcr.report(self.alpha, self.outpath + "alpha_" + str(self.counter).zfill(3) + ".map")

    # #-Rout total discharge
    # Q = pcr.kinematic(self.FlowDir, self.QRAold, q, self.alpha, self.beta, self.nrTimeSlices, self.dT, self.dX)
    # # newstate, Q = pcr.kinwavestate, pcr.kinwaveflux(self.FlowDir, self.QRAold, q, self.alpha, self.beta, self.nrTimeSlices, self.dT, self.dX)
    # pcr.report(Q, self.outpath + "Q_" + str(self.counter).zfill(3) + ".map")
    # # exit()

    # #-Store discharge for next time step
    # self.QRAold = Q
    
    #-Report discharge
    self.reporting.reporting(self, pcr, 'QallRAtot', Q)
    
    if self.mm_rep_FLAG == 1 and self.QTOT_mm_FLAG == 1:
        self.QTOTSubBasinTSS.sample(((Q * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) * 1000)

                  
    return Q