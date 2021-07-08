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

print('dynamic wave routing module imported')

# def ROUT(pcr, q, oldq, flowdir, kx):
#     rr = (q * 0.001 * pcr.cellarea()) / (24*3600)
#     ra = pcr.accuflux(flowdir, rr)
#     ra = (1 - kx) * ra + kx * oldq
#     return ra

#-init routing processes
def init(self, pcr, config, np):
    #-Read dynamic wave parameters
    self.nrTimeSlices = config.getint('ROUTING', 'NrTimeSlices')

    self.input.input(self, config, pcr, 'channelRoughness', 'ROUTING', 'channelManning')
    self.input.input(self, config, pcr, 'channelDepth', 'ROUTING', 'channelDepth')
    self.channelBottomLevel = self.DEM - self.channelDepth
    self.input.input(self, config, pcr, 'channelBottomWidth', 'ROUTING', 'channelWidth')
    self.input.input(self, config, pcr, 'channelForm', 'ROUTING', 'channelForm')
    self.input.input(self, config, pcr, 'floodplainWidth', 'ROUTING', 'floodplainWidth')
    self.input.input(self, config, pcr, 'channelDepthBC', 'ROUTING', 'channelDepthBC')

    #-Calculate timestep used in the model (sec)
    self.timeStepInSeconds = (3600 * 24)

    #-Calculate channel length through cell (m)
    hillslope = pcr.atan(self.Slope)
    temp = pcr.abs(pcr.abs(pcr.abs(pcr.scalar(self.FlowDir) - 5) - 3) - 1)
    euclideanDistance = pcr.celllength() * pcr.ifthenelse(temp == 1, pcr.scalar(self.clone) * np.sqrt(2), pcr.scalar(self.clone))
    self.channelLength = euclideanDistance / pcr.cos(hillslope)

    # self.structures = pcr.boolean(0)
    # self.structureA = pcr.scalar(0)
    # self.structureB = pcr.scalar(0)
    # self.structureCrestLevel = pcr.scalar(0)
    self.structures = pcr.readmap(self.inpath + config.get('ROUTING', 'structures'))
    self.structureA = pcr.readmap(self.inpath + config.get('ROUTING', 'structureA'))
    self.structureB = pcr.readmap(self.inpath + config.get('ROUTING', 'structureB'))
    self.structureCrestLevel = pcr.readmap(self.inpath + config.get('ROUTING', 'structureCrest'))
  
#-initial conditions routing
def initial(self, pcr, config):
    #-initial routed total runoff
    self.input.input(self, config, pcr, 'QRAold', 'ROUT_INIT', 'QRA_init', 0)

    #-initial channel water depth
    self.input.input(self, config, pcr, 'Hold', 'ROUT_INIT', 'H_init', 0)

#-dynamic routing processes
def dynamic(self, pcr, TotR):
    #-Transform inflow (TotR) from mm to m3/s
    q = (TotR / 1000 * pcr.cellarea()) / (3600 * 24) 

    #-change water depth at boundary to boundary condition map value
    self.Hold = pcr.ifthenelse(self.channelDepthBC > 0, self.channelDepthBC, self.Hold)
    pcr.report(self.Hold, self.outpath + "Hold.map")

    #-Rout total discharge
    # newstate, Q = pcr.kinwavestate, pcr.kinwaveflux(self.FlowDir, self.QRAold, q, self.alpha, self.beta, self.nrTimeSlices, self.dT, self.dX)
    Q = pcr.dynamicwaveq(self.FlowDir, q, self.Hold, self.channelBottomLevel, self.channelRoughness, \
        self.channelLength, self.channelBottomWidth, self.channelDepth, self.channelForm, self.floodplainWidth, \
        self.timeStepInSeconds, self.nrTimeSlices, \
        self.structures, self.structureA, self.structureB, self.structureCrestLevel)
    H = pcr.dynamicwaveh(self.FlowDir, q, self.Hold, self.channelBottomLevel, self.channelRoughness, \
        self.channelLength, self.channelBottomWidth, self.channelDepth, self.channelForm, self.floodplainWidth, \
        self.timeStepInSeconds, self.nrTimeSlices, \
        self.structures, self.structureA, self.structureB, self.structureCrestLevel)
    pcr.report(Q, self.outpath + "Q_" + str(self.counter).zfill(3) + ".map")
    pcr.report(H, self.outpath + "H_" + str(self.counter).zfill(3) + ".map")
    # exit()

    #-determine flow area
    A = H * self.channelBottomWidth + H * self.channelForm
    pcr.report(A, self.outpath + "A_" + str(self.counter).zfill(3) + ".map")

    #-determine flow velocity
    u = Q / A
    pcr.report(u, self.outpath + "u_" + str(self.counter).zfill(3) + ".map")

    #-Store discharge for next time step
    self.QRAold = Q

    #-
    self.Hold = H
    
    #-Report discharge
    self.reporting.reporting(self, pcr, 'QallRAtot', Q)
    
    if self.mm_rep_FLAG == 1 and self.QTOT_mm_FLAG == 1:
        self.QTOTSubBasinTSS.sample(((Q * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) * 1000)

                  
    return Q