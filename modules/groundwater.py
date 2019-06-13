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

print('groundwater module imported')

#-Function to calculate groundwater recharge
def GroundWaterRecharge(pcr, deltagw, gwrecharge, subperc, glacperc):
    gwseep = (1 - pcr.exp(-1 / deltagw)) * (subperc + glacperc)
    gwrecharge = (pcr.exp(-1 / deltagw) * gwrecharge) + gwseep
    return gwrecharge

#-Function to calculate baseflow
def BaseFlow(pcr, gw, baser, gwrecharge, basethresh, alphagw):
    baser = pcr.ifthenelse(gw <= basethresh, 0, (baser * pcr.exp(-alphagw) + gwrecharge * (1 - pcr.exp(-alphagw))))
    return baser

#-Function to calculate the groundwater height, taken from the bottom of the gw layer (zero reference)
def HLevel(pcr, Hgw, alphagw, gwrecharge, yield_gw):
    Hgw = (Hgw * pcr.exp(-alphagw)) + ((gwrecharge * (1 - pcr.exp(-alphagw))) / (800 * yield_gw * alphagw))
    return Hgw

#-init groundwater processes
def init(self, pcr, config):
    pars = ['GwDepth','GwSat','deltaGw','BaseThresh','alphaGw','YieldGw']
    for i in pars:
        try:
            setattr(self, i, pcr.readmap(self.inpath + config.get('GROUNDW_PARS',i)))
        except:
            setattr(self, i, config.getfloat('GROUNDW_PARS',i))

#-initial conditions groundwater
def initial(self, pcr, config):
    #-initial groundwater recharge
    try:
        self.GwRecharge = config.getfloat('GROUNDW_INIT','GwRecharge')
    except:
        self.GwRecharge = pcr.readmap(self.inpath + config.get('GROUNDW_INIT','GwRecharge'))
    #-initial baseflow
    try:
        self.BaseR = config.getfloat('GROUNDW_INIT','BaseR')
    except:
        self.BaseR = pcr.readmap(self.inpath + config.get('GROUNDW_INIT','BaseR'))
    #-initial groundwater storage
    try:
        self.Gw = config.getfloat('GROUNDW_INIT','Gw')
    except:
        self.Gw = pcr.readmap(self.inpath + config.get('GROUNDW_INIT','Gw'))
    #-initial groundwater level
    try:
        self.H_gw = config.getfloat('GROUNDW_INIT','H_gw')
    except:
        self.H_gw = pcr.readmap(self.inpath + config.get('GROUNDW_INIT','H_gw'))
    self.H_gw = pcr.max((self.RootDepthFlat + self.SubDepthFlat + self.GwDepth)/1000 - self.H_gw, 0)

#-dynamic groundwater processes
def dynamic(self, pcr, ActSubPerc, GlacPerc):
    # GwOld = self.Gw
    #-Groundwater recharge
    self.GwRecharge = self.groundwater.GroundWaterRecharge(pcr,	self.deltaGw, self.GwRecharge, ActSubPerc, GlacPerc)
    #-Report groundwater recharge
    self.reporting.reporting(self, pcr, 'TotGwRechargeF', self.GwRecharge)
    #-Update groundwater storage
    self.Gw = self.Gw + self.GwRecharge
    #-Baseflow
    self.BaseR = self.groundwater.BaseFlow(pcr, self.Gw, self.BaseR, self.GwRecharge, self.BaseThresh, self.alphaGw)
    #-Update groundwater storage
    self.Gw = self.Gw - self.BaseR
    #-Report groundwater storage
    self.reporting.reporting(self, pcr, 'StorGroundW', self.Gw * (1-self.openWaterFrac))
    #-Correct for open-water fraction
    self.BaseR = self.BaseR * (1-self.openWaterFrac)
    #-Report Baseflow
    self.reporting.reporting(self, pcr, 'TotBaseRF', self.BaseR)
    #-Calculate groundwater level
    self.H_gw = self.groundwater.HLevel(pcr, self.H_gw, self.alphaGw, self.GwRecharge, self.YieldGw)
    #-Report groundwater level
    self.reporting.reporting(self, pcr, 'GWL', ((self.SubDepthFlat + self.RootDepthFlat + self.GwDepth)/1000 - self.H_gw)*-1)
