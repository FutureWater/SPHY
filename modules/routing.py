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

print('routing module imported')

def ROUT(pcr, q, oldq, flowdir, kx):
    rr = (q * 0.001 * pcr.cellarea()) / (24*3600)
    ra = pcr.accuflux(flowdir, rr)
    ra = (1 - kx) * ra + kx * oldq
    return ra

#-init routing processes
def init(self, pcr, config):
    self.FlowDir = pcr.readmap(self.inpath + config.get('ROUTING','flowdir'))
    try:
        self.kx = pcr.readmap(self.inpath + config.get('ROUTING','kx'))
    except:
        self.kx = config.getfloat('ROUTING','kx')
  
#-initial conditions routing
def initial(self, pcr, config):
    #-initial routed total runoff
    try:
        self.QRAold = config.getfloat('ROUT_INIT','QRA_init')
    except:
        try:
            self.QRAold = pcr.readmap(self.inpath + config.get('ROUT_INIT','QRA_init'))
        except:
            self.QRAold = 0
    #-initial routed runoff	for the individual components
    pars = ['RootR', 'RootD', 'Rain', 'Snow', 'Glac', 'Base']
    for i in pars:
        try:
            setattr(self, i + 'RAold', pcr.readmap(self.inpath + config.get('ROUT_INIT', i + 'RA_init')))
            setattr(self, i + 'RA_FLAG', True)
        except:
            try:
                setattr(self, i + 'RAold', config.getfloat('ROUT_INIT', i + 'RA_init'))
                setattr(self, i + 'RA_FLAG', True)
            except:
                setattr(self, i + 'RA_FLAG', False)
    
#-dynamic routing processes
def dynamic(self, pcr, TotR):
    #-Rout total runoff
    Q = self.routing.ROUT(pcr, TotR, self.QRAold, self.FlowDir, self.kx)
    self.QRAold = Q
    self.reporting.reporting(self, pcr, 'QallRAtot', Q)
    if self.mm_rep_FLAG == 1 and self.QTOT_mm_FLAG == 1:
        self.QTOTSubBasinTSS.sample(((Q * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) * 1000)

    #-Routing of surface runoff, root drainage, rain, snow, glacier and baseflow
    pars = ['RootR', 'RootD', 'Rain', 'Snow', 'Glac', 'Base']
    for i in pars:
        if getattr(self, i + 'RA_FLAG') == 1:
            try:
                ParsRA = self.routing.ROUT(pcr, getattr(self, i + 'R'), getattr(self, i + 'RAold'), self.FlowDir, self.kx)
            except:
                ParsRA = self.routing.ROUT(pcr, eval(i + 'R'), getattr(self, i + 'RAold'), self.FlowDir, self.kx)
            setattr(self, i + 'RAold', ParsRA)
            self.reporting.reporting(self, pcr, i + 'RAtot', ParsRA)
            if self.mm_rep_FLAG == 1 and getattr(self, 'Q' + i.upper() + '_mm_FLAG') == 1:
                setattr(self, 'Q' + i.upper() + 'SubBasinTSS.sample', ((ParsRA * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) *1000)
                    
    return Q