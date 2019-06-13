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

print('Reservoir module imported')

#-Advanced reservoir     
def QAdv(self, pcr):
    DayNo = self.timecalc.julian(self)[0]
    #-determine if it is flood or dry season
    S1 = pcr.ifthenelse(self.ResFlStart < self.ResFlEnd, pcr.ifthenelse(DayNo>=self.ResFlStart, pcr.ifthenelse(DayNo<=self.ResFlEnd, pcr.boolean(1), pcr.boolean(0)), pcr.boolean(0)),\
                        pcr.ifthenelse(DayNo>=self.ResFlEnd, pcr.ifthenelse(DayNo>=self.ResFlStart, pcr.boolean(1), pcr.boolean(0)), pcr.ifthenelse(DayNo<=self.ResFlEnd, \
                        pcr.ifthenelse(DayNo<=self.ResFlStart, pcr.boolean(1), pcr.boolean(0)), pcr.boolean(0))))

    S_avail = pcr.max(self.StorRES - self.ResPVOL, 0)
    Q = pcr.max(pcr.ifthenelse(S1, self.ResMaxFl * S_avail / (self.ResEVOL - self.ResPVOL),\
        self.ResDemFl * S_avail / (self.ResEVOL - self.ResPVOL)), self.StorRES - self.ResEVOL)
    return Q

#-Simple reservoir    
def QSimple(self, pcr):
    Q = pcr.max(pcr.min(self.ResKr * self.StorRES * (self.StorRES / self.ResSmax)**self.ResB, self.StorRES), self.StorRES - self.ResSmax)
    return Q
    
#-Calculates reservoir outflow and the fraction to release, depending on the type of reservoir (simple or advanced)    
def QRes(self, pcr):
    if self.ResSimple and self.ResAdvanced:
        Qout = pcr.ifthenelse(self.ResFunc==1, QSimple(self, pcr), pcr.ifthenelse(self.ResFunc==2,\
        QAdv(self, pcr), 0))
    elif self.ResSimple:
        Qout = pcr.ifthenelse(self.ResFunc==1, QSimple(self, pcr), 0)
    else:
        Qout = pcr.ifthenelse(self.ResFunc==2, QAdv(self, pcr), 0)
    
    return Qout

#-init processes reservoirs
def init(self, pcr, config):
    #-set the option to calculate the reservoir inflow, outflow and storage per component
    pars = ['RootR','RootD','Rain','Snow','Glac','Base']
    for i in pars:
        var = 'Rep' + i + '_FLAG'
        setattr(self, var, config.getint('REPORTING', var))

    pcr.setglobaloption('matrixtable')
    # nominal map with reservoir IDs
    self.ResID = pcr.cover(pcr.readmap(self.inpath + config.get('RESERVOIR','ResId')), 0)
    # boolean map with stations that are not reservoirs
    self.LocationsNoRes = pcr.ifthenelse(pcr.pcrand(pcr.scalar(self.Locations) > 0, pcr.scalar(self.ResID) == 0), pcr.boolean(1), pcr.boolean(0))

    # lookup table with operational scheme to use (simple or advanced)
    ResFunc_Tab = self.inpath + config.get('RESERVOIR', 'ResFuncStor')
    # Reservoir function
    self.ResFunc = pcr.cover(pcr.lookupscalar(ResFunc_Tab, 1, self.ResID), 0)
    try:
        # lookup table with coefficients for simple reservoirs
        ResSimple_Tab = self.inpath + config.get('RESERVOIR', 'ResSimple')
        # Read coefficients for simple reservoirs
        self.ResKr = pcr.lookupscalar(ResSimple_Tab, 1, self.ResID)
        self.ResB = pcr.lookupscalar(ResSimple_Tab, 2, self.ResID)
        self.ResSmax = pcr.lookupscalar(ResSimple_Tab, 3, self.ResID) * 10**6 # convert to m3
        self.ResSimple = True
    except:
        self.ResSimple = False
    try:
        # lookup table with coefficients for advanced reservoirs
        ResAdvanced_Tab = self.inpath + config.get('RESERVOIR', 'ResAdv')
        # Read coefficients for advanced reservoirs
        self.ResEVOL = pcr.lookupscalar(ResAdvanced_Tab, 1, self.ResID) * 10**6 # convert to m3
        self.ResPVOL = pcr.lookupscalar(ResAdvanced_Tab, 2, self.ResID) * 10**6 # convert to m3
        self.ResMaxFl = pcr.lookupscalar(ResAdvanced_Tab, 3, self.ResID) * 10**6 # convert to m3/d
        self.ResDemFl = pcr.lookupscalar(ResAdvanced_Tab, 4, self.ResID) * 10**6 # convert to m3/d
        self.ResFlStart = pcr.lookupscalar(ResAdvanced_Tab, 5, self.ResID)
        self.ResFlEnd = pcr.lookupscalar(ResAdvanced_Tab, 6, self.ResID)
        self.ResAdvanced = True
    except:
        self.ResAdvanced = False
    pcr.setglobaloption('columntable')

#-initial conditions reservoirs
def initial(self, pcr, config):
    ResStor_Tab = self.inpath + config.get('RESERVOIR', 'ResFuncStor')
    ResStor = pcr.cover(pcr.lookupscalar(ResStor_Tab, 2, self.ResID), 0) * 10**6  # convert to m3
    try:
        self.StorRES = self.StorRES + ResStor
        #-Qfrac for reservoir cells should be zero, else 1
        self.QFRAC = pcr.ifthenelse(self.ResID != 0, pcr.scalar(0), self.QFRAC)
    except:
        self.StorRES = ResStor
        #-Qfrac for reservoir cells should be zero, else 1
        self.QFRAC = pcr.ifthenelse(self.ResID != 0, pcr.scalar(0), 1)

#-initial conditions reporting reservoirs
def initial_reporting(self, pcr, pcrm):
    self.ResInTSS = pcrm.TimeoutputTimeseries("ResInTSS", self, self.ResID, noHeader=True)
    self.ResOutTSS = pcrm.TimeoutputTimeseries("ResOutTSS", self, self.ResID, noHeader=True)
    self.ResStorTSS = pcrm.TimeoutputTimeseries("ResStorTSS", self, self.ResID, noHeader=True)
    self.ResETaTSS = pcrm.TimeoutputTimeseries("ResETaTSS", self, self.ResID, noHeader=True)
    self.ResInCalTSS = pcrm.TimeoutputTimeseries("ResInCalTSS", self, self.ResID, noHeader=True)
    #-set reporting of water balances for individual components
    pars = ['RootR','RootD','Rain','Snow','Glac','Base']
    for i in pars:
        if eval('self.' + i + 'RA_FLAG') and getattr(self, 'Rep' + i + '_FLAG'):
            setattr(self, 'Res' + i + 'InTSS', pcrm.TimeoutputTimeseries('Res' + i + 'InTSS', self, self.ResID, noHeader=True))
            setattr(self, 'Res' + i + 'OutTSS', pcrm.TimeoutputTimeseries('Res' + i + 'OutTSS', self, self.ResID, noHeader=True))
            setattr(self, 'Res' + i + 'StorTSS', pcrm.TimeoutputTimeseries('Res' + i + 'StorTSS', self, self.ResID, noHeader=True))
