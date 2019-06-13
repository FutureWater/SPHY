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

print('Lake module imported')

#-Function that updates the lake storage and lake level given a measured lake level. If no lake
# level is measured, then the actual storage is not updated with a measured level. The function
# returns the updated storage and lake level
def UpdateLakeHStore(self, pcr, pcrm):
    #-buffer actual storage
    OldStorage = self.StorRES
    #-Check if measured lake levels area available
    try:
        LakeLevel = pcr.readmap(pcrm.generateNameT(self.LLevel, self.counter))
        Level = True
    except:
        Level = False
    if Level:
        #-update lake storage according to measured water level
        self.StorRES = pcr.ifthenelse(self.UpdateLakeLevel, pcr.ifthenelse(pcr.defined(LakeLevel), pcr.ifthenelse(self.LakeSH_Func==1,\
            self.LakeSH_exp_a * pcr.exp(self.LakeSH_exp_b * LakeLevel), pcr.ifthenelse(self.LakeSH_Func==2, self.LakeSH_pol_a1 \
            * LakeLevel + self.LakeSH_pol_b, pcr.ifthenelse(self.LakeSH_Func==3, (self.LakeSH_pol_a2 * LakeLevel**2) + \
            self.LakeSH_pol_a1 * LakeLevel + self.LakeSH_pol_b, (self.LakeSH_pol_a3 * LakeLevel**3) + (self.LakeSH_pol_a2 \
            * LakeLevel**2) + (self.LakeSH_pol_a1 * LakeLevel + self.LakeSH_pol_b)))), self.StorRES), self.StorRES)
        # prevent storage becoming negative for whatever reason
        self.StorRES = pcr.max(self.StorRES, 0)
        #-Update the lake level based on the storage for lakes where no levels are measured
        LakeLevel = pcr.ifthenelse(self.UpdateLakeLevel, pcr.ifthenelse(pcr.defined(LakeLevel), LakeLevel, \
            pcr.ifthenelse(self.LakeHS_Func==1, self.LakeHS_exp_a * pcr.exp(self.LakeHS_exp_b * self.StorRES), pcr.ifthenelse(self.LakeHS_Func==2, self.LakeHS_pol_a1 * \
            self.StorRES + self.LakeHS_pol_b, pcr.ifthenelse(self.LakeHS_Func==3, (self.LakeHS_pol_a2 * self.StorRES**2) + \
            self.LakeHS_pol_a1 * self.StorRES + self.LakeHS_pol_b, (self.LakeHS_pol_a3 * self.StorRES**3) + (self.LakeHS_pol_a2 *\
            self.StorRES**2) + self.LakeHS_pol_a1 * self.StorRES + self.LakeHS_pol_b)))), pcr.ifthenelse(self.LakeHS_Func==1, \
            self.LakeHS_exp_a * pcr.exp(self.LakeHS_exp_b * self.StorRES), pcr.ifthenelse(self.LakeHS_Func==2, self.LakeHS_pol_a1 * \
            self.StorRES + self.LakeHS_pol_b, pcr.ifthenelse(self.LakeHS_Func==3, (self.LakeHS_pol_a2 * self.StorRES**2) + \
            self.LakeHS_pol_a1 * self.StorRES + self.LakeHS_pol_b, (self.LakeHS_pol_a3 * self.StorRES**3) + (self.LakeHS_pol_a2 *\
            self.StorRES**2) + self.LakeHS_pol_a1 * self.StorRES + self.LakeHS_pol_b))))

    else:
        # if no lake level map is available, then calculate the h based on storages
        LakeLevel = pcr.ifthenelse(self.LakeHS_Func==1, self.LakeHS_exp_a * pcr.exp(self.LakeHS_exp_b * self.StorRES), \
            pcr.ifthenelse(self.LakeHS_Func==2, self.LakeHS_pol_a1 * self.StorRES + self.LakeHS_pol_b, pcr.ifthenelse(\
            self.LakeHS_Func==3, (self.LakeHS_pol_a2 * self.StorRES**2) + self.LakeHS_pol_a1 * self.StorRES + self.LakeHS_pol_b,\
            (self.LakeHS_pol_a3 * self.StorRES**3) + (self.LakeHS_pol_a2 * self.StorRES**2) + self.LakeHS_pol_a1 * self.StorRES +\
            self.LakeHS_pol_b)))
    self.StorRES = pcr.ifthenelse(self.LakeID != 0, self.StorRES, OldStorage)
    return LakeLevel, self.StorRES

#-function that calculates the fraction of lake storage that is available for routing, and the lake outflow
def QLake(self, pcr, LakeLevel):
    Qout = pcr.ifthenelse(self.LakeQH_Func==1, self.LakeQH_exp_a * pcr.exp(self.LakeQH_exp_b * LakeLevel), pcr.ifthenelse(\
        self.LakeQH_Func==2, self.LakeQH_pol_a1 * LakeLevel + self.LakeQH_pol_b, pcr.ifthenelse(self.LakeQH_Func==3, \
        (self.LakeQH_pol_a2 * LakeLevel**2) + self.LakeQH_pol_a1 * LakeLevel + self.LakeQH_pol_b, (self.LakeQH_pol_a3 * \
        LakeLevel**3) + (self.LakeQH_pol_a2 * LakeLevel**2) + self.LakeQH_pol_a1 * LakeLevel + self.LakeQH_pol_b)))
    Qout = pcr.max(0, Qout)
    Qout = Qout * 3600 * 24  #-convert to m3/d
    Qout = pcr.cover(Qout, 0) #-for non-lake cells, Qout is zero
    return Qout

#-init processes lake
def init(self, pcr, config):
    #-set the option to calculate the lake inflow, outflow and storage per component
    pars = ['RootR','RootD','Rain','Snow','Glac','Base']
    for i in pars:
        var = 'Rep' + i + '_FLAG'
        setattr(self, var, config.getint('REPORTING', var))

    pcr.setglobaloption('matrixtable')
    # nominal map with lake IDs
    self.LakeID = pcr.cover(pcr.readmap(self.inpath + config.get('LAKE','LakeId')), 0)
    # lookup table with function for each lake (exp, 1-order poly, 2-order poly, 3-order poly)
    LakeFunc_Tab = self.inpath + config.get('LAKE', 'LakeFunc')
    # lookup table with Qh-coeficients for each lake
    LakeQH_Tab = self.inpath + config.get('LAKE', 'LakeQH')
    # lookup table with Sh-coeficients for each lake
    LakeSH_Tab = self.inpath + config.get('LAKE', 'LakeSH')
    # lookup table with hS-coeficients for each lake
    LakeHS_Tab = self.inpath + config.get('LAKE', 'LakeHS')
    # create lake coefficient maps
    self.LakeQH_Func = pcr.lookupnominal(LakeFunc_Tab, 1, self.LakeID)
    self.LakeSH_Func = pcr.lookupnominal(LakeFunc_Tab, 2, self.LakeID)
    self.LakeHS_Func = pcr.lookupnominal(LakeFunc_Tab, 3, self.LakeID)
    # Read QH coefficients
    self.LakeQH_exp_a = pcr.lookupscalar(LakeQH_Tab, 1, self.LakeID)
    self.LakeQH_exp_b = pcr.lookupscalar(LakeQH_Tab, 2, self.LakeID)
    self.LakeQH_pol_b = pcr.lookupscalar(LakeQH_Tab, 3, self.LakeID)
    self.LakeQH_pol_a1 = pcr.lookupscalar(LakeQH_Tab, 4, self.LakeID)
    self.LakeQH_pol_a2 = pcr.lookupscalar(LakeQH_Tab, 5, self.LakeID)
    self.LakeQH_pol_a3 = pcr.lookupscalar(LakeQH_Tab, 6, self.LakeID)
    # Read SH coefficients
    self.LakeSH_exp_a = pcr.lookupscalar(LakeSH_Tab, 1, self.LakeID)
    self.LakeSH_exp_b = pcr.lookupscalar(LakeSH_Tab, 2, self.LakeID)
    self.LakeSH_pol_b = pcr.lookupscalar(LakeSH_Tab, 3, self.LakeID)
    self.LakeSH_pol_a1 = pcr.lookupscalar(LakeSH_Tab, 4, self.LakeID)
    self.LakeSH_pol_a2 = pcr.lookupscalar(LakeSH_Tab, 5, self.LakeID)
    self.LakeSH_pol_a3 = pcr.lookupscalar(LakeSH_Tab, 6, self.LakeID)
    # Read HS coefficients
    self.LakeHS_exp_a = pcr.lookupscalar(LakeHS_Tab, 1, self.LakeID)
    self.LakeHS_exp_b = pcr.lookupscalar(LakeHS_Tab, 2, self.LakeID)
    self.LakeHS_pol_b = pcr.lookupscalar(LakeHS_Tab, 3, self.LakeID)
    self.LakeHS_pol_a1 = pcr.lookupscalar(LakeHS_Tab, 4, self.LakeID)
    self.LakeHS_pol_a2 = pcr.lookupscalar(LakeHS_Tab, 5, self.LakeID)
    self.LakeHS_pol_a3 = pcr.lookupscalar(LakeHS_Tab, 6, self.LakeID)
    #-read water level maps and parameters if available
    try:
        self.UpdateLakeLevel = pcr.readmap(self.inpath + config.get('LAKE','updatelakelevel'))
        self.LLevel = self.inpath + config.get('LAKE','LakeFile')
        print('measured lake levels will be used to update lake storage')
    except:
        pass
    pcr.setglobaloption('columntable')

#-initial conditions lakes
def initial(self, pcr, config):
    LakeStor_Tab = self.inpath + config.get('LAKE', 'LakeStor')
    self.StorRES = pcr.cover(pcr.lookupscalar(LakeStor_Tab, 1, self.LakeID), 0) * 10**6  # convert to m3
    #-Qfrac for lake cells should be zero, else 1
    self.QFRAC = pcr.ifthenelse(self.LakeID != 0, pcr.scalar(0), 1)

#-initial conditions reporting lakes
def initial_reporting(self, pcr, pcrm):
    self.LakeInTSS = pcrm.TimeoutputTimeseries("LakeInTSS", self, self.LakeID, noHeader=True)
    self.LakeOutTSS = pcrm.TimeoutputTimeseries("LakeOutTSS", self, self.LakeID, noHeader=True)
    self.LakeStorTSS = pcrm.TimeoutputTimeseries("LakeStorTSS", self, self.LakeID, noHeader=True)
    #-set reporting of water balances for individual components
    pars = ['RootR','RootD','Rain','Snow','Glac','Base']
    for i in pars:
        if eval('self.' + i + 'RA_FLAG') and getattr(self, 'Rep' + i + '_FLAG'):
            setattr(self, 'Lake' + i + 'InTSS', pcrm.TimeoutputTimeseries('Lake' + i + 'InTSS', self, self.ResID, noHeader=True))
            setattr(self, 'Lake' + i + 'OutTSS', pcrm.TimeoutputTimeseries('Lake' + i + 'OutTSS', self, self.ResID, noHeader=True))
            setattr(self, 'Lake' + i + 'StorTSS', pcrm.TimeoutputTimeseries('Lake' + i + 'StorTSS', self, self.ResID, noHeader=True))
