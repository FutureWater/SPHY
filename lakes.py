# The Spatial Processes in HYdrology (SPHY) model:
# A spatially distributed hydrological model that calculates soil-water and
# cryosphere processes on a cell-by-cell basis.
#
# Copyright (C) 2013  FutureWater
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
#
# Email: info@futurewater.nl

#-Authorship information-###################################################################
__authors__ = "W. Terink, A. Lutz, G. Simons, W. Immerzeel and P. Droogers"
__copyright__ = "FutureWater"
__license__ = "GPL"
__version__ = "2.1"
__email__ = "info@futurewater.nl"
__date__ ='1 January 2017'
############################################################################################

print 'Lake module imported'

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