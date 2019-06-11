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

print 'Reservoir module imported'

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
    Q = pcr.max(self.ResKr * self.StorRES * (self.StorRES / self.ResSmax)**1.5, self.StorRES - self.ResSmax)
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