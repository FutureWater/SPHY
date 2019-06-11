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
__version__ = "2.2"
__email__ = "info@futurewater.nl"
__date__ ='23 December 2017'
############################################################################################

#-Function to report the output
def REPM(self, pcr, tot, var, fname, outops, TSS=False, MAP=False):
    if outops == 'Day':
        if TSS:
            TSS.sample(var)
        if MAP:
            self.report(var, self.outpath + fname)
        tot = 0
    elif outops == 'Month':
        dim = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        if self.calendar.isleap(self.curdate.year):
            dim[1] = 29
        else:
            dim[1] = 28
        tot = tot + var
        if self.curdate.day == dim[self.curdate.month-1]:
            if TSS:
                TSS.sample(tot)
            if MAP:
                self.report(tot, self.outpath + fname + 'M')
            tot = 0
    elif outops == 'Year':
        if self.calendar.isleap(self.curdate.year):
            ydays = 366
        else:
            ydays = 365
        tot = tot + var
        if self.timecalc.julian(self)[0] == ydays:
            if TSS:
                TSS.sample(tot)
            if MAP:
                self.report(tot, self.outpath + fname + 'Y')
            tot = 0
    else:
        tot = tot + var
        if self.curdate == self.enddate:
            pcr.report(tot, self.outpath + fname + '.map')
            tot = 0
    return tot
 
#-Function to initialise the reporting
def reporting(self, pcr, tot, var):
    for outops in ['Day','Month','Year','Final']:
        try:
            TSS = eval('self.' + tot + '_' + outops + 'TS')
            try:
                MAP = eval('self.' + tot + '_' + outops + '_map')
                setattr(self, tot + '_'+outops, REPM(self, pcr, eval('self.'+tot+'_'+outops), var, eval('self.'+tot+'_fname'), outops, TSS, MAP))
            except:
                setattr(self, tot + '_'+outops, REPM(self, pcr, eval('self.'+tot+'_'+outops), var, eval('self.'+tot+'_fname'), outops, TSS))
        except:
            try:
                MAP = eval('self.' + tot + '_' + outops + '_map')
                setattr(self, tot + '_'+outops, REPM(self, pcr, eval('self.'+tot+'_'+outops), var, eval('self.'+tot+'_fname'), outops, False, MAP))
            except:
                pass    
    
