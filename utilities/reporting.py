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


#-Function to report the output
def REPM(self, pcr, tot, var, fname, outops, TSS=False, MAP=False, AVG=False):
    dim = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if self.calendar.isleap(self.curdate.year):
        dim[1] = 29

    if self.calendar.isleap(self.curdate.year):
        ydays = 366
    else:
        ydays = 365

    if outops == 'Day':
        if TSS:
            TSS.sample(var)
        if MAP:
            self.report(var, self.outpath + fname)
        tot = 0

    elif outops == 'Month':
        tot = tot + var
        if self.curdate.day == dim[self.curdate.month-1]:
            if TSS:
                TSS.sample(tot)
            if MAP:
                self.report(tot, self.outpath + fname + 'M')
            if AVG:
                self.report(tot / dim[self.curdate.month-1], self.outpath + fname + 'M')
            tot = 0

    elif outops == 'Year':
        if self.timecalc.julian(self)[0] != ydays:
            tot = tot + var
        else:
            tot = tot + var
            if TSS:
                TSS.sample(tot)
            if MAP:
                self.report(tot, self.outpath + fname + 'Y')
            if AVG:
                self.report(tot / ydays, self.outpath + fname + 'Y')
            tot = 0

    elif outops == 'MonthSum':
        if self.curdate.day != dim[self.curdate.month-1] and self.curdate.year >= self.startYear + self.spinUpYears:
            tot[self.curdate.month-1] = tot[self.curdate.month-1] + var
        if self.curdate.day == dim[self.curdate.month-1] and self.curdate.year == self.endYear:
            pcr.report(tot[self.curdate.month-1] / (self.simYears), self.outpath + fname + 'SumM' + str(self.curdate.month).zfill(2) + '.map')

    elif outops == 'YearSum':
        if self.curdate.year >= self.startYear + self.spinUpYears:
            tot = tot + var
        if self.curdate == self.enddate:
            pcr.report(tot / (self.simYears), self.outpath + fname + 'SumY.map')

    elif outops == 'MonthAvg':
        if self.curdate.day != dim[self.curdate.month-1] and self.curdate.year >= self.startYear + self.spinUpYears:
            tot[self.curdate.month-1] = tot[self.curdate.month-1] + var / dim[self.curdate.month-1]
        elif self.curdate.day == dim[self.curdate.month-1] and self.curdate.year == self.endYear:
            pcr.report(tot[self.curdate.month-1] / (self.simYears), self.outpath + fname + 'AvgM' + str(self.curdate.month).zfill(2) + '.map')

    elif outops == 'YearAvg':
        if self.curdate.year >= self.startYear + self.spinUpYears:
            tot = tot + var / ydays
        if self.curdate == self.enddate:
            pcr.report(tot / (self.simYears), self.outpath + fname + 'AvgY.map')

    else:
        if self.curdate != self.enddate:
            tot = tot + var
        else:
            pcr.report(tot, self.outpath + fname + '.map')
            tot = 0
    return tot

#-Function to initialise the reporting
def reporting(self, pcr, tot, var):
    for outops in ['Day','Month','Year','Final','MonthSum','YearSum','MonthAvg','YearAvg']:
        try:
            TSS = eval('self.' + tot + '_' + outops + 'TS')
        except:
            TSS = False
        try:
            MAP = eval('self.' + tot + '_' + outops + '_map')
        except:
            MAP = False
        try:
            AVG = eval('self.' + tot + '_' + outops + '_avg')
        except:
            AVG = False
        if TSS or MAP or AVG:
            setattr(self, tot + '_'+outops, REPM(self, pcr, eval('self.'+tot+'_'+outops), var, eval('self.'+tot+'_fname'), outops, TSS, MAP, AVG))

#-read reporting from csv file
def initial(self, pcr, csv, pcrm):
    #-set reporting options and read initial values
    with open(self.inpath + self.RepTab, 'r') as f:
        next(f) # skip headings
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            if row[0][:1] != '#':
                i = row[0]
                mapoutops = row[1]
                avgoutops = row[2]
                TSoutops = row[3]
                if mapoutops == 'NONE' and avgoutops == 'NONE' and TSoutops == 'NONE':
                    print(i + ' will NOT be reported')
                else:
                    print(i + ' will be reported')
                    fname = row[4]
                    setattr(self, i+'_fname', fname)
                    setattr(self, i, 0.)  # use this instead of the commented part above, because it is more logical to always zero as initial condition for reporting
                    if mapoutops != 'NONE':
                        mapoutops = mapoutops.split("+")
                        for j in mapoutops:
                            if j == 'D':
                                setattr(self, i+'_Day', eval('self.'+i))
                                setattr(self, i+'_Day_map', 1)
                            elif j == 'M':
                                setattr(self, i+'_Month', eval('self.'+i))
                                setattr(self, i+'_Month_map', 1)
                            elif j == 'Y':
                                setattr(self, i+'_Year', eval('self.'+i))
                                setattr(self, i+'_Year_map', 1)
                            elif j == 'MS':
                                setattr(self, i+'_MonthSum', {m: self.DEM * 0 for m in range(0, 13, 1)})
                                setattr(self, i+'_MonthSum_map', 1)
                            elif j == 'YS':
                                setattr(self, i+'_YearSum', eval('self.'+i))
                                setattr(self, i+'_YearSum_map', 1)
                            else:
                                setattr(self, i+'_Final', eval('self.'+i))
                                setattr(self, i+'_Final_map', 1)
                    if avgoutops != 'NONE':
                        avgoutops = avgoutops.split("+")
                        for j in avgoutops:
                            if j == 'M':
                                setattr(self, i+'_Month', eval('self.'+i))
                                setattr(self, i+'_Month_avg', 1)
                            elif j == 'Y':
                                setattr(self, i+'_Year', eval('self.'+i))
                                setattr(self, i+'_Year_avg', 1)
                            elif j == 'MA':
                                setattr(self, i+'_MonthAvg', {m: self.DEM * 0 for m in range(0, 13, 1)})
                                setattr(self, i+'_MonthAvg_avg', 1)
                            elif j == 'YA':
                                setattr(self, i+'_YearAvg', eval('self.'+i))
                                setattr(self, i+'_YearAvg_avg', 1)
                    if TSoutops != 'NONE':
                        TSoutops = TSoutops.split("+")
                        for j in TSoutops:
                            if j == 'D':
                                setattr(self, i+'_Day', eval('self.'+i))
                                setattr(self, i+'_DayTS', eval('pcrm.TimeoutputTimeseries("'+fname+'DTS'+'", self, self.Locations, noHeader=False)'))
                            elif j == 'M':
                                setattr(self, i+'_Month', eval('self.'+i))
                                setattr(self, i+'_MonthTS', eval('pcrm.TimeoutputTimeseries("'+fname+'MTS'+'", self, self.Locations, noHeader=False)'))
                            elif j == 'Y':
                                setattr(self, i+'_Year', eval('self.'+i))
                                setattr(self, i+'_YearTS', eval('pcrm.TimeoutputTimeseries("'+fname+'YTS'+'", self, self.Locations, noHeader=False)'))