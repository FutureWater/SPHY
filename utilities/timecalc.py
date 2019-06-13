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


#-Function to return the julian day of the year
def julian(self):
    y= self.curdate.year
    start= self.datetime.datetime(y,1,1).toordinal()
    current= self.curdate.toordinal()
    day= current-start+1
    return day, 1

#-Function to calculate the number of timesteps for the model run
def timesteps(self):
    nrTimeSteps = (self.enddate - self.startdate).days + 1
    print('Running SPHY for '+str(self.startdate.day)+'-'+str(self.startdate.month)+'-'+str(self.startdate.year)+' through '+str(self.enddate.day)+'-'+str(self.enddate.month)+'-'+str(self.enddate.year))
    print('with '+str(nrTimeSteps)+' daily timesteps')
    return nrTimeSteps