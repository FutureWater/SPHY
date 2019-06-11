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

# Advanced routing that is used for reservoirs or lakes
print 'Advanced routing module for lakes and reservoirs imported'

#-Function to rout the specific runoff
def ROUT(self, pcr, rvolume, qold, qout, sres):
    # Calculate the discharge Q (m3/d)
    Q = pcr.accufractionflux(self.FlowDir, rvolume, self.QFRAC)
    # Re-calculate Q, based on qold en kx, and assign Qout for cells being lake/reservoir
    Q = pcr.ifthenelse(self.QFRAC==0, qout, (1-self.kx) * Q + (qold*24*3600) * self.kx)
    # Only calculate inflow for lake/reservoir cells
    Qin = pcr.ifthenelse(self.QFRAC==0, pcr.upstream(self.FlowDir, Q), 0)
    sres = sres - qout + Qin
    Q = Q / (24 * 3600)  #-only convert Q to m3/s
    return sres, Q, Qin