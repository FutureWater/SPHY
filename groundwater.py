# The Spatial Processes in HYdrology (SPHY) model:
# A spatially distributed hydrological model that calculates soil-water and
# cryosphere processes on a cell-by-cell basis.
#
# Copyright (C) 2013  Wilco Terink
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
# Email: w.terink@futurewater.nl OR terinkw@gmail.com

#-Authorship information-###################################################################
__author__ = "Wilco Terink"
__copyright__ = "Wilco Terink"
__license__ = "GPL"
__version__ = "2.0"
__email__ = "w.terink@futurewater.nl, terinkw@gmail.com"
__date__ ='1 January 2017'
############################################################################################

print 'groundwater module imported'

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


