# Tool to determine soil hydraulic properties from texture and organic matter input maps,
# based on the pedotransfer function from Saxton & Rawls (2006)
# Copyright (C) 2016-2019 Joris Eekhout / Spanish National Research Council (CEBAS-CSIC)
# Email: jeekhout@cebas.csic.es
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


#-function to calculate the wilting point (1500 kPa, %)
def Dry(pcr, self, sand, clay, OM, bulk):
    #-1500 kPa moisture, first solution, %v (Eq. 1)
    dry_first = -0.024 * sand + 0.487 * clay + 0.006 * OM + 0.005 * sand * OM - 0.013 * clay * OM + 0.068 * sand * clay + 0.031
    #-1500 kPa moisture, %v (Eq. 1)
    dry = dry_first + (0.14 * dry_first - 0.02)
    return dry


#-function to calculate the field capacity (33 kPa, %)
def Field(pcr, self, sand, clay, OM, bulk):
    #-33 kPa moisture, first solution, %v (Eq. 2)
    field_first = -0.251 * sand + 0.195 * clay + 0.011 * OM + 0.006 * sand * OM - 0.027 * clay * OM + 0.452 * sand * clay + 0.299
    #-33 kPa moisture, normal density, %v (Eq. 2)
    field = field_first + (1.283 * field_first**2 - 0.374 * field_first - 0.015)
    return field


#-function to calculate the saturated water content (0 kPa, %)
def Sat(pcr, self, sand, clay, OM, bulk):
    #-SAT-33 kPa moisture, first solution, %v (Eq. 3)
    poros_first = 0.278 * sand + 0.034 * clay + 0.022 * OM - 0.018 * sand * OM - 0.027 * clay * OM - 0.584 * sand * clay + 0.078
    #-SAT-33 kPa moisture, normal density, %v (Eq. 3)
    poros = poros_first + (0.636 * poros_first - 0.107)
    #-Saturated moisture (0 kPa), normal density, %v (Eq. 5)
    sat = poros + self.pedotransfer.Field(pcr, self, sand, clay, OM, bulk) - 0.097 * sand + 0.043
    return sat


#-function to calculate the field capacity adjusted for density (33 kPa, %) and saturated water content adjusted for density (0 kPa, %)
def FieldAdj(pcr, self, sand, clay, OM, bulk):
    #-Normal density, g cm-3 (Eq. 6)
    density = (1-self.pedotransfer.Sat(pcr, self, sand, clay, OM, bulk))*2.65
    #-Adjusted density, g cm-3 (Eq. 7)
    density_adj = density * bulk
    #-Saturated moisture (0 kPa), adjusted density, %v (Eq. 8)
    sat_adj_dens = 1 - density_adj/2.65
    #-33 kPa moisture, adjusted density, %v (Eq. 9)
    field_adj_dens = self.pedotransfer.Field(pcr, self, sand, clay, OM, bulk) + 0.2 * (self.pedotransfer.Sat(pcr, self, sand, clay, OM, bulk) - sat_adj_dens)
    #-SAT-33 kPa moisture, adjusted density, %v (Eq. 10)
    poros_adj_dens = sat_adj_dens - field_adj_dens
    return field_adj_dens, poros_adj_dens, sat_adj_dens


#-function to calculate the saturated hydraulic conductivity (mm/day)
def KSat(pcr, self, sand, clay, OM, bulk):
    #-Calculate wilting, field capacity and porosity
    dry = self.pedotransfer.Dry(pcr, self, sand, clay, OM, bulk)
    temp = self.pedotransfer.FieldAdj(pcr, self, sand, clay, OM, bulk)
    field_adj_dens = temp[0]
    poros_adj_dens = temp[1]
    #-Inverse of B (Eq. 18)
    lamda = (pcr.log10(field_adj_dens) - pcr.log10(dry)) / (pcr.log10(1500) - pcr.log10(33))
    #-Saturated conductivity (matric soil), mm/day (Eq. 16)
    ksat = (1930 * (poros_adj_dens)**(3 - lamda)) * 24
    return ksat

#-function to calculate the wilting point based on previous pedotransfer functions
def Wilt(pcr, self, np):
    #-Determine exponents B and A for logarithmic function
    B = (np.log(1500) - np.log(33)) / (pcr.ln(self.RootFieldMap) - pcr.ln(self.RootDryMap))
    A = pcr.exp(np.log(33) + B * pcr.ln(self.RootFieldMap))

    #-Wilting point based on logarithmic function
    wilt = (100 / A)**(-1/B)
    return wilt

#-init pedotransfer processes
def init(self, pcr, config, np):
    self.RootSandMap = pcr.readmap(self.inpath + config.get('PEDOTRANSFER','RootSandMap')) / 100
    self.RootClayMap = pcr.readmap(self.inpath + config.get('PEDOTRANSFER','RootClayMap')) / 100
    self.RootSiltMap = 1 - self.RootSandMap - self.RootClayMap
    self.RootOMMap = pcr.readmap(self.inpath + config.get('PEDOTRANSFER','RootOMMap'))
    try:
        self.RootBulkMap = config.getfloat('PEDOTRANSFER','RootBulkMap')
    except:
        self.RootBulkMap = pcr.readmap(self.inpath + config.get('PEDOTRANSFER','RootBulkMap'))

    self.SubSandMap = pcr.readmap(self.inpath + config.get('PEDOTRANSFER','SubSandMap')) / 100
    self.SubClayMap = pcr.readmap(self.inpath + config.get('PEDOTRANSFER','SubClayMap')) / 100
    self.SubOMMap = pcr.readmap(self.inpath + config.get('PEDOTRANSFER','SubOMMap'))
    try:
        self.SubBulkMap = config.getfloat('PEDOTRANSFER','SubBulkMap')
    except:
        self.SubBulkMap = pcr.readmap(self.inpath + config.get('PEDOTRANSFER','SubBulkMap'))

    self.RootDryMap = self.pedotransfer.Dry(pcr, self, self.RootSandMap, self.RootClayMap, self.RootOMMap, self.RootBulkMap)
    temp = self.pedotransfer.FieldAdj(pcr, self, self.RootSandMap, self.RootClayMap, self.RootOMMap, self.RootBulkMap)
    self.RootFieldMap = temp[0] * self.RootFieldFrac
    self.RootSatMap = temp[2] * self.RootSatFrac
    self.RootFieldMap = pcr.min(self.RootSatMap - 0.0001, self.RootFieldMap)
    self.RootWiltMap = self.pedotransfer.Wilt(pcr, self, np) * self.RootWiltFrac
    self.RootDryMap = self.RootDryMap * self.RootDryFrac
    self.RootKsat = self.pedotransfer.KSat(pcr, self, self.RootSandMap, self.RootClayMap, self.RootOMMap, self.RootBulkMap) * self.RootKsatFrac
    self.RootDrainVel = self.RootKsat * self.Slope

    temp = self.pedotransfer.FieldAdj(pcr, self, self.SubSandMap, self.SubClayMap, self.SubOMMap, self.SubBulkMap)
    self.SubFieldMap = temp[0]
    self.SubSatMap = temp[2]
    self.SubKsat = self.pedotransfer.KSat(pcr, self, self.SubSandMap, self.SubClayMap, self.SubOMMap, self.SubBulkMap)
