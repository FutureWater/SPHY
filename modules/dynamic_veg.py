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

print('dynamic vegetation module imported')

#-Function that returns leaf area index (LAI)
def LAI(pcr, ndvi, fpar_max, fpar_min, lai_max, ndvi_min, ndvi_max):
    SR = (1 + ndvi)/(1 - ndvi)
    SR_max = (1 + ndvi_max)/(1 - ndvi_max)
    SR_min = (1 + ndvi_min)/(1 - ndvi_min)
    FPAR = pcr.min((SR - SR_min) / (SR_max - SR_min) * (fpar_max - fpar_min), 0.95)
    LAI = lai_max * pcr.log10(1-FPAR)/pcr.log10(1-fpar_max)
    return LAI

#-Function that returns the maximum storage (Smax)
def Smax(self, pcr, LAI):
    if self.Smax_model == 1:
        Smax = 0.935 + 0.498*LAI - 0.00575*(LAI**2)
    elif self.Smax_model == 2:
        Smax = self.Sinf * (1 - pcr.exp(-self.Kappa * LAI))
    return Smax

#-Function that returns crop factor (Kc)
def cropCoefficient(pcr, ndvi, ndvi_min, ndvi_max, kc_min, kc_max):
    Kc = kc_min + (kc_max - kc_min) * pcr.max(pcr.min((ndvi - ndvi_min)/(ndvi_max - ndvi_min), 1), 0)
    return Kc

#-Function that returns the interception, precipitation throughfall, and remaining storage
def Inter_function(pcr, S, Smax, Etr):
    PreT = pcr.max(0, S - Smax)
    S = S - PreT
    Int = pcr.min(1.5 * Etr, S)
    S = S - Int
    return Int, PreT, S 


#-init processes dynamic vegetation
def init(self, pcr, config):
    #-set the ndvi map series to be read
    self.ndvi = self.inpath + config.get('DYNVEG', 'NDVI')
    #-read the vegetation parameters
    LAImax_table = self.inpath + config.get('DYNVEG', 'LAImax')
    self.LAImax = pcr.lookupscalar(LAImax_table, self.LandUse)
    #-read Smax model
    self.Smax_model = config.getint('DYNVEG', 'Smax_model')
    #-Read additional parameters in case of Smax model 2
    if self.Smax_model == 2:
        Smax_table = config.get('DYNVEG', 'Smax_table')
        if Smax_table == "":
            self.Sinf = config.getfloat('DYNVEG', 'Sinf')
            self.Kappa = config.getfloat('DYNVEG', 'Kappa')
        else:
            Smax_table = self.inpath + config.get('DYNVEG', 'Smax_table')
            pcr.setglobaloption('matrixtable')
            self.Sinf = pcr.lookupscalar(Smax_table, 1, self.LandUse)
            self.Kappa = pcr.lookupscalar(Smax_table, 2, self.LandUse)
            pcr.setglobaloption('columntable')

    pars = ['NDVImax','NDVImin','NDVIbase','KCmax','KCmin','FPARmax','FPARmin']
    for i in pars:
        try:
            setattr(self, i, pcr.readmap(self.inpath + config.get('DYNVEG', i)))
        except:
            setattr(self, i, config.getfloat('DYNVEG', i))

#-initial conditions dynamic vegetation
def initial(self, pcr):
    #-initial canopy storage
    self.Scanopy = 0
    #-initial ndvi if first map is not provided
    self.ndviOld = pcr.scalar((self.NDVImax + self.NDVImin)/2)
    #-set initial kc value to one, if kc map is not available for first timestep
    self.KcOld = pcr.scalar(1)

#-dynamic processes dynamic vegetation
def dynamic(self, pcr, pcrm, np, Precip, ETref):
    #-try to read the ndvi map series. If not available, then use ndvi old
    try:
        # ndvi = pcr.readmap(pcrm.generateNameT(self.ndvi, self.counter))
        try:
            ndvi = pcr.readmap(pcrm.generateNameT(self.ndvi, self.counter))
        except:
            ndvi = pcr.readmap(pcrm.generateNameT(self.ndvi, self.curdate.timetuple().tm_yday))
    except:
        ndvi = self.ndviOld
    self.ndviOld = ndvi
    #-fill missing ndvi values with average
    ndvi = pcr.cover(ndvi, pcr.areaaverage(ndvi, self.clone))
    #-set maximum value to 0.999
    ndvi = pcr.min(ndvi, 0.999)
    #-Report ndvi
    self.reporting.reporting(self, pcr, 'NDVI', ndvi)

    # #-calculate the vegetation parameters
    # vegoutput = self.dynamic_veg.Veg_function(self, pcr, ndvi, self.FPARmax, self.FPARmin, self.LAImax, self.NDVImin, self.NDVImax, self.KCmin, self.KCmax)
    # #-Kc
    # self.Kc = vegoutput[0]
    # #-LAI
    # self.LAI = vegoutput[2]

    #-Determine crop coefficient (Kc)
    self.Kc = self.dynamic_veg.cropCoefficient(pcr, ndvi, self.NDVImin, self.NDVImax, self.KCmin, self.KCmax)
    #-Determine LAI
    self.LAI = LAI(pcr, ndvi, self.FPARmax, self.FPARmin, self.LAImax, self.NDVImin, self.NDVImax)
    #-report leaf area index
    self.reporting.reporting(self, pcr, 'LAI', self.LAI)

    #-Update canopy storage
    self.Scanopy = self.Scanopy + Precip
    #-Determine maximum storage (Smax)
    Smax = self.dynamic_veg.Smax(self, pcr, self.LAI)
    #-interception and effective precipitation
    intercep = self.dynamic_veg.Inter_function(pcr, self.Scanopy, Smax, ETref)
    #-interception
    Int = intercep[0]
    Int = Int * (1-self.openWaterFrac)
    #-report interception corrected for fraction
    self.reporting.reporting(self, pcr, 'TotIntF', Int * (1-self.GlacFrac))
    #-effective precipitation
    Precip = intercep[1]
    #-Report effective precipitation corrected for fraction
    self.reporting.reporting(self, pcr, 'TotPrecEF', Precip * (1-self.GlacFrac))
    #-canopy storage
    self.Scanopy = intercep[2]
    #-Report effective precipitation corrected for fraction
    self.reporting.reporting(self, pcr, 'StorCanop', self.Scanopy * (1-self.openWaterFrac))

    return Precip
