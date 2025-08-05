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

print('glacier module imported')

import numpy as np
#-Function to calculate melt from clean ice or debris covered glaciers
def GlacCDMelt(pcr, temp, ddf, glacfrac):
    glacdmelt = pcr.max(0, temp) * ddf * glacfrac
    return glacdmelt

#-Total glacier melt
def GMelt(glaccimelt, glacdcmelt):
    glacmelt = glaccimelt + glacdcmelt
    return glacmelt

#-Function to calculate runoff from glaciers
def GlacR(glacf, gmelt, glacfrac):
    glacr = glacf * gmelt * glacfrac
    return glacr

#-Function to calculate glacier percolation to groundwater
def GPerc(glacf, gmelt, glacfrac):
    gperc = (1 - glacf) * gmelt * glacfrac
    return gperc


#-init glacier processes
def init(self, pcr, config, pd, np, os):
    #-Table with glaciers properties (U_ID,MOD_ID,GLAC_ID,MOD_H,GLAC_H,DEBRIS,FRAC_GLAC)
    self.GlacTable = pd.read_csv(os.path.join(self.inpath, config.get('GLACIER', 'GlacTable')))
    cols = pd.DataFrame(columns=['MOD_T', 'GLAC_T', 'Prec_GLAC', 'Rain_GLAC', 'Snow_GLAC', 'PotSnowMelt_GLAC',\
                                'ActSnowMelt_GLAC', 'OldSnowStore_GLAC', 'SnowStore_GLAC', 'OldSnowWatStore_GLAC', 'SnowWatStore_GLAC',\
                                'MaxSnowWatStore_GLAC', 'OldTotalSnowStore_GLAC', 'TotalSnowStore_GLAC',\
                                'SnowR_GLAC', 'GlacMelt', 'GlacR', 'GlacPerc'])
    self.GlacTable = pd.concat([self.GlacTable, cols], axis=1).fillna(0)
    #-sort on MOD_ID column
    self.GlacTable.sort_values(by='MOD_ID', inplace=True)
    #-read GlacID flag
    self.GlacID_flag = config.getint('GLACIER','GlacID_flag')
    #-read Glacier variables and split
    self.GlacVars = config.get('GLACIER','GlacVars').split(',')
    #-read memory error variable
    self.GlacID_memerror = config.getint('GLACIER','GlacID_memerror')
    #-Initialize reporting per glacier ID
    if self.GlacID_flag:
        self.glacid = sorted(self.GlacTable['GLAC_ID'].unique())  #-get the unique glacier ids
        if self.GlacID_memerror == 0:  #-with no memory error we can store the pandas dataframes in the computer's memory
            drange = pd.date_range(self.startdate, self.enddate, freq='D')
            for p in self.GlacVars: #-make panda dataframes for each variable to report
                setattr(self, p + '_Table', pd.DataFrame(index = drange, columns=self.glacid,dtype=np.float32))  #-create table for each variable to report
    
    #-Read the glacier maps
    self.ModelID = pcr.pcr2numpy(pcr.readmap(os.path.join(self.inpath, config.get('GLACIER', 'ModelID'))), self.MV) #-Nominal
    self.ModelID_1d = self.ModelID.flatten()    #-1 dim array with model cell IDs
    SelModelID = pd.unique(self.GlacTable['MOD_ID'])  #-model id cells for which to extract temperature, precip, etc. (=cells that have glaciers)
    #-Create keys for glacier cells (index in ModelID_1d where cell has glacier)
    self.GlacierKeys= np.ones(self.ModelID_1d.shape)* self.MV
    n = np.arange(self.ModelID_1d.size)
    iCnt= 0
    for ID in SelModelID:
        if ID in self.ModelID_1d:
            key= n[self.ModelID_1d == ID]
            self.GlacierKeys[iCnt] = key
        iCnt += 1
    self.GlacierKeys = self.GlacierKeys[self.GlacierKeys != self.MV]
    self.GlacierKeys = [int(float(x)) for x in self.GlacierKeys] #-make the index an integer
    #-Read the glacier parameters
    pars = ['DDFG','DDFDG','GlacF']
    for	i in pars:
        try:
            setattr(self, i, pcr.readmap(self.inpath + config.get('GLACIER',i)))
        except:
            setattr(self, i, config.getfloat('GLACIER',i))
    #-Lapse rate for temperature
    self.TLapse_table = pd.read_csv(self.inpath + config.get('GLACIER','TLapse'), header=None, index_col=0, sep=' ', skipinitialspace =True)
    #-Map with glacier IDs
    self.GlacID = pcr.readmap(self.inpath + config.get('GLACIER','GlacID'))
    #-Check if glacier retreat should be calculated
    self.GlacRetreat = config.getint('GLACIER','GlacRetreat')
    if self.GlacRetreat == 1:
        #-Get date for updating the glacier fraction (once a year)
        GlacUpdate = config.get('GLACIER', 'GlacUpdate').split(',')
        self.GlacUpdate = {}
        self.GlacUpdate['month'] = int(GlacUpdate[1])
        self.GlacUpdate['day'] = int(GlacUpdate[0])

#-initial glacier processes
def initial(self, pcr, pd):
    self.GlacTable.set_index('MOD_ID', inplace=True)
    #-initial SnowStore on glacier
    try:
        SnowStore_1d = pcr.pcr2numpy(self.SnowStore, self.MV).flatten()
        S = pd.DataFrame(data={'SnowStore_GLAC': SnowStore_1d[self.GlacierKeys]}, index=self.ModelID_1d[self.GlacierKeys])
        self.GlacTable.update(S)
        S = None; SnowStore_1d = None; del S, SnowStore_1d
    except:
        self.GlacTable['SnowStore_GLAC'] = self.SnowStore
    #-initial SnowWatstore on glacier
    try:
        SnowWatStore_1d = pcr.pcr2numpy(self.SnowWatStore, self.MV).flatten()
        S = pd.DataFrame(data={'SnowWatStore_GLAC': SnowWatStore_1d[self.GlacierKeys]}, index=self.ModelID_1d[self.GlacierKeys])
        self.GlacTable.update(S)
        S = None; SnowWatStore_1d = None; del S, SnowWatStore_1d
    except:
        self.GlacTable['SnowWatStore_GLAC'] = self.SnowWatStore
    self.GlacTable['TotalSnowStore_GLAC']  = self.GlacTable['SnowStore_GLAC'] + self.GlacTable['SnowWatStore_GLAC']
    #-Create table for accumulating glacier melt over certain period
    self.GlacTable['AccuGlacMelt'] = 0.
    #-initialize the glacier fraction
    self.GlacFrac = np.zeros(self.ModelID_1d.shape)
    self.GlacFrac[self.GlacierKeys] = self.GlacTable.groupby(self.GlacTable.index).sum()['FRAC_GLAC']
    self.GlacFrac = self.GlacFrac.reshape(self.ModelID.shape)
    self.GlacFrac = pcr.numpy2pcr(pcr.Scalar, self.GlacFrac, self.MV)
    self.GlacFrac = pcr.ifthen(self.clone, self.GlacFrac)  #-only use values where clone is True
    pcr.report(self.GlacFrac, self.outpath + 'GlacFrac_' + self.curdate.strftime('%Y%m%d') + '.map')
    # 1-D Masks for debris and clean ice
    self.CImask = self.GlacTable['DEBRIS'] == 0
    self.DBmask = np.invert(self.CImask)
    #-If dynamic glacier retreat is on, then make spatial maps of initial glacier area, volume, and depth
    if self.GlacRetreat == 1:
        #-Calculate the  area and depth per model ID
        GlacTable_MODid = self.GlacTable.loc[:,['FRAC_GLAC', 'ICE_DEPTH']]
        GlacTable_MODid['ICE_DEPTH'] = GlacTable_MODid['ICE_DEPTH'] * GlacTable_MODid['FRAC_GLAC']
        GlacTable_MODid['AREA'] = GlacTable_MODid['FRAC_GLAC'] * self.cellArea
        GlacTable_MODid = GlacTable_MODid.groupby(GlacTable_MODid.index).sum()
        #-Report pcraster map of glacier area
        glacArea = np.zeros(self.ModelID_1d.shape)
        glacArea[self.GlacierKeys] = GlacTable_MODid['AREA']
        glacArea = glacArea.reshape(self.ModelID.shape)
        glacArea = pcr.numpy2pcr(pcr.Scalar, glacArea, self.MV)
        glacArea = pcr.ifthen(self.clone, glacArea)  #-only use values where clone is True
        pcr.report(glacArea, self.outpath + 'GlacArea_' + self.curdate.strftime('%Y%m%d') + '.map')
        #-Report pcraster map of glacier depth
        iceDepth = np.zeros(self.ModelID_1d.shape)
        iceDepth[self.GlacierKeys] = GlacTable_MODid['ICE_DEPTH']
        iceDepth = iceDepth.reshape(self.ModelID.shape)
        iceDepth = pcr.numpy2pcr(pcr.Scalar, iceDepth, self.MV)
        iceDepth = pcr.ifthen(self.clone, iceDepth)  #-only use values where clone is True
        pcr.report(iceDepth, self.outpath + 'iceDepth_' + self.curdate.strftime('%Y%m%d') + '.map')
        self.oldIceDepth = iceDepth * 1000 # in mm
        #-Report pcraster map of glacier volume
        pcr.report(iceDepth * glacArea, self.outpath + 'vIce_' + self.curdate.strftime('%Y%m%d') + '.map')
        #-Delete variables that are not needed
        glacArea = None; del glacArea; #iceDepth = None; del iceDepth;
        GlacTable_MODid = None; del GlacTable_MODid

#-dynamic glacier processes
def dynamic(self, pcr, pd, Temp, Precip):
    #-1 dim array of Tavg map
    T_1d = pcr.pcr2numpy(Temp, self.MV).flatten()
    T = pd.DataFrame(data={'MOD_T': T_1d[self.GlacierKeys]}, index=self.ModelID_1d[self.GlacierKeys])
    #-update table with model cel temperature
    self.GlacTable.update(T)
    T = None; T_1d = None; del T, T_1d
    #-lapse temperature for glaciers
    self.GlacTable['GLAC_T'] = self.GlacTable['MOD_T'] - (self.GlacTable['MOD_H'] - self.GlacTable['GLAC_H']) * float(self.TLapse_table.loc[self.curdate.month])
    #-1 dim array of Precip map
    P_1d = pcr.pcr2numpy(Precip, self.MV).flatten()
    P = pd.DataFrame(data={'Prec_GLAC': P_1d[self.GlacierKeys]}, index=self.ModelID_1d[self.GlacierKeys])
    #-update table with model cel precipitation
    self.GlacTable.update(P)
    P = None; P_1d = None; del P, P_1d
    #-Snow and rain differentiation
    self.GlacTable['Rain_GLAC'] = 0; self.GlacTable['Snow_GLAC'] = 0;
    mask = self.GlacTable['GLAC_T'] >= self.Tcrit
    self.GlacTable.loc[mask, 'Rain_GLAC'] = self.GlacTable.loc[mask, 'Prec_GLAC']
    self.GlacTable.loc[np.invert(mask), 'Snow_GLAC'] = self.GlacTable.loc[np.invert(mask), 'Prec_GLAC']
    #-Set the melting temperature (=>0)
    Tmelt = np.maximum(self.GlacTable['GLAC_T'], 0)
    #-Snow melt
    self.GlacTable['PotSnowMelt_GLAC'] = Tmelt * self.DDFS
    self.GlacTable['ActSnowMelt_GLAC'] = np.minimum(self.GlacTable['SnowStore_GLAC'], self.GlacTable['PotSnowMelt_GLAC'])
    #-Update snow store
    self.GlacTable['OldSnowStore_GLAC'] = self.GlacTable['SnowStore_GLAC']
    self.GlacTable['SnowStore_GLAC'] = self.GlacTable['SnowStore_GLAC'] + self.GlacTable['Snow_GLAC'] - self.GlacTable['ActSnowMelt_GLAC']
    self.GlacTable.loc[self.GlacTable['GLAC_T'] < 0., 'SnowStore_GLAC'] = self.GlacTable.loc[self.GlacTable['GLAC_T'] < 0., 'SnowStore_GLAC'] + \
        self.GlacTable.loc[self.GlacTable['GLAC_T'] < 0., 'SnowWatStore_GLAC']
    #-Caclulate the maximum amount of water that can be stored in snowwatstore
    self.GlacTable['MaxSnowWatStore_GLAC'] = self.SnowSC * self.GlacTable['SnowStore_GLAC']
    self.GlacTable['OldSnowWatStore_GLAC'] = self.GlacTable['SnowWatStore_GLAC']
    #-Calculate the actual amount of water stored in snowwatstore
    self.GlacTable['SnowWatStore_GLAC'] = np.minimum(self.GlacTable['MaxSnowWatStore_GLAC'], self.GlacTable['SnowWatStore_GLAC'] +\
        self.GlacTable['ActSnowMelt_GLAC'] + self.GlacTable['Rain_GLAC'])
    self.GlacTable.loc[self.GlacTable['GLAC_T'] < 0., 'SnowWatStore_GLAC'] = 0
    #-Changes in total water storage in snow (SnowStore and SnowWatStore)
    self.GlacTable['OldTotalSnowStore_GLAC'] = self.GlacTable['TotalSnowStore_GLAC']
    self.GlacTable['TotalSnowStore_GLAC'] = self.GlacTable['SnowStore_GLAC'] + self.GlacTable['SnowWatStore_GLAC']
    #-Snow runoff
    mask = (self.GlacTable['SnowWatStore_GLAC'] == self.GlacTable['MaxSnowWatStore_GLAC']) & (self.GlacTable['OldSnowStore_GLAC'] > 0.) #-mask with true where SnowWatStore_GLAC == MaxSnowWatStore_GLAC
    self.GlacTable.loc[mask,'SnowR_GLAC'] = self.GlacTable.loc[mask, 'ActSnowMelt_GLAC'] + self.GlacTable.loc[mask, 'Rain_GLAC'] - \
        (self.GlacTable.loc[mask, 'SnowWatStore_GLAC'] - self.GlacTable.loc[mask, 'OldSnowWatStore_GLAC'])
    self.GlacTable.loc[np.invert(mask), 'SnowR_GLAC'] = 0
    mask = None; del mask

    #-Glacier melt
    self.GlacTable['GlacMelt'] = 0   #-first set to 0 then update hereafter
    #-Masks for full glacier melt (=no snow melt in timestep) and partial glacier melt (=where snowpack is fully melted within timestep)
    partialMelt = (self.GlacTable['OldSnowStore_GLAC'] > 0.) & (self.GlacTable['SnowStore_GLAC'] == 0)  #-mask with true where snowpack has melted within timestep: for these cells also partial glacier melt
    fullMelt = (self.GlacTable['OldSnowStore_GLAC'] == 0) & (self.GlacTable['SnowStore_GLAC'] == 0)
    #-Melt from Clean Ice Glaciers
    mask = (partialMelt & self.CImask)  #-mask for Clean Ice glacier and partial melt
    self.GlacTable.loc[mask, 'GlacMelt'] = np.maximum(self.GlacTable.loc[mask, 'GLAC_T'] - (self.GlacTable.loc[mask, 'OldSnowStore_GLAC'] / self.DDFS), 0) * self.DDFG
    mask = (fullMelt & self.CImask)  #-mask for Clean Ice glacier and full melt
    self.GlacTable.loc[mask, 'GlacMelt'] = Tmelt.loc[mask] * self.DDFG
    #-Melt from Debris Covered Glaciers
    mask = (partialMelt & self.DBmask)  #-mask for Debris covered Glacier and partial melt
    self.GlacTable.loc[mask, 'GlacMelt'] = np.maximum(self.GlacTable.loc[mask, 'GLAC_T'] - (self.GlacTable.loc[mask, 'OldSnowStore_GLAC'] / self.DDFS), 0) * self.DDFDG
    mask = (fullMelt & self.DBmask)  #-mask for Debris covered Glacier and full melt
    self.GlacTable.loc[mask, 'GlacMelt'] = Tmelt.loc[mask] * self.DDFDG
    #-Accumulate glacier melt
    self.GlacTable['AccuGlacMelt'] = self.GlacTable['AccuGlacMelt'] + self.GlacTable['GlacMelt']
    #-Glacier runoff
    mask = self.GlacTable['OldSnowStore_GLAC'] == 0  #-only add rain to glacmelt when there was no snowpack at beginning to time-step
    #-Glacier runoff consisting of melt and rainfall
    self.GlacTable.loc[mask, 'GlacR'] = self.GlacF * (self.GlacTable.loc[mask, 'GlacMelt'] + self.GlacTable.loc[mask, 'Rain_GLAC'])
    #-Glacier percolation consisting of melt and rainfall
    self.GlacTable.loc[mask, 'GlacPerc'] = (1-self.GlacF) * (self.GlacTable.loc[mask, 'GlacMelt'] + self.GlacTable.loc[mask, 'Rain_GLAC'])
    #-Glacier runoff consisting of melt only
    self.GlacTable.loc[np.invert(mask), 'GlacR'] = self.GlacF * self.GlacTable.loc[np.invert(mask), 'GlacMelt']
    #-Glacier percolation consisting of melt only
    self.GlacTable.loc[np.invert(mask), 'GlacPerc'] = (1-self.GlacF) * self.GlacTable.loc[np.invert(mask), 'GlacMelt']
    mask = None; del mask

    #-Aggregation for model grid ID
    GlacTable_MODid = self.GlacTable.loc[:,['Rain_GLAC', 'Snow_GLAC', 'ActSnowMelt_GLAC', 'SnowStore_GLAC',\
                            'SnowWatStore_GLAC', 'TotalSnowStore_GLAC', 'SnowR_GLAC', 'GlacMelt', 'GlacR', 'GlacPerc']]
    GlacTable_MODid = GlacTable_MODid.multiply(self.GlacTable['FRAC_GLAC'], axis='index') #-Multiply with the glacier fraction
    GlacTable_MODid = GlacTable_MODid.groupby(GlacTable_MODid.index).sum() #-Summarize by model ID
    #-report back to model ID
    #-Rainfall on glacier
    Rain_GLAC = np.zeros(self.ModelID_1d.shape)
    Rain_GLAC[self.GlacierKeys] = GlacTable_MODid['Rain_GLAC']
    Rain_GLAC = Rain_GLAC.reshape(self.ModelID.shape)
    Rain_GLAC = pcr.numpy2pcr(pcr.Scalar, Rain_GLAC, self.MV)
    #-Snowfall on glacier
    Snow_GLAC = np.zeros(self.ModelID_1d.shape)
    Snow_GLAC[self.GlacierKeys] = GlacTable_MODid['Snow_GLAC']
    Snow_GLAC = Snow_GLAC.reshape(self.ModelID.shape)
    Snow_GLAC = pcr.numpy2pcr(pcr.Scalar, Snow_GLAC, self.MV)
    #-Act snowmelt from glacier
    ActSnowMelt_GLAC = np.zeros(self.ModelID_1d.shape)
    ActSnowMelt_GLAC[self.GlacierKeys] = GlacTable_MODid['ActSnowMelt_GLAC']
    ActSnowMelt_GLAC = ActSnowMelt_GLAC.reshape(self.ModelID.shape)
    ActSnowMelt_GLAC = pcr.numpy2pcr(pcr.Scalar, ActSnowMelt_GLAC, self.MV)
#             #-Snowstore on glacier
#             SnowStore_GLAC = np.zeros(self.ModelID_1d.shape)
#             SnowStore_GLAC[self.GlacierKeys] = GlacTable_MODid['SnowStore_GLAC']
#             SnowStore_GLAC = SnowStore_GLAC.reshape(self.ModelID.shape)
#             SnowStore_GLAC = pcr.numpy2pcr(pcr.Scalar, SnowStore_GLAC, self.MV)
#             #-SnowWatStore on glacier
#             SnowWatStore_GLAC = np.zeros(self.ModelID_1d.shape)
#             SnowWatStore_GLAC[self.GlacierKeys] = GlacTable_MODid['SnowWatStore_GLAC']
#             SnowWatStore_GLAC = SnowWatStore_GLAC.reshape(self.ModelID.shape)
#             SnowWatStore_GLAC = pcr.numpy2pcr(pcr.Scalar, SnowWatStore_GLAC, self.MV)
    #-TotalSnowStore on glacier
    self.TotalSnowStore_GLAC = np.zeros(self.ModelID_1d.shape)
    self.TotalSnowStore_GLAC[self.GlacierKeys] = GlacTable_MODid['TotalSnowStore_GLAC']
    self.TotalSnowStore_GLAC = self.TotalSnowStore_GLAC.reshape(self.ModelID.shape)
    self.TotalSnowStore_GLAC = pcr.numpy2pcr(pcr.Scalar, self.TotalSnowStore_GLAC, self.MV)
    #-SnowR from glacier
    SnowR_GLAC = np.zeros(self.ModelID_1d.shape)
    SnowR_GLAC[self.GlacierKeys] = GlacTable_MODid['SnowR_GLAC']
    SnowR_GLAC = SnowR_GLAC.reshape(self.ModelID.shape)
    SnowR_GLAC = pcr.numpy2pcr(pcr.Scalar, SnowR_GLAC, self.MV)
    #-Glacier melt
    GlacMelt = np.zeros(self.ModelID_1d.shape)
    GlacMelt[self.GlacierKeys] = GlacTable_MODid['GlacMelt']
    GlacMelt = GlacMelt.reshape(self.ModelID.shape)
    GlacMelt = pcr.numpy2pcr(pcr.Scalar, GlacMelt, self.MV)
    #-Report glacier melt
    self.reporting.reporting(self, pcr, 'TotGlacMelt', GlacMelt)
    if self.mm_rep_FLAG == 1 and (self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1):
        self.GMeltSubBasinTSS.sample(pcr.catchmenttotal(GlacMelt, self.FlowDir) / pcr.catchmenttotal(1, self.FlowDir))
    #-Glacier runoff
    GlacR = np.zeros(self.ModelID_1d.shape)
    GlacR[self.GlacierKeys] = GlacTable_MODid['GlacR']
    GlacR = GlacR.reshape(self.ModelID.shape)
    GlacR = pcr.numpy2pcr(pcr.Scalar, GlacR, self.MV)

    #-Report glacier runoff
    self.reporting.reporting(self, pcr, 'TotGlacR', GlacR)
    #-Glacier percolation
    GlacPerc = np.zeros(self.ModelID_1d.shape)
    GlacPerc[self.GlacierKeys] = GlacTable_MODid['GlacPerc']
    GlacPerc = GlacPerc.reshape(self.ModelID.shape)
    GlacPerc = pcr.numpy2pcr(pcr.Scalar, GlacPerc, self.MV)
    #-Report glacier percolation to groundwater
    self.reporting.reporting(self, pcr, 'TotGlacPerc', GlacPerc)
    GlacTable_MODid = None; del GlacTable_MODid
    
    return Rain_GLAC, Snow_GLAC, ActSnowMelt_GLAC, SnowR_GLAC, GlacMelt, GlacPerc, GlacR


def dynamic_reporting(self, pcr, pd, np):
    #-Reporting for each glacier ID
    if self.GlacID_flag:
        GlacTable_GLACid = self.GlacTable.loc[:, self.GlacVars]

        #-Muliply with fraction, summarize, and divide by total fraction to get glacier ID average
        GlacTable_GLACid = GlacTable_GLACid.multiply(self.GlacTable['FRAC_GLAC'], axis='index')  #-Multiply with fraction
        GlacTable_GLACid['GLAC_ID'] = self.GlacTable['GLAC_ID']; GlacTable_GLACid['FRAC_GLAC'] = self.GlacTable['FRAC_GLAC']  #-Add GLAC_ID column and FRAC_GLAC column
        GlacTable_GLACid = GlacTable_GLACid.groupby('GLAC_ID').sum() #-Summarize by Glacier ID
        FracSum = GlacTable_GLACid['FRAC_GLAC'] #-Get summed glacier fraction
        GlacTable_GLACid = GlacTable_GLACid.div(FracSum, axis='index') #-Divide by total fraction for glacier weighted average (frac==1)
        GlacTable_GLACid['FRAC_GLAC'] = FracSum; FracSum = None; del FracSum
        GlacTable_GLACid = GlacTable_GLACid.transpose()  #-Transpose the glacier id table (ID as columns and vars as index)
        GlacTable_GLACid.fillna(0., inplace=True);

        #-Fill Glacier variable tables for reporting
        if self.GlacID_memerror == 0:
            for v in self.GlacVars:
                vv = getattr(self, v + '_Table'); vv.loc[self.curdate,:] = GlacTable_GLACid.loc[v,:].astype(np.float32)
            v = None; vv = None; del v, vv; GlacTable_GLACid = None; del GlacTable_GLACid
            if self.curdate == self.enddate: #-do the reporting at the final model time-step
                for v in self.GlacVars:
                    eval('self.' + v + '_Table.to_csv("'  + self.outpath + v + '.csv")')
        else:
            df = pd.DataFrame(columns=self.glacid, dtype=np.float32)
            for v in self.GlacVars:
                #df = pd.DataFrame(columns=self.glacid, dtype=float)
                df.loc[self.curdate,:] = GlacTable_GLACid.loc[v,:].astype(np.float32)
                if self.curdate == self.startdate:
                    df.to_csv(self.outpath + v + '.csv', mode='w')
                else:
                    df.to_csv(self.outpath + v + '.csv', mode='a', header=False) #-no header for time-step > 1
            df = None; del df; GlacTable_GLACid = None; del GlacTable_GLACid

    #-Check if glacier retreat should be calculated
    if self.GlacRetreat == 1 and self.curdate.month == self.GlacUpdate['month'] and self.curdate.day == self.GlacUpdate['day']:
        self.dateAfterUpdate = self.curdate + self.datetime.timedelta(days=1)
        #-Create a table with fields used for updating the glacier fraction at defined update date
        GlacFracTable = self.GlacTable.loc[:,['U_ID', 'GLAC_ID','FRAC_GLAC','ICE_DEPTH']]
        #-Set the initial ice volumes and snow store
        GlacFracTable['V_ice_t0'] = GlacFracTable['FRAC_GLAC'] * GlacFracTable['ICE_DEPTH'] * self.cellArea
        GlacFracTable['TotalSnowStore_GLAC'] = self.GlacTable['TotalSnowStore_GLAC']
        GlacFracTable['AccuGlacMelt'] = self.GlacTable['AccuGlacMelt']
        GlacFracTable['dMelt'] = GlacFracTable['TotalSnowStore_GLAC'] - GlacFracTable['AccuGlacMelt']
        GlacFracTable['dMelt'] = GlacFracTable['dMelt'] / 1000 * GlacFracTable['FRAC_GLAC'] * self.cellArea  #-convert to m3
        #-Drop unnecessary columns
        GlacFracTable.drop(['TotalSnowStore_GLAC', 'AccuGlacMelt'], axis=1, inplace=True)

        #-Mask to determine ablation and accumulation UIDs
        ablMask = GlacFracTable['dMelt'] < 0.
        #-Set the ablation and accumulation in the corresponding fields
        GlacFracTable['Accumulation'] = 0.
        GlacFracTable['Ablation'] = 0.
        GlacFracTable.loc[np.invert(ablMask),'Accumulation'] = GlacFracTable.loc[np.invert(ablMask), 'dMelt']
        GlacFracTable.loc[ablMask,'Ablation'] = GlacFracTable.loc[ablMask, 'dMelt']
        #-Set the ice volumes for the ablation UIDs
        GlacFracTable['V_ice_ablation'] = 0.
        GlacFracTable.loc[ablMask,'V_ice_ablation'] = GlacFracTable.loc[ablMask,'V_ice_t0']
        #-Calculate totals per Glacier ID
        GlacID_grouped = GlacFracTable.groupby('GLAC_ID').sum()
        GlacID_grouped = GlacID_grouped.loc[:,['dMelt', 'Accumulation', 'V_ice_ablation']]

        #-Calculate total delta Melt (dMelt), accumulation, and ice volumes (of ablation cells) for each glacier ID
        GlacFracTable['dMelt_group'] = 0.
        GlacFracTable['Accumulation_group'] = 0.
        GlacFracTable['V_ice_ablation_group'] = 0.

        for index, row in GlacID_grouped.iterrows():
            GlacFracTable.loc[GlacFracTable['GLAC_ID'] == index,'dMelt_group'] = row['dMelt']
            GlacFracTable.loc[GlacFracTable['GLAC_ID'] == index,'Accumulation_group'] = row['Accumulation']
            GlacFracTable.loc[GlacFracTable['GLAC_ID'] == index,'V_ice_ablation_group'] = row['V_ice_ablation']
        #-Remove GlacID_grouped table
        GlacID_grouped = None; del GlacID_grouped
        #-Mask for determining if redistribution is negative (remove ice from accumulation cells) or positive (add ice to ablation cells)
        negDistMask = (GlacFracTable['dMelt_group'] < 0.) & (GlacFracTable['dMelt'] >= 0.)
        posDistMask = (GlacFracTable['dMelt_group'] < 0.) & (GlacFracTable['dMelt'] < 0.)
        #-Calculate the ice redistribution
        GlacFracTable['Ice_redist'] = 0.
        GlacFracTable.loc[negDistMask,'Ice_redist'] = -GlacFracTable.loc[negDistMask, 'Accumulation']
        GlacFracTable.loc[posDistMask,'Ice_redist'] = GlacFracTable.loc[posDistMask,'V_ice_ablation'] / GlacFracTable.loc[posDistMask,'V_ice_ablation_group'] *\
            GlacFracTable.loc[posDistMask,'Accumulation_group']
        #-Update ice volume
        GlacFracTable['V_ice_t1'] = GlacFracTable['V_ice_t0'] + GlacFracTable['Accumulation'] + GlacFracTable['Ablation'] + GlacFracTable['Ice_redist']
        #-Remove distribution masks
        negDistMask = None; del negDistMask; posDistMask = None; del posDistMask
        #-Remove unnecessary columns
        GlacFracTable.drop(['dMelt', 'Accumulation', 'Ablation','V_ice_ablation','dMelt_group',\
                'Accumulation_group', 'V_ice_ablation_group', 'Ice_redist'], axis=1, inplace=True)
        #-Calculate where updated ice volume becomes negative and postitive
        GlacFracTable['V_ice_negative'] = np.minimum(0., GlacFracTable['V_ice_t1'])
        GlacFracTable['V_ice_positive'] = np.maximum(0., GlacFracTable['V_ice_t1'])
        #-Calculate totals per Glacier ID
        GlacID_grouped = GlacFracTable.groupby('GLAC_ID').sum()
        GlacID_grouped = GlacID_grouped.loc[:,['V_ice_negative', 'V_ice_positive']]
        GlacFracTable['V_ice_negative_group'] = 0.
        GlacFracTable['V_ice_positive_group'] = 0.
        for index, row in GlacID_grouped.iterrows():
            GlacFracTable.loc[GlacFracTable['GLAC_ID'] == index,'V_ice_negative_group'] = row['V_ice_negative']
            GlacFracTable.loc[GlacFracTable['GLAC_ID'] == index,'V_ice_positive_group'] = row['V_ice_positive']
        #-Remove GlacID_grouped table
        GlacID_grouped = None; del GlacID_grouped
        #-Calculate the ice redistribution
        GlacFracTable['Ice_redist'] = GlacFracTable['V_ice_positive'] / GlacFracTable['V_ice_positive_group'] * GlacFracTable['V_ice_negative_group']
        GlacFracTable['Ice_redist'].fillna(0., inplace=True)
        #-Remove unnecessary columns
        GlacFracTable.drop(['V_ice_negative', 'V_ice_positive', 'V_ice_negative_group', 'V_ice_positive_group'], axis=1, inplace=True)
        #-Update ice volume
        GlacFracTable['V_ice_t2'] = np.maximum(0., GlacFracTable['V_ice_t1'] + GlacFracTable['Ice_redist'])
        #-Update ice thickness
        GlacFracTable['ICE_DEPTH_new'] = GlacFracTable['V_ice_t2'] / (GlacFracTable['FRAC_GLAC'] * self.cellArea)
        noIceMask = (GlacFracTable['ICE_DEPTH_new'] <=0.)  # it can be that melt is greater than total available ice: that results in balance error
        #-Update glacier fraction
        GlacFracTable['FRAC_GLAC_new'] = GlacFracTable['FRAC_GLAC']
        GlacFracTable.loc[noIceMask,'FRAC_GLAC_new'] = 0.
        noIceMask = None; del noIceMask

        #-Update the glactable
        self.GlacTable['FRAC_GLAC'] = GlacFracTable['FRAC_GLAC_new']
        self.GlacTable['ICE_DEPTH'] = GlacFracTable['ICE_DEPTH_new']
        self.GlacTable['SnowStore_GLAC'] = 0.
        self.GlacTable['SnowWatStore_GLAC'] = 0.
        self.GlacTable['TotalSnowStore_GLAC'] = 0.

        #-remove SnowWatStore_GLAC from total snowstore
        self.TotalSnowStore = self.TotalSnowStore - self.TotalSnowStore_GLAC

        #-Set accumulated glacier melt to zero as initial condition for next period
        self.GlacTable['AccuGlacMelt'] = 0.
        #-Remove the GlacFracTable
        GlacFracTable = None; del GlacFracTable

    #-Write glacier table, spatial maps of glacier fraction, ice depth, and average ice depth per glacier on day after update day
    if self.GlacRetreat == 1 and self.curdate == self.dateAfterUpdate:
        #-Write glac table to csv (can be used as initial setting for new run)
        glac_csv = self.GlacTable.loc[:,['U_ID','GLAC_ID','MOD_H','GLAC_H','DEBRIS','FRAC_GLAC','ICE_DEPTH']]
        glac_csv.insert(1, 'MOD_ID', self.GlacTable.index)
        glac_csv.to_csv(self.outpath + 'glacTable_' + self.curdate.strftime('%Y%m%d') + '.csv', index=False)
        glac_csv = None; del glac_csv

        #-Calculate average model ID ice depth and total model glacier fraction
        GlacTable_MODid = self.GlacTable.loc[:,['FRAC_GLAC', 'ICE_DEPTH']]
        GlacTable_MODid['ICE_DEPTH'] = GlacTable_MODid['ICE_DEPTH'] * GlacTable_MODid['FRAC_GLAC']
        GlacTable_MODid['AREA'] = GlacTable_MODid['FRAC_GLAC'] * self.cellArea
        GlacTable_MODid = GlacTable_MODid.groupby(GlacTable_MODid.index).sum()
        GlacTable_MODid.fillna(0., inplace=True)

        #-Report updated glacier fraction map
        self.GlacFrac = np.zeros(self.ModelID_1d.shape)
        self.GlacFrac[self.GlacierKeys] = GlacTable_MODid['FRAC_GLAC']
        self.GlacFrac = self.GlacFrac.reshape(self.ModelID.shape)
        self.GlacFrac = pcr.numpy2pcr(pcr.Scalar, self.GlacFrac, self.MV)
        self.GlacFrac = pcr.ifthen(self.clone, self.GlacFrac)  #-only use values where clone is True
        pcr.report(self.GlacFrac, self.outpath + 'GlacFrac_' + self.curdate.strftime('%Y%m%d') + '.map')

        #-Report pcraster map of glacier area
        glacArea = np.zeros(self.ModelID_1d.shape)
        glacArea[self.GlacierKeys] = GlacTable_MODid['AREA']
        glacArea = glacArea.reshape(self.ModelID.shape)
        glacArea = pcr.numpy2pcr(pcr.Scalar, glacArea, self.MV)
        glacArea = pcr.ifthen(self.clone, glacArea)  #-only use values where clone is True
        pcr.report(glacArea, self.outpath + 'GlacArea_' + self.curdate.strftime('%Y%m%d') + '.map')
        #-Report pcraster map of glacier depth
        iceDepth = np.zeros(self.ModelID_1d.shape)
        iceDepth[self.GlacierKeys] = GlacTable_MODid['ICE_DEPTH']
        iceDepth = iceDepth.reshape(self.ModelID.shape)
        iceDepth = pcr.numpy2pcr(pcr.Scalar, iceDepth, self.MV)
        iceDepth = pcr.ifthen(self.clone, iceDepth)  #-only use values where clone is True
        pcr.report(iceDepth, self.outpath + 'iceDepth_' + self.curdate.strftime('%Y%m%d') + '.map')
        #-Report pcraster map of glacier volume
        pcr.report(iceDepth * glacArea, self.outpath + 'vIce_' + self.curdate.strftime('%Y%m%d') + '.map')
        #-Delete variables that are not needed
        glacArea = None; del glacArea; iceDepth = None; del iceDepth; GlacTable_MODid = None; del GlacTable_MODid
