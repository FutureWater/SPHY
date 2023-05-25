# The Spatial Processes in HYdrology (SPHY) model:
# A spatially distributed hydrological model 
# Copyright (C) 2013-2023  FutureWater
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

# Advanced routing that is used for reservoirs or lakes
print('Advanced routing module for lakes and reservoirs imported')

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

#-init advanced routing processes
def init(self, pcr, config):
    self.FlowDir = pcr.readmap(self.inpath + config.get('ROUTING','flowdir'))
    try:
        self.kx = pcr.readmap(self.inpath + config.get('ROUTING','kx'))
    except:
        self.kx = config.getfloat('ROUTING','kx')

#-initial conditions advanced routing
def initial(self, pcr, config):
    #-initial routed total runoff
    try:
        self.QRAold = config.getfloat('ROUT_INIT','QRA_init')
    except:
        try:
            self.QRAold = pcr.readmap(self.inpath + config.get('ROUT_INIT','QRA_init'))
        except:
            self.QRAold = 0
    #-initial routed runoff	for the individual components
    pars = ['RootR', 'RootD', 'Rain', 'Snow', 'Glac', 'Base']
    for i in pars:
        try:
            setattr(self, i + 'RAold', pcr.readmap(self.inpath + config.get('ROUT_INIT', i + 'RA_init')))
            setattr(self, i + 'RA_FLAG', True)
        except:
            try:
                setattr(self, i + 'RAold', config.getfloat('ROUT_INIT', i + 'RA_init'))
                setattr(self, i + 'RA_FLAG', True)
            except:
                setattr(self, i + 'RA_FLAG', False)

    #-initial storage in lakes/reservoirs of individual flow components
    pars = ['RootRRA','RootDRA','RainRA','SnowRA','GlacRA','BaseRA']
    for i in pars:
        column = pars.index(i)  # identify column to be read from lake or reservoir table
        try: #-try to sum the storages read from the lake and reservoir tables if both thse modules are used
            setattr(self, i + 'stor', (pcr.cover(pcr.lookupscalar(LakeStor_Tab, column + 2, self.LakeID), 0) + \
                                    pcr.cover(pcr.lookupscalar(ResStor_Tab, column + 3, self.ResID), 0)) * 10**6)
            if eval('self.' + i + '_FLAG'):
                setattr(self, i + '_FLAG', True)
            else:
                setattr(self, i + '_FLAG', False)
        except:
            try: #-try to read the storages from the lake table
                setattr(self, i + 'stor', pcr.cover(pcr.lookupscalar(LakeStor_Tab, column + 2, self.LakeID), 0) * 10**6)
                if eval('self.' + i + '_FLAG'):
                    setattr(self, i + '_FLAG', True)
                else:
                    setattr(self, i + '_FLAG', False)
            except: #-try to read the storages from the reservoir table
                try:
                    setattr(self, i + 'stor', pcr.cover(pcr.lookupscalar(ResStor_Tab, column + 3, self.ResID), 0) * 10**6)
                    if eval('self.' + i + '_FLAG'):
                        setattr(self, i + '_FLAG', True)
                    else:
                        setattr(self, i + '_FLAG', False)
                except:
                    setattr(self, i + '_FLAG', False)

#-dynamic processes advanced routing
def dynamic(self, pcr, pcrm, config, TotR, ETOpenWater, PrecipTot):
    #-Update storage in lakes/reservoirs (m3) with specific runoff
    self.StorRES = self.StorRES + pcr.ifthenelse(self.QFRAC==0, 0.001 * pcr.cellarea() * TotR, 0)
    OldStorage = self.StorRES
    #-Calculate lake/reservoir outflow volumes
    if self.LakeFLAG ==1 and self.ResFLAG ==1:
        tempvar = self.lakes.UpdateLakeHStore(self, pcr, pcrm)
        LakeLevel = tempvar[0]
        self.StorRES = tempvar[1]
        LakeQ = self.lakes.QLake(self, pcr, LakeLevel)
        ResQ = self.reservoirs.QRes(self, pcr)
        Qout = pcr.ifthenelse(self.ResID != 0, ResQ, pcr.ifthenelse(self.LakeID!=0, LakeQ, 0))
    elif self.LakeFLAG ==1:
        tempvar = self.lakes.UpdateLakeHStore(self, pcr, pcrm)
        LakeLevel = tempvar[0]
        self.StorRES = tempvar[1]
        Qout = self.lakes.QLake(self, pcr, LakeLevel)
    else:
        Qout = self.reservoirs.QRes(self, pcr)

    if self.travelTimeFLAG == 0:
        #-Calculate volume available for routing (=outflow lakes/reservoir + cell specific runoff)
        RunoffVolume = pcr.upstream(self.FlowDir, Qout) + pcr.ifthenelse(self.QFRAC==0, 0, 0.001 * pcr.cellarea() * TotR)
        #-Routing of total flow
        tempvar = self.advanced_routing.ROUT(self, pcr, RunoffVolume, self.QRAold, Qout, self.StorRES)
        self.StorRES = tempvar[0]
        Q = tempvar[1]
        Qin = tempvar[2]
        self.QRAold = Q
    else:
        #-Update channelStorage by adding the total runoff (m3)
        self.channelStorage += pcr.upstream(self.FlowDir, Qout) + pcr.ifthenelse(self.QFRAC==0, 0, TotR * 0.001 * pcr.cellarea())

        #-Determine flow velocity (m/day)
        flowVelocity, hydraulicRadius = self.travel_time_routing.flowVelocity(self, pcr, self.waterDepth)

        #-Determine channel discharge (m3/s) for first iteration step
        Q = pcr.accutraveltimefractionflux(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity), self.QFRAC) / self.dT

        # Only calculate inflow for lake/reservoir cells
        Qin = pcr.ifthenelse(self.QFRAC==0, pcr.upstream(self.FlowDir, Q * self.dT), 0)

        #-Determine flow velocity for first iteration (m/s)
        u = pcr.max(flowVelocity, 1e-0) / self.dT
        
        #-Determine water depth for first iteration (m)
        self.waterDepth = self.travel_time_routing.waterDepth(pcr, Q, u, self.channelDepth, self.channelWidth, self.floodplainWidth)

        #-Determine flow velocity (m/s), water depth (m) and resulting discharge (m3/s) through iteration
        Q, u, hydraulicRadius = self.travel_time_routing.flow_velocity_iteration(self, pcr, Q)

        #-store flow velocity
        self.flowVelocity = u

        #-Only calculate inflow for lake/reservoir cells
        Qin = pcr.ifthenelse(self.QFRAC==0, pcr.upstream(self.FlowDir, Q * self.dT), 0)
        self.StorRES = self.StorRES - Qout + Qin

        #-Determine Q out in m3/s and store in reservoir cells of Q map
        Q = pcr.ifthenelse(self.ResID != 0, Qout / (3600 * 24), Q)

        #-Determine flow velocity in m/day
        flowVelocity = u * self.dT

        #-Update channelStorage after routing (m3)
        self.channelStorage = pcr.accutraveltimefractionstate(self.FlowDir, self.channelStorage, pcr.max(0.0, flowVelocity), self.QFRAC) 

    #-In case open water evaporation is used, correct reeservoir storage for evaporation and precipitation input
    if self.ETOpenWaterFLAG == 1:
        #-determine actual evapotranspiration per reservoir in m3/day
        ETaRES = pcr.ifthenelse(self.StorRES > 0, pcr.min((pcr.areatotal(ETOpenWater * pcr.cellarea() * self.openWaterFrac, self.openWaterNominal) * 0.001), self.StorRES), 0)
        #-Determine total precipitation as input for fraction of open water
        PrecipRES = pcr.ifthenelse(self.StorRES > 0, (pcr.areatotal(PrecipTot * pcr.cellarea() * self.openWaterFrac, self.openWaterNominal) * 0.001), 0)
        #-update storage by subtracting the actual evapotranspiration per reservoir
        self.StorRES = self.StorRES - ETaRES + PrecipRES

    #-Report discharge
    self.reporting.reporting(self, pcr, 'QallRAtot', Q)

    #-report flux in mm
    if self.mm_rep_FLAG == 1 and self.QTOT_mm_FLAG == 1:
        self.QTOTSubBasinTSS.sample(((Q * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) * 1000)
    #-report lake and reservoir waterbalance
    if self.LakeFLAG == 1 and config.getint('REPORTING', 'Lake_wbal') ==1:
        self.LakeInTSS.sample(Qin)
        self.LakeOutTSS.sample(Qout)
        self.LakeStorTSS.sample(self.StorRES)
    if self.ResFLAG == 1 and config.getint('REPORTING', 'Res_wbal') ==1:
        self.ResInTSS.sample(Qin)
        self.ResOutTSS.sample(Qout)
        self.ResStorTSS.sample(self.StorRES)
        if self.ETOpenWaterFLAG:
            self.ResETaTSS.sample(ETaRES)
            self.ResInCalTSS.sample(Qin - ETaRES)

    #-Routing of individual contributers
    #-Snow routing
    if self.SnowRA_FLAG == 1 and self.SnowFLAG == 1:
        self.SnowRAstor = self.SnowRAstor + pcr.ifthenelse(self.QFRAC==0, SnowR * 0.001 * pcr.cellarea(), 0)
        cQfrac = pcr.cover(self.SnowRAstor / OldStorage, 0)
        cQout = cQfrac * Qout
        cRunoffVolume = pcr.upstream(self.FlowDir, cQout) + pcr.ifthenelse(self.QFRAC==0, 0, 0.001 * pcr.cellarea() * SnowR)
        tempvar = self.advanced_routing.ROUT(self, pcr, cRunoffVolume, self.SnowRAold, cQout, self.SnowRAstor)
        self.SnowRAstor = tempvar[0]
        SnowRA = tempvar[1]
        cQin = tempvar[2]
        self.SnowRAold = SnowRA
        self.reporting.reporting(self, pcr, 'SnowRAtot', SnowRA)
        if self.mm_rep_FLAG == 1 and self.QSNOW_mm_FLAG == 1:
            self.QSNOWSubBasinTSS.sample(((SnowRA * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) *1000)
        #-report lake and reservoir waterbalance
        if self.LakeFLAG == 1 and config.getint('REPORTING', 'Lake_wbal') ==1:
            self.LakeSnowInTSS.sample(cQin)
            self.LakeSnowOutTSS.sample(cQout)
            self.LakeSnowStorTSS.sample(self.SnowRAstor)
        if self.ResFLAG == 1 and config.getint('REPORTING', 'Res_wbal') ==1:
            self.ResSnowInTSS.sample(cQin)
            self.ResSnowOutTSS.sample(cQout)
            self.ResSnowStorTSS.sample(self.SnowRAstor)
    #-Rain routing
    if self.RainRA_FLAG == 1:
        self.RainRAstor = self.RainRAstor + pcr.ifthenelse(self.QFRAC==0, RainR * 0.001 * pcr.cellarea(), 0)
        cQfrac = pcr.cover(self.RainRAstor / OldStorage, 0)
        cQout = cQfrac * Qout
        cRunoffVolume = pcr.upstream(self.FlowDir, cQout) + pcr.ifthenelse(self.QFRAC==0, 0, 0.001 * pcr.cellarea() * RainR)
        tempvar = self.advanced_routing.ROUT(self, pcr, cRunoffVolume, self.RainRAold, cQout, self.RainRAstor)
        self.RainRAstor = tempvar[0]
        RainRA = tempvar[1]
        cQin = tempvar[2]
        self.RainRAold = RainRA
        self.reporting.reporting(self, pcr, 'RainRAtot', RainRA)
        if self.mm_rep_FLAG == 1 and self.QRAIN_mm_FLAG == 1:
            self.QRAINSubBasinTSS.sample(((RainRA * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) *1000)
        #-report lake and reservoir waterbalance
        if self.LakeFLAG == 1 and config.getint('REPORTING', 'Lake_wbal') ==1:
            self.LakeRainInTSS.sample(cQin)
            self.LakeRainOutTSS.sample(cQout)
            self.LakeRainStorTSS.sample(self.RainRAstor)
        if self.ResFLAG == 1 and config.getint('REPORTING', 'Res_wbal') ==1:
            self.ResRainInTSS.sample(cQin)
            self.ResRainOutTSS.sample(cQout)
            self.ResRainStorTSS.sample(self.RainRAstor)
    #-Glacier routing
    if self.GlacRA_FLAG == 1 and self.GlacFLAG == 1:
        self.GlacRAstor = self.GlacRAstor + pcr.ifthenelse(self.QFRAC==0, GlacR * 0.001 * pcr.cellarea(), 0)
        cQfrac = pcr.cover(self.GlacRAstor / OldStorage, 0)
        cQout = cQfrac * Qout
        cRunoffVolume = pcr.upstream(self.FlowDir, cQout) + pcr.ifthenelse(self.QFRAC==0, 0, 0.001 * pcr.cellarea() * GlacR)
        tempvar = self.advanced_routing.ROUT(self, pcr, cRunoffVolume, self.GlacRAold, cQout, self.GlacRAstor)
        self.GlacRAstor = tempvar[0]
        GlacRA = tempvar[1]
        cQin = tempvar[2]
        self.GlacRAold = GlacRA
        self.reporting.reporting(self, pcr, 'GlacRAtot', GlacRA)
        if self.mm_rep_FLAG == 1 and self.QGLAC_mm_FLAG == 1:
            self.QGLACSubBasinTSS.sample(((GlacRA * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) *1000)
        #-report lake and reservoir waterbalance
        if self.LakeFLAG == 1 and config.getint('REPORTING', 'Lake_wbal') ==1:
            self.LakeGlacInTSS.sample(cQin)
            self.LakeGlacOutTSS.sample(cQout)
            self.LakeGlacStorTSS.sample(self.GlacRAstor)
        if self.ResFLAG == 1 and config.getint('REPORTING', 'Res_wbal') ==1:
            self.ResGlacInTSS.sample(cQin)
            self.ResGlacOutTSS.sample(cQout)
            self.ResGlacStorTSS.sample(self.GlacRAstor)
    #-Baseflow routing
    if self.BaseRA_FLAG == 1:
        self.BaseRAstor = self.BaseRAstor + pcr.ifthenelse(self.QFRAC==0, self.BaseR * 0.001 * pcr.cellarea(), 0)
        cQfrac = pcr.cover(self.BaseRAstor / OldStorage, 0)
        cQout = cQfrac * Qout
        cRunoffVolume = pcr.upstream(self.FlowDir, cQout) + pcr.ifthenelse(self.QFRAC==0, 0, 0.001 * pcr.cellarea() * self.BaseR)
        tempvar = self.routing.ROUT(self, pcr, cRunoffVolume, self.BaseRAold, cQout, self.BaseRAstor)
        self.BaseRAstor = tempvar[0]
        BaseRA = tempvar[1]
        cQin = tempvar[2]
        self.BaseRAold = BaseRA
        self.reporting.reporting(self, pcr, 'BaseRAtot', BaseRA)
        if self.mm_rep_FLAG == 1 and self.QBASE_mm_FLAG == 1:
            self.QBASESubBasinTSS.sample(((BaseRA * 3600 * 24) / pcr.catchmenttotal(pcr.cellarea(), self.FlowDir)) *1000)
        #-report lake and reservoir waterbalance
        if self.LakeFLAG == 1 and config.getint('REPORTING', 'Lake_wbal') ==1:
            self.LakeBaseInTSS.sample(cQin)
            self.LakeBaseOutTSS.sample(cQout)
            self.LakeBaseStorTSS.sample(self.BaseRAstor)
        if self.ResFLAG == 1 and config.getint('REPORTING', 'Res_wbal') ==1:
            self.ResBaseInTSS.sample(cQin)
            self.ResBaseOutTSS.sample(cQout)
            self.ResBaseStorTSS.sample(self.BaseRAstor)
    
    if self.travelTimeFLAG == 1:
        return Q, u, hydraulicRadius
    else:
        return Q