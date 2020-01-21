# Sediment transport module that determines sediment flux and sediment yield at the stations,
# considering reservoir sedimentation when the reservoir module is used
# Copyright (C) 2015-2019 Joris Eekhout / Spanish National Research Council (CEBAS-CSIC)
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

#- Equations for simulating sediment transport
print('Sediment transport module imported')

#-Sediment transport capacity
def TC(self, pcr, runoff):
    #-runoff discharge per unit width in m2/day
    q = (runoff/1000) * pcr.celllength()
    #-transport capacity calculation ton/day
    TC = self.TC_k * q**self.TC_beta * self.Slope**self.TC_gamma
    #-only apply TC to cells with runoff (>0) and in the channel (upstream area > upstream_km2)
    TC_upstream = pcr.ifthen(pcr.pcrand(self.Channel == 1, runoff > 0), TC)
    #-determine average transport capacity per subcatchment
    TCSubcatchment = pcr.areaaverage(TC_upstream, self.subcatchmentRes)
    #-set all other cells to 0
    TCSubcatchment = pcr.cover(TCSubcatchment, 0)
    return TC
    # return TCSubcatchment

#-Sediment transport
def SedTrans(self, pcr, np, sed, TC):
    #-determine sediment transport without reservoirs
    if self.ResFLAG == 0:
        #-rout sediment based on transport capacity
        sedimentFlux = pcr.accucapacityflux(self.FlowDir, sed, TC)
        sedimentYield = pcr.accucapacitystate(self.FlowDir, sed, TC)
        sedDep = sedimentYield

    #-determine sediment transport with reservoirs
    else:
        #-store sed in sedTrans to be used in routing algorithm
        sedTrans = sed

        #-initiate empty map to be used in for-loop
        sedimentYield = self.DEM * 0
        sedimentFlux = self.DEM * 0
        subFinished = self.DEM * 0
        sedDep = self.DEM * 0

        #-increase the transport capacity in the reservoir cells such that all sediment is transported towards the end of the reservoir
        TC = pcr.ifthenelse(pcr.scalar(self.SedReservoirs) > 0, 1e10, TC)

        #-Determine total soil erosion
        sedTotal = np.sum(pcr.pcr2numpy(sed, 0))

        #-Only rout sediment when soil erosion > 0
        if (sedTotal > 0):
            #-loop through the catchments, rout sediment and determine sedimentation in reservoirs based on trapping efficiency
            for step in range(max(self.subcatchmentSteps)+1): #-for-loop through the steps
                #-determine the subcatchments for this step
                subcatchments = np.where(self.subcatchmentSteps == step)[0]

                #-rout sediment based on transport capacity
                sedTransCapFlux = pcr.accucapacityflux(self.FlowDir, sedTrans, TC)
                sedTransCapState = pcr.accucapacitystate(self.FlowDir, sedTrans, TC)

                #-initiate empty map to be used in for-loop
                stepBool = self.DEM * 0

                #-for-loop through the subcatchments per step
                for subcatchment in subcatchments:
                    #-create boolean map with location of reservoir
                    reservoirBool = pcr.scalar(pcr.ifthenelse(self.ResSedID == int(self.subcatchmentOrder[subcatchment]), pcr.scalar(1), pcr.scalar(0)))

                    #-extract routed sediment value at the reservoir from sedTransCapFlux
                    reservoirFluxTC = pcr.ifthen(reservoirBool == 1, sedTransCapFlux)

                    #-store trapped sediment in sedimentYield (multiply routed sediment value with trapping efficiency to be stored in reservoir cell)
                    sedimentYield = pcr.ifthenelse(reservoirBool == 1, reservoirFluxTC * self.TrappingEff, sedimentYield)

                    #-update subFinished and give subcatchment cells value 1
                    subFinished = pcr.ifthenelse(pcr.scalar(self.subcatchmentRes) == int(self.subcatchmentOrder[subcatchment]), pcr.scalar(1), subFinished)

                    #-update sedTrans, set all subcatchment cells to 0
                    sedTrans = pcr.ifthenelse(pcr.scalar(subFinished) == 1, 0, sedTrans)
                    #-add reservoir outflow to cell downstream of reservoir (multiply routed sediment value with outflow efficiency)
                    sedTrans = sedTrans + pcr.upstream(self.FlowDir, pcr.ifthenelse(reservoirBool == 1, reservoirFluxTC * self.OutflowEff, pcr.scalar(0)))

                    #-create boolean map with location of reservoir
                    stepBool = stepBool + pcr.scalar(pcr.ifthenelse(self.subcatchmentRes == int(self.subcatchmentOrder[subcatchment]), self.subcatchmentRes == int(self.subcatchmentOrder[subcatchment]), pcr.boolean(0)))

                # store sedTransCapFlux in sedimentFlux
                sedimentFlux = sedTransCapFlux * subFinished + sedimentFlux

                #-rout sediment based on transport capacity
                sedDep = sedDep + sedTransCapState * stepBool

            #-rout sediment based on transport capacity
            sedTransCapFlux = pcr.accucapacityflux(self.FlowDir, sedTrans, TC)

            # store sedTransCapFlux in sedimentFlux
            sedimentFlux = sedTransCapFlux * (1 - subFinished) + sedimentFlux
    
    return sedimentYield, sedDep, sedimentFlux


#-init processes
def init(self, pcr, config, csv, np):
    #-init processes when reservoir module is used
    if self.ResFLAG == 1:
        #-nominal map with reservoir IDs
        self.ResSedID = pcr.cover(self.ResID, 0)
    else:
        self.ResSedID = self.Locations

    #-determine upstream area map
    self.UpstreamArea = pcr.accuflux(self.FlowDir, 1) * pcr.cellarea() / 10**6

    #-determine upstream area smaller than upstream_km2 and define channel cells based on upstream area
    self.Upstream_km2 = config.getfloat('SEDIMENT_TRANS', 'upstream_km2')
    self.Channel = self.UpstreamArea > self.Upstream_km2

    #-determine average slope per stream based on stream order and averaged per subcatchment
    self.Basin = pcr.subcatchment(self.FlowDir, self.ResSedID)
    self.StreamOrder = pcr.streamorder(self.FlowDir)
    self.Streams = pcr.scalar(self.Channel) * pcr.scalar(self.Basin) * 10 + pcr.scalar(self.Channel) * pcr.scalar(self.StreamOrder)
    self.SlopeStreams = pcr.areaaverage(self.Slope, pcr.nominal(self.Streams))
    self.SlopeStreams = pcr.ifthenelse(self.Channel == 1, self.SlopeStreams, self.Slope)

    #-read transport capacity parameters
    self.TC_beta = config.getfloat('SEDIMENT_TRANS', 'TC_beta')
    self.TC_gamma = config.getfloat('SEDIMENT_TRANS', 'TC_gamma')

    #-init processes when reservoir module is used
    if self.ResFLAG == 1:
        #-nominal map with reservoir IDs and extent
        if self.ETOpenWaterFLAG == 1:
            self.SedReservoirs = pcr.cover(self.openWaterNominal, 0)
        else:
            self.SedReservoirs = pcr.readmap(self.inpath + config.get('RESERVOIR', 'reservoirs'))
        self.SedReservoirs = pcr.cover(self.SedReservoirs, 0)

        #-read table with the trapping efficiency per reservoir
        self.TrapEffTab = self.inpath + config.get('SEDIMENT_TRANS', 'TrapEffTab')
        self.TrappingEff = pcr.cover(pcr.lookupscalar(self.TrapEffTab, self.ResSedID), 0)

        #-construct map where all cells have 1 and only the reservoir cells have trapping efficiency value obtained from the table
        self.OutflowEff = pcr.cover(1-pcr.lookupscalar(self.TrapEffTab, self.ResSedID), 1)

        #-determine subcatchment map
        self.subcatchmentRes = pcr.subcatchment(self.FlowDir, self.ResSedID)

        #-read reservoir order for sediment transport and add the values to self.subcatchmentOrder and self.subcatchmentSteps
        self.ResOrder = config.get('SEDIMENT_TRANS', 'ResOrder')
        self.subcatchmentOrder = []
        self.subcatchmentSteps = []

        #-loop through the rows of the text file
        with open(self.inpath + self.ResOrder, 'rt') as f:
            next(f) # skip headings
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                self.subcatchmentOrder = np.append(self.subcatchmentOrder, int(row[0]))
                self.subcatchmentOrder = self.subcatchmentOrder.astype(np.int)
                self.subcatchmentSteps = np.append(self.subcatchmentSteps, int(row[1]))
                self.subcatchmentSteps = self.subcatchmentSteps.astype(np.int)
        
        #-loop through the steps and define map with step per subcatchment
        self.subcatchmentStepsMap = self.DEM * 0
        for step in range(max(self.subcatchmentSteps)+1):
            #-determine the subcatchments for this step
            subcatchments = np.where(self.subcatchmentSteps == step)[0]

            #-for-loop through the subcatchments per step
            for subcatchment in subcatchments:
                #-update subFinished and give subcatchment cells value 1
                self.subcatchmentStepsMap = pcr.ifthenelse(pcr.scalar(self.subcatchmentRes) == int(self.subcatchmentOrder[subcatchment]), pcr.scalar(step + 1), self.subcatchmentStepsMap)
        #-set all subcatchments with value 0 to the maximum step 
        self.subcatchmentStepsMap = pcr.ifthenelse(self.subcatchmentStepsMap == 0, int(max(self.subcatchmentSteps)+1), self.subcatchmentStepsMap)
        #-map with step value for stations that are not reservoirs
        self.LocationsNoResSteps = pcr.scalar(self.LocationsNoRes) * self.subcatchmentStepsMap

    #-determine roughness factor for mmf only
    if self.SedModel == 2:
        #-Read flag if channels should be excluded from the detachment by runoff calculation
        self.manningChannelsFLAG = config.getint('SEDIMENT_TRANS', 'manningChannelFLAG')

        #-read manning value for channels
        if self.manningChannelsFLAG == 1:
            self.manningChannel = config.getfloat('SEDIMENT_TRANS', 'manningChannel')

        #-Determine flow velocity for transport capacity calculation
        self.n_veg_TC = self.mmf.manningVegetation(self, pcr, self.d_TC, self.Diameter, self.NoElements)
        self.n_veg_TC = pcr.ifthenelse(self.NoVegetation == 1, 0, self.n_veg_TC)
        self.n_veg_TC = pcr.ifthenelse(self.NoErosion == 1, 0, self.n_veg_TC)
        self.n_veg_TC = pcr.ifthenelse(self.n_table > 0, self.n_table, self.n_veg_TC)
        self.n_TC = (self.n_soil**2 + self.n_veg_TC**2)**0.5
        #-set manning value of channels to predetermined value
        if self.manningChannelsFLAG == 1:
            self.n_TC = pcr.ifthenelse(self.Channel == 1, self.manningChannel, self.n_TC)
        self.v_TC = self.mmf.FlowVelocity(self, pcr, self.n_TC, self.d_TC)

        #-Determine flow velocity after harvest, manning for tilled conditions is used
        if self.harvest_FLAG:
            self.n_veg_TC_harvest = self.mmf.manningVegetation(self, pcr, self.d_field, self.Diameter_harvest, self.NoElements_harvest)
            self.n_veg_TC_harvest = pcr.ifthenelse(self.Tillage_harvest == 1, 0, self.n_veg_field_harvest)
            self.n_TC_harvest = (self.n_soil**2 + self.n_veg_TC_harvest**2)**0.5
            #-set manning value of channels to predetermined value
            if self.manningChannelsFLAG == 1:
                self.n_TC_harvest = pcr.ifthenelse(pcr.pcrand(self.Channel == 1, self.n_TC_harvest > 0), self.manningChannel, self.n_TC_harvest)
            self.v_TC_harvest = self.mmf.FlowVelocity(self, pcr, self.n_TC_harvest, self.d_TC)

        #-Determine flow velocity for bare soil conditions (reference conditions)
        self.v_b = self.mmf.FlowVelocity(self, pcr, self.n_bare, self.d_bare)

        #-Determine roughness factor for transport capacity calculation
        self.roughnessFactor = self.v_TC / self.v_b

#-initial conditions sediment transport
def initial(self, pcr, config):
    try:
        self.SYieldR = pcr.readmap(self.inpath + config.get('SEDIMENT_TRANS', 'Sed_init'))
    except:
        self.SYieldR = config.getfloat('SEDIMENT_TRANS', 'Sed_init')

#-dynamic sediment transport processes musle
def dynamic_musle(self, pcr):
    #-transport capacity
    TC = self.sediment_transport.TC(self, pcr, (Q * 3600 * 24) / pcr.cellarea() * 1000)

    #-report the transport capacity per subcatchment
    self.reporting.reporting(self, pcr, 'TC', TC)

    #-determine sediment yield at reservoirs
    tempvar = self.sediment_transport.SedTrans(self, pcr, np, sed, TC)
    sedimentYield = tempvar[0]
    sedDep = tempvar[1]

    #-report the deposition in channel cells
    self.reporting.reporting(self, pcr, 'SedDep', sedDep)

    #-report sediment yield in the reservoirs
    self.reporting.reporting(self, pcr, 'SYieldRA', sedimentYield)

#-dynamic sediment transport processes mmf
def dynamic_mmf(self, pcr, Runoff, np, G):
    #-change the flow factor for harvested areas to actual and tillage conditions
    if self.harvest_FLAG == 1:
        self.roughnessFactorUpdate = pcr.ifthenelse(self.Harvested == 1, self.v_TC_harvest / self.v_b, self.roughnessFactor)
    else:
        self.roughnessFactorUpdate = self.roughnessFactor
    
    #-determine transport capacity
    TC = self.mmf.TransportCapacity(self, pcr, self.roughnessFactorUpdate, self.RootClayMap + self.RootSiltMap + self.RootSandMap, Runoff)

    #-report the transport capacity
    self.reporting.reporting(self, pcr, 'TC', TC)

    #-determine sediment yield at stations
    sedYield, sedDep, sedFlux = self.sediment_transport.SedTrans(self, pcr, np, G * pcr.cellarea() / 1000, TC)

    #-report the sediment deposition by transport capacity (ton/day)
    self.reporting.reporting(self, pcr, 'SedDep', sedDep)

    #-report sediment yield in the stations (ton/day)
    self.reporting.reporting(self, pcr, 'SedYld', sedYield)

    #-report sediment flux in the stations (ton/day)
    self.reporting.reporting(self, pcr, 'SedFlux', sedFlux)