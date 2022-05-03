# Conservation module
# Copyright (C) 2017-2019 Joris Eekhout / Spanish National Research Council (CEBAS-CSIC)
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


print('conservation module imported')


#-init processes
def init(self, pcr, config):
    try:
        #-set conservation flag to 1
        self.conservationFLAG = 1

        #-read conservation measures map
        self.input.input(self, config, pcr, 'conservationMeasures', 'EROSION', 'conservationMeasures', 0)

        #-read table with conservation input parameters per conservation measure class
        pcr.setglobaloption('matrixtable')
        CONSERVATION_table = self.inpath + config.get('EROSION', 'CONSERVATION_table')
        self.NoElements_conservation = pcr.lookupscalar(CONSERVATION_table, 1, self.conservationMeasures)
        self.Diameter_conservation = pcr.lookupscalar(CONSERVATION_table, 2, self.conservationMeasures)
        self.n_table_conservation = pcr.lookupscalar(CONSERVATION_table, 3, self.conservationMeasures)
        pcr.setglobaloption('columntable')

        #-Determine flow velocity for conservation measures
        self.n_veg_TC_conservation = self.roughness.manningVegetation(self.d_field, self.Diameter_conservation, self.NoElements_conservation)
        self.n_veg_TC_conservation = pcr.ifthenelse(self.n_table_conservation > 0, self.n_table_conservation, self.n_veg_TC_conservation)
        self.n_TC_conservation = (self.n_soil**2 + self.n_veg_TC_conservation**2)**0.5
        self.v_TC_conservation = self.mmf.FlowVelocity(self, pcr, self.n_TC_conservation, self.d_TC)
    except:
        #-set conservation flag to 0
        self.conservationFLAG = 0