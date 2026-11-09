# Morphodynamics module that determines morphological change in rills and channels and
# reservoir sedimentation when the reservoir module is used.
# Copyright (C) 2021-2023 Joris Eekhout / Spanish National Research Council (CEBAS-CSIC)
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

print('Morphodynamics module imported')


#-init morphodynamics processes
def init(self, pcr, config):
    if self.travelTimeFLAG == 1:
        #-read input parameters
        self.input.input(self, config, pcr, 'layerFraction', 'MORPHODYNAMICS', 'layerFraction', 0)
        
        #-set gravel fraction to 0 for hillslope
        self.RootGravelMap = 0

        #-For loop over the sediment classes
        for sedimentClass in self.sedimentClasses:
            #-read channel material percentage
            self.input.input(self, config, pcr, 'channel' + sedimentClass, 'MORPHODYNAMICS', 'channel' + sedimentClass, 0)
            
            #-assign channel material to channels and use soil texture fraction for hillslopes
            setattr(self, "bedFraction" + sedimentClass, pcr.ifthenelse(self.channelHillslope == 1, getattr(self, "channel" + sedimentClass) / 100, getattr(self, "Root" + sedimentClass + "Map")))

            #-if reservoir is used
            if self.ResFLAG == 1:
                #-set channel depth in reservoirs to 0
                setattr(self, "bedFraction" + sedimentClass, pcr.ifthenelse(pcr.scalar(self.Reservoirs) > 0, 0, getattr(self, "bedFraction" + sedimentClass)))

        #-redefine grain sizes and fractions
        d1, f1 = self.deltaClay, self.bedFractionClay
        d2, f2 = self.deltaSilt, self.bedFractionSilt
        d3, f3 = self.deltaSand, self.bedFractionSand
        d4, f4 = self.deltaGravel, self.bedFractionGravel

        #-calculate cumulative fractions
        C1 = f1
        C2 = C1 + f2
        C3 = C2 + f3

        #-prevent division by zero during interpolation
        f2_safe = pcr.max(f2, 1e-6)
        f3_safe = pcr.max(f3, 1e-6)
        f4_safe = pcr.max(f4, 1e-6)

        #-define function to calculate percentile diameter
        def calculate_percentile_diameter(P_target, C1, C2, C3, f2_safe, f3_safe, f4_safe, d1, d2, d3, d4):
            #-interpolate the exact size depending on which bin P_target falls into
            D_Gravel = d3 + ((P_target - C3) / f4_safe) * (d4 - d3)
            D_Sand   = d2 + ((P_target - C2) / f3_safe) * (d3 - d2)
            D_Silt   = d1 + ((P_target - C1) / f2_safe) * (d2 - d1)
            D_Clay   = d1 

            #-map the correct interpolation to the pixel
            D_map = pcr.ifthenelse(C3 < P_target, D_Gravel,
                        pcr.ifthenelse(C2 < P_target, D_Sand,
                        pcr.ifthenelse(C1 < P_target, D_Silt, D_Clay)))
            
            return D_map
        
        #-calculate the D50 and D90 maps
        self.D50 = calculate_percentile_diameter(0.50, C1, C2, C3, f2_safe, f3_safe, f4_safe, d1, d2, d3, d4)
        self.D90 = calculate_percentile_diameter(0.90, C1, C2, C3, f2_safe, f3_safe, f4_safe, d1, d2, d3, d4)

        #-determine the active layer thickness (m)
        global_active_thickness = self.layerFraction * self.D90
        
        #-determine channel (and rill) area (channel width times cell length) (m2)
        channel_area = self.channelWidth * pcr.celllength()
        
        #-determine bulk density (kg/m3)
        bulk_density = self.rho_s * 0.6

        #-convert to mass (tons)
        self.activeLayer = (global_active_thickness * channel_area * bulk_density) * 1e-3

#-initial morphodynamics processes
def initial(self, pcr):
    if self.travelTimeFLAG == 1:
        #-For loop over the sediment classes
        for sedimentClass in self.sedimentClasses:
            #-assign bed thickness to all sediment classes
            setattr(self, "sedimentStore" + sedimentClass, 0)