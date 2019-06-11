# Spatial Processes in HYdrology (SPHY) model
<i>Version 2.2</i>

The Spatial Processes in Hydrology (SPHY) model is a hydrological modeling tool suitable for a wide range of water resource management applications. SPHY is a state-of-the-art, easy to use, robust tool, that can be applied for operational as well as strategic decision support. The SPHY model has been applied and tested in various studies ranging from real-time soil moisture predictions in flat lands, to operational reservoir inflow forecasting applications in mountainous catchments, irrigation scenarios in the Nile Basin, and detailed climate change impact studies in the snow- and glacier-melt dominated the Himalayan region.

<img src="https://github.com/FutureWater/SPHY/blob/SPHY2.0/SPHY_concepts.jpg" alt="SPHY model concepts" height="500" width="500">

<b>Changes with respect to version 2.1</b></br>
Glaciers in version 2.1 and 2.0 were not mass conserving. In previous versions they were implemented as a fixed mass generating glacier melt using a degree day factor. This new release accounts for rainfall and snowfall onto the glacier, accumulation and melt of snow, and redistribution of ice from the accumulation to the ablation zone. These modifications allow glaciers to retreat over time if the melt rate is higher than the accumulation rate. Details and an example application can be found in <a href="https://github.com/FutureWater/SPHY/blob/SPHY2.2/SPHY2.2 mass conserving glacier module.pptx">Concepts and application of the mass conserving glacier module</a>.

<img src="https://github.com/FutureWater/SPHY/blob/SPHY2.2/glacier_mass_balance.jpg" alt="Example of annual glacier mass balance" height="400" width="400">

<b>Documentation</b>
<ul>

<li><a href="http://www.geosci-model-dev.net/8/2009/2015/gmd-8-2009-2015.pdf" target="_blank">Terink, W., A.F. Lutz, G.W.H. Simons, W.W. Immerzeel, P. Droogers. 2015. SPHY v2.0: Spatial Processes in HYdrology. Geoscientific Model Development 8: 2009-2034.</a></li>

<li><a href="https://github.com/FutureWater/SPHY/blob/SPHY2.1/SPHY_manualV6.pdf" target="_blank">Terink, W., A.F. Lutz, W.W. Immerzeel. 2015. SPHY v2.0: Spatial Processes in HYdrology. Model theory, installation, and data preparation. FutureWater Report 142.</a></li>

<li><a href="https://github.com/FutureWater/SPHY/blob/SPHY2.1/SPHY_case_studies.pdf" target="_blank">Terink, W., A.F. Lutz, G.W.H. Simons, W.W. Immerzeel. 2015. SPHY: Spatial Processes in HYdrology. Case-studies for training. FutureWater Report 144.</a></li>

<li><a href="https://github.com/FutureWater/SPHY/blob/SPHY2.1/SPHY_GUIs.pdf" target="_blank">Terink, W., A.F. Lutz, W.W. Immerzeel. 2015. SPHY: Spatial Processes in HYdrology. Graphical User-Interfaces (GUIs). FutureWater Report 143.</a></li>

<li><a href="https://github.com/FutureWater/SPHY/blob/SPHY2.1/SPHY_reservoir_module.pdf" target="_blank">Terink, W., P. Droogers, G.W.H. Simons. 2015. Reservoir module in SPHY. Implemented in SPHY v2.1 FutureWater Report 136.</a></li>

<li><a href="https://github.com/FutureWater/SPHY/blob/SPHY2.2/SPHY2.2 mass conserving glacier module.pptx">Concepts and application of the mass conserving glacier module.</a></li>

</ul>

<b>Copyright</b></br>
Copyright (C) 2018 FutureWater. The Spatial Processes in HYdrology (SPHY) model is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <a href="http://www.gnu.org/licenses/" target="_blank">http://www.gnu.org/licenses/</a>.

Contact:
info@futurewater.nl
