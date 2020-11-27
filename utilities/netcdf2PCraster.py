# Tool to force the SPHY model with NetCDF files
# Copyright (C) 2018-2019 Joris Eekhout / Spanish National Research Council (CEBAS-CSIC)
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

import subprocess
import netCDF4 as nc 
import os
import numpy as np
# import sys
# np.set_printoptions(threshold=sys.maxsize)
import csv
from scipy.interpolate import griddata
from scipy import spatial
from scipy.spatial import distance
from pyproj import Proj, transform
from math import *

#-file cache to minimize/reduce opening/closing files
filecache = dict()

#-initial processes to determine the x and coordinates of the model grid and netcdf grid
def netcdf2pcrInit(self, pcr, forcing):
    #-define input and ouput projections
    if getattr(self, forcing + 'InProj') == "rotated":
        inProj = Proj(init="epsg:4326")
    else:
        inProj = Proj(init = getattr(self, forcing + 'InProj'))
    outProj = Proj(init = getattr(self, forcing + 'OutProj'))

    #-get the attributes of cloneMap
    attributeClone = getMapAttributesALL(self.clonefile)
    cellsizeClone = attributeClone['cellsize']
    rowsClone = attributeClone['rows']
    colsClone = attributeClone['cols']
    xULClone = attributeClone['xUL']
    yULClone = attributeClone['yUL']

    #-read netcdf file
    f = nc.Dataset(getattr(self, forcing + 'NC'))
    filecache[getattr(self, forcing + 'NC')] = f

    #-get coordinates of upper right and lower left corners of model grid
    xURClone = xULClone + cellsizeClone * colsClone
    yURClone = yULClone
    xLLClone = xULClone
    yLLClone = yULClone - cellsizeClone * rowsClone
    yLRClone = yLLClone
    xLRClone = xURClone

    #-transform coordinates to netcdf projection coordinates
    xULCloneInput,yULCloneInput = transform(outProj, inProj, xULClone, yULClone)
    xURCloneInput,yURCloneInput = transform(outProj, inProj, xURClone, yURClone)
    xLRCloneInput,yLRCloneInput = transform(outProj, inProj, xLRClone, yLRClone)
    xLLCloneInput,yLLCloneInput = transform(outProj, inProj, xLLClone, yLLClone)

    #-determine netcdf cell size and subset coordinates to model domain
    if getattr(self, forcing + 'InProj') == "rotated":
        #-get coordinates from netcdf file
        xrot = f.variables[getattr(self, forcing + 'VarX')][:]
        yrot = f.variables[getattr(self, forcing + 'VarY')][:]

        #-get coordinates of north pole
        npLat = f.variables[getattr(self, forcing + 'VarX')].grid_north_pole_latitude
        npLon = f.variables[getattr(self, forcing + 'VarX')].grid_north_pole_longitude

        #-transform x and y coordinates to grid
        xrot,yrot = np.meshgrid(xrot, yrot)

        #-transform rotated grid to lat,lon-coordinates
        xLatLon = xrot * 0
        yLatLon = yrot * 0
        for idx, row in enumerate(xrot):
            for idy, val in enumerate(row):
                x, y = rotated_grid_transform((xrot[idx, idy], yrot[idx, idy]), 2, (npLon, npLat))
                xLatLon[idx, idy] = x
                yLatLon[idx, idy] = y

        #-transform x,y-coordinates to 2d array
        xyLatLon = [np.array(xLatLon).flatten(), np.array(yLatLon).flatten()]
        xyLatLon = list(map(list, zip(*xyLatLon)))

        #-function to find closest node from list of coordinates
        def closest_node(node, nodes):
            closest_index = distance.cdist([node], nodes).argmin()
            indices = np.where(np.ma.getdata(xLatLon) == nodes[closest_index][0])
            indices = np.array(indices).flatten().astype('int32')
            return indices.tolist()

        #-get indices of corner points clone map from netcdf
        indicesUL = closest_node((xULCloneInput,yULCloneInput), xyLatLon)
        indicesLL = closest_node((xLLCloneInput,yLLCloneInput), xyLatLon)
        indicesUR = closest_node((xURCloneInput,yURCloneInput), xyLatLon)
        indicesLR = closest_node((xLRCloneInput,yLRCloneInput), xyLatLon)

        #-determine indices of the corners of netcdf grid corresponding to model grid (+ buffer)
        xyUL = max(min(indicesUL[0], indicesLL[0], indicesUR[0], indicesLR[0]) - 2, 0)
        xyLL = min(max(indicesUL[0], indicesLL[0], indicesUR[0], indicesLR[0]) + 2, xLatLon.shape[0] - 1)
        xyUR = max(min(indicesUL[1], indicesLL[1], indicesUR[1], indicesLR[1]) - 2, 0)
        xyLR = min(max(indicesUL[1], indicesLL[1], indicesUR[1], indicesLR[1]) + 2, xLatLon.shape[1] - 1)

        #-determine x,y-coordinates corresponding to model grid (+ buffer) from netcdf grid
        x = xLatLon[xyUL:(xyLL + 1), xyUR:(xyLR + 1)]
        y = yLatLon[xyUL:(xyLL + 1), xyUR:(xyLR + 1)]

    else:
        #-get cell size, number of rows and columns and upper left corner coordinates from netcdf grid
        cellsizeInput = f.variables[getattr(self, forcing + 'VarY')][1]- f.variables[getattr(self, forcing + 'VarY')][0]
        cellsizeInput = float(cellsizeInput)

        #-determine x-coordinates corresponding to model grid (+ buffer) from netcdf grid
        xIdxSta = np.argmin(abs(f.variables[getattr(self, forcing + 'VarX')][:] - (min(xULCloneInput, xLLCloneInput) - 2 * cellsizeInput)))
        xIdxEnd = np.argmin(abs(f.variables[getattr(self, forcing + 'VarX')][:] - (max(xURCloneInput, xLRCloneInput) + 2 * cellsizeInput)))
        x = f.variables[getattr(self, forcing + 'VarX')][xIdxSta:(xIdxEnd + 1)]

        #-determine y-coordinates corresponding to model grid (+ buffer) from netcdf grid
        yIdxEnd = np.argmin(abs(f.variables[getattr(self, forcing + 'VarY')][:] - (max(yULCloneInput, yURCloneInput) + 2 * cellsizeInput)))
        yIdxSta = np.argmin(abs(f.variables[getattr(self, forcing + 'VarY')][:] - (min(yLLCloneInput, yLRCloneInput) - 2 * cellsizeInput)))
        y = f.variables[getattr(self, forcing + 'VarY')][yIdxSta:(yIdxEnd + 1)]

        #-transform x and y coordinates to grid
        x,y = np.meshgrid(x, y)

    #-project x and y coordinates to model grid projection
    x,y = transform(inProj, outProj, x, y)

    #-transform x and y coordinates to arrays
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()

    #-determine model grid x and y coordinates and save in grid
    xi = np.arange(xULClone + cellsizeClone * 0.5, (xULClone + cellsizeClone * 0.5) + colsClone * cellsizeClone, cellsizeClone)
    yi = np.arange((yULClone + cellsizeClone * 0.5) - rowsClone * cellsizeClone, yULClone + cellsizeClone * 0.5, cellsizeClone)
    yi = np.flipud(yi)
    xi,yi = np.meshgrid(xi,yi)

    #-determine x,y-coordinates of netcdf file and model domain and indices of netcdf corresponding to model domain
    setattr(self, forcing + 'x', x)
    setattr(self, forcing + 'y', y)
    setattr(self, forcing + 'xi', xi)
    setattr(self, forcing + 'yi', yi)
    if getattr(self, forcing + 'InProj') == "rotated":
        setattr(self, forcing + 'xyUL', xyUL)
        setattr(self, forcing + 'xyLL', xyLL)
        setattr(self, forcing + 'xyUR', xyUR)
        setattr(self, forcing + 'xyLR', xyLR)
    else:
        setattr(self, forcing + 'xIdxSta', xIdxSta)
        setattr(self, forcing + 'xIdxEnd', xIdxEnd)
        setattr(self, forcing + 'yIdxSta', yIdxSta)
        setattr(self, forcing + 'yIdxEnd', yIdxEnd)

#-function to interpolate netcdf gridded data to model grid
def netcdf2pcrDynamic(self, pcr, forcing): #ncFile, varName, dateInput, method, factor, x, y, xi, yi, xIdxSta, xIdxEnd, yIdxSta, yIdxEnd):
    #-read netcdf file
    f = nc.Dataset(getattr(self, forcing + 'NC'))
    filecache[getattr(self, forcing + 'NC')] = f
    
    #-get index from netcdf corresponding with current date
    idx = int(nc.date2index(self.curdate, f.variables['time'], select ='exact'))

    #-get raw netcdf gridded data from netcdf, transform to array and multiply with factor
    if getattr(self, forcing + 'InProj') == "rotated":
        z = f.variables[getattr(self, forcing + 'VarName')][idx, getattr(self, forcing + 'xyUL'):(getattr(self, forcing + 'xyLL') + 1), getattr(self, forcing + 'xyUR'):(getattr(self, forcing + 'xyLR') + 1)]
    else:
        z = f.variables[getattr(self, forcing + 'VarName')][idx, getattr(self, forcing + 'yIdxSta'):(getattr(self, forcing + 'yIdxEnd') + 1), getattr(self, forcing + 'xIdxSta'):(getattr(self, forcing + 'xIdxEnd') + 1)]
    z = np.asarray(z).ravel()
    with np.errstate(invalid='ignore'): # surpress error message when there are already nans in the z array
        z = np.where(z<=-9999, np.nan, z) * getattr(self, forcing + 'Factor')
    
    #-remove nans from arrays
    x = getattr(self, forcing + 'x')[~np.isnan(z)]
    y = getattr(self, forcing + 'y')[~np.isnan(z)]
    z = z[~np.isnan(z)]
    
    #-interpolate with method (linear or cubic)
    zi = griddata((x, y), z, (getattr(self, forcing + 'xi'), getattr(self, forcing + 'yi')), method=getattr(self, forcing + 'Method'))
    zi = np.where(np.isnan(zi), -9999, zi)

    #-convert to PCRaster Python map
    output = pcr.numpy2pcr(pcr.Scalar, zi, -9999)

    return output

#-function to interpolate netcdf gridded data to model grid
def netcdf2pcrTimeIdx(self, pcr, forcing): #ncFile, varName, dateInput, method, factor, x, y, xi, yi, xIdxSta, xIdxEnd, yIdxSta, yIdxEnd):
    #-read netcdf file
    f = nc.Dataset(getattr(self, forcing + 'NC'))
    filecache[getattr(self, forcing + 'NC')] = f
    
    #-get index from netcdf corresponding with current date
    setattr(self, forcing + 'TimeIdx', nc.date2index(self.curdate, f.variables['time'], select ='exact'))

#-function to get map attributes from clone map
def getMapAttributesALL(cloneMap):
    cOut,err = subprocess.Popen(str('mapattr -p %s ' %(cloneMap)), stdout=subprocess.PIPE, stderr=open(os.devnull), shell=True).communicate()
    cellsize = float(cOut.split()[7])
    mapAttr = {'cellsize': float(cellsize)        ,\
               'rows'    : float(cOut.split()[3]) ,\
               'cols'    : float(cOut.split()[5]) ,\
               'xUL'     : float(cOut.split()[17]),\
               'yUL'     : float(cOut.split()[19])}
    co = None; cOut = None; err = None
    del co; del cOut; del err
    return mapAttr 

#-function to read the config and 
def getConfigNetcdf(self, config, forcing, section):
    setattr(self, forcing + 'NC', config.get(section, forcing + 'Netcdf'))
    setattr(self, forcing + 'NetcdfInput', config.get(section, forcing + 'NetcdfInput').split(','))
    netcdfInput = getattr(self, forcing + 'NetcdfInput')
    setattr(self, forcing + 'VarName', netcdfInput[0])
    setattr(self, forcing + 'VarX', netcdfInput[1])
    setattr(self, forcing + 'VarY', netcdfInput[2])
    setattr(self, forcing + 'Method', netcdfInput[3])
    setattr(self, forcing + 'Factor', float(netcdfInput[4]))
    setattr(self, forcing + 'InProj', netcdfInput[5])
    setattr(self, forcing + 'OutProj', netcdfInput[6])


#-function to transform rotated lat-lon to regular lat-lon
def rotated_grid_transform(grid_in, option, SP_coor):
    lon = grid_in[0]
    lat = grid_in[1]

    if option == 2:
        lon = -lon
        lat = -lat

    lon = (lon*pi)/180 # Convert degrees to radians
    lat = (lat*pi)/180

    SP_lon = SP_coor[0]
    SP_lat = SP_coor[1]

    theta = 90+SP_lat # Rotation around y-axis
    phi = SP_lon # Rotation around z-axis

    theta = (theta*pi)/180
    phi = (phi*pi)/180 # Convert degrees to radians

    x = cos(lon)*cos(lat) # Convert from spherical to cartesian coordinates
    y = sin(lon)*cos(lat)
    z = sin(lat)

    if option == 1: # Regular -> Rotated

        x_new = cos(theta)*cos(phi)*x + cos(theta)*sin(phi)*y + sin(theta)*z
        y_new = -sin(phi)*x + cos(phi)*y
        z_new = -sin(theta)*cos(phi)*x - sin(theta)*sin(phi)*y + cos(theta)*z

    else:  # Rotated -> Regular

        phi = -phi
        theta = -theta

        x_new = cos(theta)*cos(phi)*x + sin(phi)*y + sin(theta)*cos(phi)*z
        y_new = -cos(theta)*sin(phi)*x + cos(phi)*y - sin(theta)*sin(phi)*z
        z_new = -sin(theta)*x + cos(theta)*z



    lon_new = atan2(y_new,x_new) # Convert cartesian back to spherical coordinates
    lat_new = asin(z_new)

    lon_new = (lon_new*180)/pi # Convert radians back to degrees
    lat_new = (lat_new*180)/pi

    if option == 1: # Regular -> Rotated
        lon_new = -lon_new
        lat_new = -lat_new
    
    return lon_new, lat_new
