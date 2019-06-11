# The Spatial Processes in HYdrology (SPHY) model:
# A spatially distributed hydrological model that calculates soil-water and
# cryosphere processes on a cell-by-cell basis.
#
# Copyright (C) 2013  FutureWater
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
# Email: info@futurewater.nl

#-Authorship information-###################################################################
__authors__ = "W. Terink, A. Lutz, G. Simons, W. Immerzeel and P. Droogers"
__copyright__ = "FutureWater"
__license__ = "GPL"
__version__ = "2.0"
__email__ = "info@futurewater.nl"
__date__ ='1 January 2017'
############################################################################################

print 'Reservoirs and advanced routing module imported'

#-Update storage (volume available for routing (m3))
def Storage(pcr, storact, baser, rainr, glacr, snowr):
    storact = storact + (0.001 * pcr.cellarea() * (baser + rainr + glacr + snowr))
    return storact

#-Function to calculate the fraction to rout and the updated storage (remainder after routing)
def QFRACSTOR(pcr, storact, lakeid, updatellevel, qh_func, qh_exp_a, qh_exp_b, qh_pol_b, qh_pol_a1, qh_pol_a2, \
                qh_pol_a3, hs_func, hs_exp_a, hs_exp_b, hs_pol_b, hs_pol_a1, hs_pol_a2, hs_pol_a3, \
                sh_func, sh_exp_a, sh_exp_b, sh_pol_b, sh_pol_a1, sh_pol_a2, sh_pol_a3, h, level=False):
    #-buffer old storage
    oldstorage = storact
    if level:  # if a lake level is measured, update the storage with measured lake level
        storact = pcr.ifthenelse(updatellevel==0, storact, pcr.ifthenelse(pcr.defined(h), pcr.ifthenelse(sh_func==1, \
            sh_exp_a * pcr.exp(sh_exp_b * h), pcr.ifthenelse(sh_func==2, (sh_pol_a1 * h) + sh_pol_b , \
            pcr.ifthenelse(sh_func==3, (sh_pol_a2 * h**2) + (sh_pol_a1 * h) + sh_pol_b, pcr.ifthenelse(sh_func==4, \
            (sh_pol_a3 * h**3) + (sh_pol_a2 * h**2) + (sh_pol_a1 * h) + sh_pol_b, storact)))), storact))
        h = pcr.ifthenelse(updatellevel==0, pcr.ifthenelse(hs_func==1, hs_exp_a * pcr.exp(hs_exp_b * storact), \
            pcr.ifthenelse(hs_func==2, (hs_pol_a1 * storact) + hs_pol_b, pcr.ifthenelse(hs_func==3, (hs_pol_a2 * storact**2) \
            + (hs_pol_a1 * storact) + hs_pol_b, (hs_pol_a3 * storact**3) + (hs_pol_a2 * storact**2) + (hs_pol_a1 * storact) \
            + hs_pol_b))), h)
    else: # if no lake level is available, then calculate the h based on storages
        h = pcr.ifthenelse(hs_func==1, hs_exp_a * pcr.exp(hs_exp_b * storact), pcr.ifthenelse(hs_func==2, \
            (hs_pol_a1 * storact) + hs_pol_b, pcr.ifthenelse(hs_func==3, (hs_pol_a2 * storact**2) + (hs_pol_a1 * storact) \
            + hs_pol_b, (hs_pol_a3 * storact**3) + (hs_pol_a2 * storact**2) + (hs_pol_a1 * storact) + hs_pol_b)))
    # Calculate the discharge
    Q = pcr.ifthenelse(qh_func==1, qh_exp_a * pcr.exp(qh_exp_b * h), pcr.ifthenelse(qh_func==2, (qh_pol_a1 * h) + qh_pol_b, \
            pcr.ifthenelse(qh_func==3, (qh_pol_a2 * h**2) + (qh_pol_a1 * h) + qh_pol_b, (qh_pol_a3 * h**3) + (qh_pol_a2 * h**2) \
            + (qh_pol_a1 * h) + qh_pol_b)))
    Q = pcr.max(0, Q)
    Q = Q * 3600 * 24 # convert to m3/d
    fracQ = pcr.min(pcr.max(Q / storact, 0), 1)
    fracQ = pcr.ifthenelse(lakeid != 0, fracQ, 1.0)  # for non-lakes fracQ == 1
    storact = pcr.ifthenelse(lakeid != 0, storact, oldstorage)
    return fracQ, storact
    
#-Function to rout the specific runoff
def ROUT(pcr, storact, fracq, oldq, flowdir, kx):
    S = pcr.accufractionstate(flowdir, storact, fracq)  
    Q = pcr.accufractionflux(flowdir, storact, fracq) 
    #-Convert Q to m3/s
    Q = Q / (24 * 3600)
    Q = (1 - kx) * Q + kx * oldq
    return S, Q
    
    
