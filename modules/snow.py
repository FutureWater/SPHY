# The Spatial Processes in HYdrology (SPHY) model:
# A spatially distributed hydrological model
# Copyright (C) 2013-2026  FutureWater
# Email: sphy@futurewater.nl
#
# Authors (alphabetical order):
# P. Droogers, J. Eekhout, A. Fernandez-Rodriguez, W. Immerzeel, S. Khanal, A. Lutz, T. Schults, G. Simons, W. Terink.
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


print("snow module imported")


# -Function to calculate the cold/warm fractions from sinusoidal hourly temperature
def ColdWarmFractions(pcr, tempmin, tempmax):
    """Compute fraction of the day below/above 0C using sinusoidal interpolation.
    T(h) = (Tmax+Tmin)/2 + (Tmax-Tmin)/2 * cos(pi * h / 12), h=1..24
    Returns (frac_cold, frac_warm) where frac_cold + frac_warm = 1
    """
    t_mid = (tempmax + tempmin) / 2
    t_amp = (tempmax - tempmin) / 2
    hours_warm = pcr.scalar(0)
    for ij in range(1, 25, 1):
        t_hour = t_mid + t_amp * pcr.cos(3.1415 * ij / 12)
        hours_warm = hours_warm + pcr.ifthenelse(
            t_hour >= 0, pcr.scalar(1), pcr.scalar(0)
        )
    frac_warm = hours_warm / 24
    frac_cold = 1 - frac_warm
    return frac_cold, frac_warm


# -Function to calculate the potential snow melt (uses positive degree-hours)
def PotSnowMelt(pcr, tempmin, tempmax, ddfs):
    """T(h) = (Tmax+Tmin)/2 + (Tmax-Tmin)/2 * cos(pi * h / 12)"""
    t_mid = (tempmax + tempmin) / 2
    t_amp = (tempmax - tempmin) / 2
    thour = pcr.scalar(0)
    for ij in range(1, 25, 1):
        th_max = pcr.max(0, t_mid + t_amp * pcr.cos(3.1415 * ij / 12))
        thour = thour + th_max
    melt = thour * ddfs / 24
    return melt


# -Function to calculate the actual snow melt
def ActSnowMelt(pcr, snowstore, potmelt):
    melt = pcr.min(snowstore, potmelt)
    return melt


# -Function that updates the snow storage (binary refreeze based on Tavg)
def SnowStoreUpdate(pcr, snowstore, snow, actmelt, temp, tcrit, snowwatstore):
    # Refreeze: if daily average temp < Tcrit, all liquid water refreezes
    # This includes both existing SWS and freshly melted water (ActMelt may be >0
    # on cold-average days when TempMax >= 0, due to sinusoidal hourly melt).
    # Without including ActMelt, the melted water disappears from the water balance.
    refreeze = pcr.ifthenelse(temp < tcrit, snowwatstore + actmelt, pcr.scalar(0))
    snowstore = snowstore + snow - actmelt + refreeze
    return snowstore, refreeze


# -Function that determines the maximum amount of water that can be stored in the snowpack
def MaxSnowWatStorage(snowsc, snowstore):
    maxsnowwatstore = snowsc * snowstore
    return maxsnowwatstore


# -Function to calculate the actual snow water storage (binary thaw based on Tavg)
def SnowWatStorage(
    pcr, temp, tcrit, maxsnowwatstore, snowwatstore, actmelt, rain, refreeze
):
    # If Tavg < Tcrit: all liquid water was refrozen, SWS = 0
    # If Tavg >= Tcrit: remove refrozen water (=0), add melt and rain, cap at max
    snowwatstore = pcr.ifthenelse(
        temp < tcrit,
        pcr.scalar(0),
        pcr.min(maxsnowwatstore, pcr.max(0, snowwatstore - refreeze + actmelt + rain)),
    )
    return snowwatstore


# -Function to calculate the total snow storage (snowstore + snowwatstore)
def TotSnowStorage(snowstore, snowwatstore, snowfrac, rainfrac):
    totalsnowstore = (snowstore + snowwatstore) * (snowfrac + rainfrac)
    return totalsnowstore


# -Function to calculate runoff from snow
def SnowR(pcr, snowwatstore, maxsnowwatstore, actmelt, rain, oldsnowwatstore, snowfrac):
    snowr = pcr.ifthenelse(
        snowwatstore == maxsnowwatstore,
        (((actmelt + rain) - (snowwatstore - oldsnowwatstore)) * snowfrac),
        0,
    )
    return snowr


# -init snow processes
def init(self, pcr, config):
    pars = ["Tcrit", "SnowSC", "DDFS", "SnowF", "SnowCth"]
    for i in pars:
        try:
            setattr(self, i, pcr.readmap(self.inpath + config.get("SNOW", i)))
        except:
            setattr(self, i, config.getfloat("SNOW", i))


# -initial snow processes
def initial(self, pcr, config):
    try:
        self.SnowStore = config.getfloat("SNOW_INIT", "SnowIni")
    except:
        self.SnowStore = pcr.readmap(self.inpath + config.get("SNOW_INIT", "SnowIni"))
    # -initial water stored in snowpack
    try:
        self.SnowWatStore = config.getfloat("SNOW_INIT", "SnowWatStore")
    except:
        self.SnowWatStore = pcr.readmap(
            self.inpath + config.get("SNOW_INIT", "SnowWatStore")
        )
    self.TotalSnowStore = self.SnowStore + self.SnowWatStore


# -dynamic snow processes
def dynamic(
    self,
    pcr,
    Temp,
    TempMin,
    TempMax,
    Precip,
    Snow_GLAC,
    ActSnowMelt_GLAC,
    SnowFrac,
    RainFrac,
    SnowR_GLAC,
):
    # -Snow and rain differentiation
    Snow = pcr.ifthenelse(Temp >= self.Tcrit, 0, Precip)
    Rain = pcr.ifthenelse(Temp < self.Tcrit, 0, Precip)
    # -Report Snow for entire cell (snow+glacier fraction)
    self.reporting.reporting(self, pcr, "TotSnow", Snow)
    self.reporting.reporting(
        self, pcr, "TotSnowF", Snow * (1 - self.GlacFrac) + Snow_GLAC
    )
    # -Snow melt (uses positive degree-hours from sinusoidal of Tmin/Tmax)
    PotSnowMelt = pcr.ifthenelse(
        TempMax < 0, 0, self.snow.PotSnowMelt(pcr, TempMin, TempMax, self.DDFS)
    )
    ActSnowMelt = self.snow.ActSnowMelt(pcr, self.SnowStore, PotSnowMelt)
    # -Report snow melt for entire cell (snow+glacier fraction)
    self.reporting.reporting(self, pcr, "TotSnowMelt", ActSnowMelt)
    self.reporting.reporting(
        self, pcr, "TotSnowMeltF", ActSnowMelt * (1 - self.GlacFrac) + ActSnowMelt_GLAC
    )
    # -Update snow store (binary refreeze: if Tavg < Tcrit, all liquid water refreezes)
    self.SnowStore, Refreeze = self.snow.SnowStoreUpdate(
        pcr, self.SnowStore, Snow, ActSnowMelt, Temp, self.Tcrit, self.SnowWatStore
    )
    # -Caclulate the maximum amount of water that can be stored in snowwatstore
    MaxSnowWatStore = self.snow.MaxSnowWatStorage(self.SnowSC, self.SnowStore)
    OldSnowWatStore = self.SnowWatStore
    # -Calculate the actual amount of water stored in snowwatstore (binary: Tavg < Tcrit → 0)
    self.SnowWatStore = self.snow.SnowWatStorage(
        pcr,
        Temp,
        self.Tcrit,
        MaxSnowWatStore,
        self.SnowWatStore,
        ActSnowMelt,
        Rain,
        Refreeze,
    )
    # -Changes in total water storage in snow (SnowStore and SnowWatStore)
    OldTotalSnowStore = self.TotalSnowStore
    self.TotalSnowStore = (
        self.snow.TotSnowStorage(self.SnowStore, self.SnowWatStore, SnowFrac, RainFrac)
        + self.TotalSnowStore_GLAC
    )  # for entire cell
    # -Report Snow storage
    self.reporting.reporting(self, pcr, "StorSnow", self.TotalSnowStore)
    # -Determine if cell is covered with snow
    SnowCover = pcr.ifthenelse(
        self.TotalSnowStore > self.SnowCth, pcr.scalar(1), pcr.scalar(0)
    )
    self.reporting.reporting(self, pcr, "SCover", SnowCover)
    self.reporting.reporting(
        self, pcr, "StorSnowW", self.SnowWatStore
    )  # sonu added note this is only SnowWatStore
    # -Snow runoff
    SnowR = (
        self.snow.SnowR(
            pcr,
            self.SnowWatStore,
            MaxSnowWatStore,
            ActSnowMelt,
            Rain,
            OldSnowWatStore,
            SnowFrac,
        )
        + SnowR_GLAC
    )  # for entire cell
    ##sonu added snow infiltration##
    SnowSoil = SnowR * self.SnowF
    SnowR = SnowR * (1 - self.SnowF)
    SnowR = SnowR * (1 - self.openWaterFrac)
    # -Report Snow runoff
    self.reporting.reporting(self, pcr, "TotSnowRF", SnowR)

    return Rain, SnowR, SnowSoil, OldTotalSnowStore
