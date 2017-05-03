import math
import numpy as np
import wxgen.util

"""
Module related to radiation calculations
"""

def swing(jday, hour, lat, lon, cloud_cover=0, pressure=1013, temperature=15):
   """
   Solar global radiation (Short-Wave INcominG solar radiation)

   All parameters can be provided in vector form for faster computation.
   Temperature and pressure is used in refraction calculation.

   Arguments:
      jday (np.array): Julian day
      hour (np.array): Hour of the day in UTC
      lat (np.array): Latitude in degrees
      lon (np.array): Longitude in degrees
      cloud_cover (np.array): cloud_cover (0-1)
      pressure (np.array): pressure (hPa)
      temperature (np.array): temperature (C)

   Returns:
      np.array: radiation in the same shape as the input fields (W/m2)

   C.Lussana - 02.05.2017 (cristianl@met.no)
   Credits: 
      routines used to compute solar position are taken from the shareware
      Linux GIS grass v5.7.0. In particular see the source code solpos00.c 
      and the directory sunmask
      7 September 2004 Cristian Lussana
   """
   shape = cloud_cover.shape

   year = 2017 * np.ones(shape)

   utime = hour
   delta = year - 1949.
   leap  = (delta / 4.).astype(int)
   daynum = jday
   julday = 32916.5 + delta * 365. + leap + daynum + utime / 24.

   # constants
   raddeg = math.pi / 180
   degrad = 180 / math.pi
   solcon = 1367.
   # compute the julian day, if needed
   # Time used in the calculation of ecliptic coordinates
   # Noon 1 JAN 2000 = 2400000 + 51545 days Julian data
   #   Michalsky, J. 1988. The Astronomical Almanac's algorithm for approximate
   #   solar position (1950-2050). Solar Energy 40 (3), pp.227-235
   ectime = julday - 51545.

   # Mean longitude
   #   Michalsky, J. 1988. The Astronomical Almanac's algorithm for approximate
   #   solar position (1950-2050). Solar Energy 40 (3), pp.227-235
   mnlong = 280.460 + 0.9856474 * ectime
   # (dump the multiplies of 360, so the answer is between 0 and 360)
   mnlong = mnlong - 360. * (mnlong / 360.).astype(int)
   mnlong[mnlong < 0] += 360.
   # Mean anomaly
   #   Michalsky, J. 1988. The Astronomical Almanac's algorithm for approximate
   #   solar position (1950-2050). Solar Energy 40 (3), pp.227-235
   mnanom = 357.528 + 0.9856003 * ectime
   # (dump the multiplies of 360, so the answer is between 0 and 360)
   mnanom = mnanom - 360. * (mnanom / 360.).astype(int)
   mnanom[mnanom < 0] += 360.
   # Ecliptic longitude
   #   Michalsky, J. 1988. The Astronomical Almanac's algorithm for approximate
   #   solar position (1950-2050). Solar Energy 40 (3), pp.227-235
   eclong = mnlong + 1.915 * np.sin(mnanom * raddeg) + 0.020 * np.sin(2. * mnanom * raddeg)
   # (dump the multiplies of 360, so the answer is between 0 and 360)
   eclong = eclong - 360. * (eclong / 360.).astype(int)
   eclong[eclong < 0] += 360.

   """ Obliquity of the ecliptic
      Michalsky, J. 1988. The Astronomical Almanac's algorithm for approximate
      solar position (1950-2050). Solar Energy 40 (3), pp.227-235
   """
   ecobli = 23.439 - 4.0e-07 * ectime

   """ Declination
      Michalsky, J. 1988. The Astronomical Almanac's algorithm for approximate
      solar position (1950-2050). Solar Energy 40 (3), pp.227-235
      (No adjustement for century non-leap years since this function is bounded by 1950 - 2050)
   """
   declin = degrad * np.arcsin( np.sin(ecobli * raddeg) * np.sin(eclong * raddeg) )

   """ Right ascension
      Michalsky, J. 1988. The Astronomical Almanac's algorithm for approximate
      solar position (1950-2050). Solar Energy 40 (3), pp.227-235
   """
   top = np.cos(raddeg * ecobli) * np.sin(raddeg * eclong)
   bottom = np.cos(raddeg * eclong)
   rascen = degrad * np.arctan2(top, bottom)
   # (make it a positive angle)
   rascen[rascen < 0] += 360.

   """ Greenwich mean sideral time
      Michalsky, J. 1988. The Astronomical Almanac's algorithm for approximate
      solar position (1950-2050). Solar Energy 40 (3), pp.227-235
   """
   gmst = 6.697375 + 0.0657098242 * ectime + utime
   # (dump the multiplies of 24, so the answer is between 0 and 24)
   gmst = gmst - 24. * (gmst / 24.).astype(int)
   gmst[gmst < 0] += 24.

   """ Local mean sideral time
      Michalsky, J. 1988. The Astronomical Almanac's algorithm for approximate
      solar position (1950-2050). Solar Energy 40 (3), pp.227-235
   """
   lmst = gmst * 15. + lon
   # (dump the multiplies of 360, so the answer is between 0 and 360)
   lmst = lmst - 360. * (lmst / 360.).astype(int)
   lmst[lmst < 0] += 360.

   """ Hour angle
      Michalsky, J. 1988. The Astronomical Almanac's algorithm for approximate
      solar position (1950-2050). Solar Energy 40 (3), pp.227-235
   """
   hrang = lmst - rascen
   # (force it between -180 and 180 degrees)
   hrang[hrang < -180] += 360.
   hrang[hrang > 180] -= 360.

   """ ETR solar zenith angle
   Reference:
      Iqbal, M. 1983. An introduction to solar radiation
      Academic Press, NY., page 3
   a bit of trigonometry
   """
   cd = np.cos(raddeg * declin)
   ch = np.cos(raddeg * hrang)
   cl = np.cos(raddeg * lat)
   sd = np.sin(raddeg * declin)
   sl = np.sin(raddeg * lat)

   cz = sd * sl + cd * cl * ch
   # watch out for the roundoff errors
   cz[cz > 1] = 1.
   cz[cz < -1] = -1.
   zenetr = np.arccos(cz) * degrad
   # limit the degrees below the horizon to 9 [+90 -> 99]
   zenetr[zenetr > 99] = 99.
   elevetr = 90. - zenetr

   """ Refraction correction, degrees
   Reference:
      Zimmerman, john C. 1981. Sun-pointing programs and their accuracy.
      SAND81-0761, Experimental System Operation Division 4271,
      Sandia National Laboratories, Albuquerque, NM.
   prestemp    # temporary pressure/temperature correction
   refcor      # temporary refraction correction
   tanelev     # tangent of the solar elevation angle
   """
   # If the sun is near zenith, the algorithm bombs; refraction near 0
   refcor = np.zeros(shape)
   # Otherwise, we have refraction
   tanelev = np.tan(raddeg * elevetr)
   I = np.where((elevetr <= 85) & (elevetr >= 5))[0]
   refcor[I] = 58.1 / tanelev[I] - 0.07 / tanelev[I]**3 + 0.000086 / tanelev[I]**5
   I = np.where((elevetr < 5) & (elevetr >= -0.575))[0]
   refcor[I] = 1735. + elevetr[I] * ( -518.2 + elevetr[I] * (103.4 + elevetr[I] * ( -12.79 + elevetr[I] * 0.711)) )
   I = np.where(elevetr < -0.575)[0]
   refcor[I] = -20.774 / tanelev[I]
   I = np.where(elevetr <= 85)[0]
   prestemp = (pressure[I] * 283.) / (1013. * (273. + temperature[I]))
   refcor[I] = refcor[I] * prestemp / 3600.

   elevref = elevetr + refcor
   # limit the degrees below the horizon to 9
   elevref[elevref < -9] = -9
   # refracted solar zenith angle
   zenref = 90. - elevref
   coszen = np.cos(raddeg * zenref)

   """ Global Radiation
   Reference:
      Holtslag, A. A. M., and A. P. Van Ulden. "A simple scheme for daytime 
      estimates of the surface fluxes from routine weather data." 
      Journal of climate and Applied Meteorology 22.4 (1983): 517-529.
   NOTE: Tab2: a1,a2 used here are the one estimated for De Bilt
   modified as described in:
   Sozzi, Roberto, Mauro Valentini, and Teodoro Georgiadis. Introduzione 
    alla turbolenza atmosferica: concetti, stime, misure. Pitagora, 2002.
   """
   sinpsi = np.sin(raddeg * elevref)
   a1 = 1041.
   a2 = -69.
   b1 = -0.75
   b2 = 3.4
   sinmin = -a2 / a1

   sunrad = np.zeros(shape)
   I = np.where(sinpsi >= sinmin)[0]
   sunrad[I] = (a1 * sinpsi[I] * np.exp(-0.057 / sinpsi[I])) *(1. + b1 * cloud_cover[I] ** b2)

   return sunrad
