/*------------------------------------------------------------------------------
* tides.c : tidal displacement corrections
*
*          Copyright (C) 2015-2017 by T.TAKASU, All rights reserved.
*
* references :
*     [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*     [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*         2003, November 2003
*     [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*         Space Technology Library, 2004
*     [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*         May 2009
*     [5] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*
* version : $Revision:$ $Date:$
* history : 2015/05/10 1.0  separated from ppp.c
*           2015/06/11 1.1  fix bug on computing days in tide_oload() (#128)
*           2017/04/11 1.2  fix bug on calling geterp() in timdedisp()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define SQR(x)      ((x)*(x))

// The following few functions are in support of dehanttideinel. This code is
// a translation of the respective iers fortran code to C and to RTKLIB, and
// the reference code is included in the lib/iers/src/ directory.

// This subroutine gives the in-phase and out-of-phase corrections induced by
// mantle anelasticity in the long period band.
static void step2lon_(const double xsta[3], double t, double xcorsta[3]) {
  static const double datdi[5][9] = {{0., 0., 0., 1., 0., .47, .23, .16, .07},
                                     {0., 2., 0., 0., 0., -.2, -.12, -.11, -.05},
                                     {1., 0., -1., 0., 0., -.11, -.08, -.09, -.04},
                                     {2., 0., 0., 0., 0., -.13, -.11, -.15, -.07},
                                     {2., 0., 0., 1., 0., -.05, -.05, -.06, -.03}};

  //  Compute the phase angles in degrees.
  double s = ((t * 1.85139e-6 - .0014663889) * t + 481267.88194) * t + 218.31664563;
  double pr = (((t * 7e-9 + 2.1e-8) * t + 3.08889e-4) * t + 1.396971278) * t;
  s += pr;
  double h = (((t * -6.54e-9 + 2e-8) * t + 3.0322222e-4) * t + 36000.7697489) * t + 280.46645;
  double p =
      (((t * 5.263e-8 - 1.24991e-5) * t - .01032172222) * t + 4069.01363525) * t + 83.35324312;
  double zns =
      (((t * 1.65e-8 - 2.13944e-6) * t - .00207561111) * t + 1934.13626197) * t + 234.95544499;
  double ps =
      (((t * -3.34e-9 - 1.778e-8) * t + 4.5688889e-4) * t + 1.71945766667) * t + 282.93734098;
  double rsta = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1] + xsta[2] * xsta[2]);
  double sinphi = xsta[2] / rsta;
  double cosphi = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1]) / rsta;
  double cosla = xsta[0] / cosphi / rsta;
  double sinla = xsta[1] / cosphi / rsta;
  // Reduce angles to between the range 0 and 360.
  s = fmod(s, 360);
  //  tau = fmod(tau, 360)
  h = fmod(h, 360);
  p = fmod(p, 360);
  zns = fmod(zns, 360);
  ps = fmod(ps, 360);
  double dr_tot = 0, dn_tot = 0.;
  for (int i = 0; i < 3; ++i) xcorsta[i] = 0;
  for (int j = 0; j < 5; ++j) {
    double thetaf = (datdi[j][0] * s + datdi[j][1] * h + datdi[j][2] * p + datdi[j][3] * zns +
                     datdi[j][4] * ps) *
                    D2R;
    double dr = datdi[j][5] * (sinphi * sinphi * 3 - 1) / 2 * cos(thetaf) +
                datdi[j][7] * (sinphi * sinphi * 3 - 1) / 2 * sin(thetaf);
    double dn = datdi[j][6] * (cosphi * sinphi * 2) * cos(thetaf) +
                datdi[j][8] * (cosphi * sinphi * 2) * sin(thetaf);
    double de = 0;
    dr_tot += dr;
    dn_tot += dn;
    xcorsta[0] += dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
    xcorsta[1] += dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
    xcorsta[2] += dr * sinphi + dn * cosphi;
  }
  for (int i = 0; i < 3; ++i) xcorsta[i] /= 1e3;
}

// This subroutine gives the in-phase and out-of-phase corrections induced by
// mantle anelasticity in the diurnal band.
static void step2diu_(const double xsta[3], double fhr, double t, double xcorsta[3]) {
  static const double datdi[31][9] = {{-3, 0, 2, 0, 0, -0.01, 0, 0, 0},
                                      {-3, 2, 0, 0, 0, -0.01, 0, 0, 0},
                                      {-2, 0, 1, -1, 0, -0.02, 0, 0, 0},
                                      {-2, 0, 1, 0, 0, -0.08, 0, -0.01, 0.01},
                                      {-2, 2, -1, 0, 0, -0.02, 0, 0, 0},
                                      {-1, 0, 0, -1, 0, -0.10, 0, 0, 0},
                                      {-1, 0, 0, 0, 0, -0.51, 0, -0.02, 0.03},
                                      {-1, 2, 0, 0, 0, 0.01, 0, 0, 0},
                                      {0, -2, 1, 0, 0, 0.01, 0, 0, 0},
                                      {0, 0, -1, 0, 0, 0.02, 0, 0, 0},
                                      {0, 0, 1, 0, 0, 0.06, 0, 0, 0},
                                      {0, 0, 1, 1, 0, 0.01, 0, 0, 0},
                                      {0, 2, -1, 0, 0, 0.01, 0, 0, 0},
                                      {1, -3, 0, 0, 1, -0.06, 0, 0, 0},
                                      {1, -2, 0, -1, 0, 0.01, 0, 0, 0},
                                      {1, -2, 0, 0, 0, -1.23, -0.07, 0.06, 0.01},
                                      {1, -1, 0, 0, -1, 0.02, 0, 0, 0},
                                      {1, -1, 0, 0, 1, 0.04, 0, 0, 0},
                                      {1, 0, 0, -1, 0, -0.22, 0.01, 0.01, 0},
                                      {1, 0, 0, 0, 0, 12.00, -0.80, -0.67, -0.03},
                                      {1, 0, 0, 1, 0, 1.73, -0.12, -0.10, 0},
                                      {1, 0, 0, 2, 0, -0.04, 0, 0, 0},
                                      {1, 1, 0, 0, -1, -0.50, -0.01, 0.03, 0},
                                      {1, 1, 0, 0, 1, 0.01, 0, 0, 0},
                                      {0, 1, 0, 1, -1, -0.01, 0, 0, 0},
                                      {1, 2, -2, 0, 0, -0.01, 0, 0, 0},
                                      {1, 2, 0, 0, 0, -0.11, 0.01, 0.01, 0},
                                      {2, -2, 1, 0, 0, -0.01, 0, 0, 0},
                                      {2, 0, -1, 0, 0, -0.02, 0, 0, 0},
                                      {3, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {3, 0, 0, 1, 0, 0, 0, 0, 0}};
  //  Compute the phase angles in degrees.
  double s = ((t * 1.85139e-6 - .0014663889) * t + 481267.88194) * t + 218.31664563;
  double tau =
      fhr * 15. + 280.4606184 + ((t * -2.58e-8 + 3.8793e-4) * t + 36000.7700536) * t + (-s);
  double pr = (((t * 7e-9 + 2.1e-8) * t + 3.08889e-4) * t + 1.396971278) * t;
  s += pr;
  double h__ = (((t * -6.54e-9 + 2e-8) * t + 3.0322222e-4) * t + 36000.7697489) * t + 280.46645;
  double p =
      (((t * 5.263e-8 - 1.24991e-5) * t - 0.01032172222) * t + 4069.01363525) * t + 83.35324312;
  double zns =
      (((t * 1.65e-8 - 2.13944e-6) * t - 0.00207561111) * t + 1934.13626197) * t + 234.95544499;
  double ps =
      (((t * -3.34e-9 - 1.778e-8) * t + 4.5688889e-4) * t + 1.71945766667) * t + 282.93734098;
  // Reduce angles to between the range 0 and 360.
  s = fmod(s, 360);
  tau = fmod(tau, 360);
  h__ = fmod(h__, 360);
  p = fmod(p, 360);
  zns = fmod(zns, 360);
  ps = fmod(ps, 360);
  double rsta = sqrt(SQR(xsta[0]) + SQR(xsta[1]) + SQR(xsta[2]));
  double sinphi = xsta[2] / rsta;
  double cosphi = sqrt(SQR(xsta[0]) + SQR(xsta[1])) / rsta;
  double cosla = xsta[0] / cosphi / rsta;
  double sinla = xsta[1] / cosphi / rsta;
  double zla = atan2(xsta[1], xsta[0]);
  // Initialize.
  for (int i = 0; i < 3; ++i) xcorsta[i] = 0.;
  for (int j = 0; j < 31; ++j) {
    // Convert from degrees to radians.
    double thetaf = (tau + datdi[j][0] * s + datdi[j][1] * h__ + datdi[j][2] * p +
                     datdi[j][3] * zns + datdi[j][4] * ps) *
                    D2R;
    double dr = datdi[j][5] * 2 * sinphi * cosphi * sin(thetaf + zla) +
                datdi[j][6] * 2 * sinphi * cosphi * cos(thetaf + zla);
    double dn = datdi[j][7] * (SQR(cosphi) - SQR(sinphi)) * sin(thetaf + zla) +
                datdi[j][8] * (SQR(cosphi) - SQR(sinphi)) * cos(thetaf + zla);
    double de = datdi[j][7] * sinphi * cos(thetaf + zla) - datdi[j][8] * sinphi * sin(thetaf + zla);
    xcorsta[0] += dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
    xcorsta[1] += dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
    xcorsta[2] += dr * sinphi + dn * cosphi;
  }
  for (int i = 0; i < 3; ++i) xcorsta[i] /= 1e3;
}

// This subroutine gives the out-of-phase corrections induced by mantle
// anelasticity in the diurnal band.
static void st1idiu_(const double xsta[3], const double xsun[3], const double xmon[3],
                     double fac2sun, double fac2mon, double xcorsta[3]) {
  const double dhi = -0.0025, dli = -7e-4;
  // Compute the normalized position vector of the IGS station.
  double rsta = norm(xsta, 3);
  double sinphi = xsta[2] / rsta;
  double cosphi = sqrt(SQR(xsta[0]) + SQR(xsta[1])) / rsta;
  double cos2phi = SQR(cosphi) - SQR(sinphi);
  double sinla = xsta[1] / cosphi / rsta;
  double cosla = xsta[0] / cosphi / rsta;
  // Compute the normalized position vector of the Moon.
  double rmon = norm(xmon, 3);
  // Compute the normalized position vector of the Sun.
  double rsun = norm(xsun, 3);
  double drsun = dhi * -3. * sinphi * cosphi * fac2sun * xsun[2] *
                 (xsun[0] * sinla - xsun[1] * cosla) / (rsun * rsun);
  double drmon = dhi * -3. * sinphi * cosphi * fac2mon * xmon[2] *
                 (xmon[0] * sinla - xmon[1] * cosla) / (rmon * rmon);
  double dnsun =
      dli * -3. * cos2phi * fac2sun * xsun[2] * (xsun[0] * sinla - xsun[1] * cosla) / (rsun * rsun);
  double dnmon =
      dli * -3. * cos2phi * fac2mon * xmon[2] * (xmon[0] * sinla - xmon[1] * cosla) / (rmon * rmon);
  double desun =
      dli * -3. * sinphi * fac2sun * xsun[2] * (xsun[0] * cosla + xsun[1] * sinla) / (rsun * rsun);
  double demon =
      dli * -3. * sinphi * fac2mon * xmon[2] * (xmon[0] * cosla + xmon[1] * sinla) / (rmon * rmon);
  double dr = drsun + drmon;
  double dn = dnsun + dnmon;
  double de = desun + demon;
  //  Compute the corrections for the station.
  xcorsta[0] = dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
  xcorsta[1] = dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
  xcorsta[2] = dr * sinphi + dn * cosphi;
}

// This subroutine gives the out-of-phase corrections induced by mantle
// anelasticity in the semi-diurnal band.
static void st1isem_(const double xsta[3], const double xsun[3], const double xmon[3],
                     double fac2sun, double fac2mon, double xcorsta[3]) {
  const double dhi = -0.0022, dli = -7e-4;
  // Compute the normalized position vector of the IGS station.
  double rsta = norm(xsta, 3);
  double sinphi = xsta[2] / rsta;
  double cosphi = sqrt(SQR(xsta[0]) + SQR(xsta[1])) / rsta;
  double sinla = xsta[1] / cosphi / rsta;
  double cosla = xsta[0] / cosphi / rsta;
  double costwola = SQR(cosla) - SQR(sinla);
  double sintwola = cosla * 2. * sinla;
  // Compute the normalized position vector of the Moon.
  double rmon = norm(xmon, 3);
  // Compute the normalized position vector of the Sun.
  double rsun = norm(xsun, 3);
  double drsun = -3.0 / 4.0 * dhi * SQR(cosphi) * fac2sun *
                 ((SQR(xsun[0]) - SQR(xsun[1])) * sintwola - xsun[0] * 2. * xsun[1] * costwola) /
                 SQR(rsun);
  double drmon = -3.0 / 4.0 * dhi * SQR(cosphi) * fac2mon *
                 ((SQR(xmon[0]) - SQR(xmon[1])) * sintwola - xmon[0] * 2. * xmon[1] * costwola) /
                 SQR(rmon);
  double dnsun = 3.0 / 2.0 * dli * sinphi * cosphi * fac2sun *
                 ((SQR(xsun[0]) - SQR(xsun[1])) * sintwola - xsun[0] * 2. * xsun[1] * costwola) /
                 SQR(rsun);
  double dnmon = 3.0 / 2.0 * dli * sinphi * cosphi * fac2mon *
                 ((SQR(xmon[0]) - SQR(xmon[1])) * sintwola - xmon[0] * 2. * xmon[1] * costwola) /
                 SQR(rmon);
  double desun = -3.0 / 2.0 * dli * cosphi * fac2sun *
                 ((SQR(xsun[0]) - SQR(xsun[1])) * costwola + xsun[0] * 2. * xsun[1] * sintwola) /
                 SQR(rsun);
  double demon = -3.0 / 2.0 * dli * cosphi * fac2mon *
                 ((SQR(xmon[0]) - SQR(xmon[1])) * costwola + xmon[0] * 2. * xmon[1] * sintwola) /
                 SQR(rmon);
  double dr = drsun + drmon;
  double dn = dnsun + dnmon;
  double de = desun + demon;
  xcorsta[0] = dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
  xcorsta[1] = dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
  xcorsta[2] = dr * sinphi + dn * cosphi;
}

// This subroutine gives the corrections induced by the latitude dependence
// given by L^1 in Mathews et al. 1991 (See References).
static void st1l1_(const double xsta[3], const double xsun[3], const double xmon[3], double fac2sun,
                   double fac2mon, double xcorsta[3]) {
  const double l1d = 0.0012, l1sd = 0.0024;
  // Compute the normalized position vector of the IGS station.
  double rsta = norm(xsta, 3);
  double sinphi = xsta[2] / rsta;
  double cosphi = sqrt(SQR(xsta[0]) + SQR(xsta[1])) / rsta;
  double sinla = xsta[1] / cosphi / rsta;
  double cosla = xsta[0] / cosphi / rsta;
  // Compute the normalized position vector of the Moon.
  double rmon = norm(xmon, 3);
  // Compute the normalized position vector of the Sun.
  double rsun = norm(xsun, 3);
  // Compute the station corrections for the diurnal band.
  double l1 = l1d;
  double dnsun =
      -l1 * SQR(sinphi) * fac2sun * xsun[2] * (xsun[0] * cosla + xsun[1] * sinla) / SQR(rsun);
  double dnmon =
      -l1 * SQR(sinphi) * fac2mon * xmon[2] * (xmon[0] * cosla + xmon[1] * sinla) / SQR(rmon);
  double desun = l1 * sinphi * (SQR(cosphi) - SQR(sinphi)) * fac2sun * xsun[2] *
                 (xsun[0] * sinla - xsun[1] * cosla) / SQR(rsun);
  double demon = l1 * sinphi * (SQR(cosphi) - SQR(sinphi)) * fac2mon * xmon[2] *
                 (xmon[0] * sinla - xmon[1] * cosla) / SQR(rmon);
  double de = 3 * (desun + demon);
  double dn = 3 * (dnsun + dnmon);
  xcorsta[0] = -de * sinla - dn * sinphi * cosla;
  xcorsta[1] = de * cosla - dn * sinphi * sinla;
  xcorsta[2] = dn * cosphi;
  // Compute the station corrections for the semi-diurnal band.
  l1 = l1sd;
  double costwola = SQR(cosla) - SQR(sinla);
  double sintwola = 2 * cosla * sinla;
  dnsun = -l1 / 2 * sinphi * cosphi * fac2sun *
          ((SQR(xsun[0]) - SQR(xsun[1])) * costwola + xsun[0] * 2 * xsun[1] * sintwola) / SQR(rsun);
  dnmon = -l1 / 2 * sinphi * cosphi * fac2mon *
          ((SQR(xmon[0]) - SQR(xmon[1])) * costwola + xmon[0] * 2 * xmon[1] * sintwola) / SQR(rmon);
  desun = -l1 / 2 * SQR(sinphi) * cosphi * fac2sun *
          ((SQR(xsun[0]) - SQR(xsun[1])) * sintwola - xsun[0] * 2 * xsun[1] * costwola) / SQR(rsun);
  demon = -l1 / 2 * SQR(sinphi) * cosphi * fac2mon *
          ((SQR(xmon[0]) - SQR(xmon[1])) * sintwola - xmon[0] * 2 * xmon[1] * costwola) / SQR(rmon);
  de = 3 * (desun + demon);
  dn = 3 * (dnsun + dnmon);
  xcorsta[0] += -de * sinla - dn * sinphi * cosla;
  xcorsta[1] += +de * cosla - dn * sinphi * sinla;
  xcorsta[2] += dn * cosphi;
}

// This subroutine computes the station tidal displacement caused by lunar and
// solar gravitational attraction (see References).
static void dehanttideinel(gtime_t tutc, const double xsta[3], const double xsun[3],
                           const double xmon[3], double dxtide[3]) {
  // Nominal second degree and third degree love numbers and shida numbers.
  const double h20 = 0.6078, l20 = 0.0847, h3 = 0.292, l3 = 0.015;
  // Scalar product of station vector with sun/moon vector.
  double rsta = sqrt(SQR(xsta[0]) + SQR(xsta[1]) + SQR(xsta[2]));
  double rsun = sqrt(SQR(xsun[0]) + SQR(xsun[1]) + SQR(xsun[2]));
  double rmon = sqrt(SQR(xmon[0]) + SQR(xmon[1]) + SQR(xmon[2]));
  double scs = dot(xsta, xsun, 3);
  double scm = dot(xsta, xmon, 3);
  double scsun = scs / rsta / rsun;
  double scmon = scm / rsta / rmon;
  // Computation of new h2 and l2.
  double cosphi = sqrt(SQR(xsta[0]) + SQR(xsta[1])) / rsta;
  double h2 = h20 - (1 - 3.0 / 2.0 * SQR(cosphi)) * 6e-4;
  double l2 = l20 + (1 - 3.0 / 2.0 * SQR(cosphi)) * 2e-4;
  // P2 term
  double p2sun = 3 * (h2 / 2 - l2) * SQR(scsun) - h2 / 2;
  double p2mon = 3 * (h2 / 2 - l2) * SQR(scmon) - h2 / 2;
  // P3 term
  double p3sun = 5.0 / 2.0 * (h3 - l3 * 3) * (SQR(scsun) * scsun) + 3.0 / 2.0 * (l3 - h3) * scsun;
  double p3mon = 5.0 / 2.0 * (h3 - l3 * 3) * (SQR(scmon) * scmon) + 3.0 / 2.0 * (l3 - h3) * scmon;
  // Term in direction of sun/moon vector.
  double x2sun = 3 * l2 * scsun;
  double x2mon = 3 * l2 * scmon;
  double x3sun = 3.0 / 2.0 * l3 * (SQR(scsun) * 5 - 1);
  double x3mon = 3.0 / 2.0 * l3 * (SQR(scmon) * 5 - 1);
  // Factors for sun/moon using iau current best estimates (see references).
  const double mass_ratio_sun = 332946.0482;
  const double mass_ratio_moon = 0.0123000371;
  const double re = 6378136.6;
  double resun = re / rsun;
  double fac2sun = mass_ratio_sun * re * SQR(resun) * resun;
  double remon = re / rmon;
  double fac2mon = mass_ratio_moon * re * SQR(remon) * remon;
  double fac3sun = fac2sun * resun;
  double fac3mon = fac2mon * remon;
  // Total displacement.
  for (int i = 0; i < 3; i++) {
    dxtide[i] = fac2sun * (x2sun * xsun[i] / rsun + p2sun * xsta[i] / rsta) +
                fac2mon * (x2mon * xmon[i] / rmon + p2mon * xsta[i] / rsta) +
                fac3sun * (x3sun * xsun[i] / rsun + p3sun * xsta[i] / rsta) +
                fac3mon * (x3mon * xmon[i] / rmon + p3mon * xsta[i] / rsta);
  }
  // Corrections for the out-of-phase part of love numbers (part h_2^(0)i
  // and l_2^(0)i )
  //
  // First, for the diurnal band.
  double xcorsta[3] = {0};
  st1idiu_(xsta, xsun, xmon, fac2sun, fac2mon, xcorsta);
  for (int i = 0; i < 3; ++i) dxtide[i] += xcorsta[i];
  // Second, for the semi-diurnal band.
  st1isem_(xsta, xsun, xmon, fac2sun, fac2mon, xcorsta);
  for (int i = 0; i < 3; ++i) dxtide[i] += xcorsta[i];
  // Corrections for the latitude dependence of love numbers (part l^(1) )
  st1l1_(xsta, xsun, xmon, fac2sun, fac2mon, xcorsta);
  for (int i = 0; i < 3; ++i) dxtide[i] += xcorsta[i];
  // Consider corrections for step 2.
  // Corrections for the diurnal band:
  //  First, we need to know the date converted in julian centuries
  double ep[6];
  time2epoch(tutc, ep);
  double fhr = ep[3] + ep[4] / 60.0 + ep[5] / 3600.0;

  // Terrestrial time.
  gtime_t tgps = utc2gpst(tutc);
  const double ep2000[] = {2000, 1, 1, 11, 59, 08.816};  // GPST of J2000.0
  double t = timediff(tgps, epoch2time(ep2000)) / 86400.0 / 36525.0;

  // Second, we can call the subroutine step2diu, for the diurnal band
  // corrections, (in-phase and out-of-phase frequency dependence):
  step2diu_(xsta, fhr, t, xcorsta);
  for (int i = 0; i < 3; ++i) dxtide[i] += xcorsta[i];
  // Corrections for the long-period band,
  // (in-phase and out-of-phase frequency dependence):
  step2lon_(xsta, t, xcorsta);
  for (int i = 0; i < 3; ++i) dxtide[i] += xcorsta[i];

#ifdef RTK_DISABLED
  // Consider corrections for step 3.
  // Uncorrect for the permanent tide.
  double sinphi = xsta[1] / rsta;
  cosphi = sqrt(SQR(xsta[0]) + SQR(xsta[1])) / rsta;
  double cosla = xsta[0] / cosphi / rsta;
  double sinla = xsta[1] / cosphi / rsta;
  double dr = -sqrt(5.0 / 4.0 / PI) * h2 * 0.3146 * (SQR(sinphi) * 3.0 / 2.0 - 0.5);
  double dn = -sqrt(5.0 / 4.0 / PI) * l2 * 0.3146 * 3 * cosphi * sinphi;
  dxtide[0] = dxtide[0] - dr * cosla * cosphi + dn * cosla * sinphi;
  dxtide[1] = dxtide[1] - dr * sinla * cosphi + dn * sinla * sinphi;
  dxtide[2] = dxtide[2] - dr * sinphi - dn * cosphi;
#endif
}

/* displacement by ocean tide loading (ref [2] 7) ----------------------------*/
static void tide_oload(gtime_t tut, const double *odisp, double *denu)
{
    const double args[][5]={
        {1.40519E-4, 2.0,-2.0, 0.0, 0.00},  // M2 semidiurnal principal lunar.
        {1.45444E-4, 0.0, 0.0, 0.0, 0.00},  // S2 semidiurnal principal solar.
        {1.37880E-4, 2.0,-3.0, 1.0, 0.00},  // N2 semidiurnal lunar elliptical.
        {1.45842E-4, 2.0, 0.0, 0.0, 0.00},  // K2 semidiurnal luni-solar declination.
        {0.72921E-4, 1.0, 0.0, 0.0, 0.25},  // K1 diurnal luni-solar declination.
        {0.67598E-4, 1.0,-2.0, 0.0,-0.25},  // O1 diurnal principal lunar.
        {0.72523E-4,-1.0, 0.0, 0.0,-0.25},  // P1 diurnal principal solar.
        {0.64959E-4, 1.0,-3.0, 1.0,-0.25},  // Q1 diurnal elliptical lunar.
        {0.53234E-5, 0.0, 2.0, 0.0, 0.00},  // Mf Lunar fortnightly.
        {0.26392E-5, 0.0, 1.0,-1.0, 0.00},  // Mm Lunar monthly.
        {0.03982E-5, 2.0, 0.0, 0.0, 0.00}   // Ssa Solar semiannual.
    };
    const double ep1975[]={1975,1,1,0,0,0};
    double ep[6],fday,days,t,t2,t3,a[5],ang,dp[3]={0};
    int i,j;
    
    trace(3,"tide_oload:\n");
    
    /* angular argument: see subroutine arg.f for reference [1] */
    time2epoch(tut,ep);
    fday=ep[3]*3600.0+ep[4]*60.0+ep[5];
    ep[3]=ep[4]=ep[5]=0.0;
    days=timediff(epoch2time(ep),epoch2time(ep1975))/86400.0+1.0;
    t=(27392.500528+1.000000035*days)/36525.0;
    t2=t*t; t3=t2*t;
    
    a[0]=fday;
    a[1]=(279.69668+36000.768930485*t+3.03E-4*t2)*D2R; /* H0 */
    a[2]=(270.434358+481267.88314137*t-0.001133*t2+1.9E-6*t3)*D2R; /* S0 */
    a[3]=(334.329653+4069.0340329577*t-0.010325*t2-1.2E-5*t3)*D2R; /* P0 */
    a[4]=2.0*PI;
    
    /* displacements by 11 constituents */
    for (i=0;i<11;i++) {
        ang=0.0;
        for (j=0;j<5;j++) ang+=a[j]*args[i][j];
        for (j=0;j<3;j++) dp[j]+=odisp[j+i*6]*cos(ang-odisp[j+3+i*6]*D2R);
    }
    denu[0]=-dp[1];
    denu[1]=-dp[2];
    denu[2]= dp[0];
    
    trace(5,"tide_oload: denu=%.3f %.3f %.3f\n",denu[0],denu[1],denu[2]);
}
/* iers mean pole (ref [7] eq.7.25) ------------------------------------------*/
static void iers_mean_pole(gtime_t tut, double *xp_bar, double *yp_bar)
{
    const double ep2000[]={2000,1,1,0,0,0};
    double y,y2,y3;
    
    y=timediff(tut,epoch2time(ep2000))/86400.0/365.25;
    
    if (y<3653.0/365.25) { /* until 2010.0 */
        y2=y*y; y3=y2*y;
        *xp_bar= 55.974+1.8243*y+0.18413*y2+0.007024*y3; /* (mas) */
        *yp_bar=346.346+1.7896*y-0.10729*y2-0.000908*y3;
    }
    else { /* after 2010.0 */
        *xp_bar= 23.513+7.6141*y; /* (mas) */
        *yp_bar=358.891-0.6287*y;
    }
}
/* displacement by pole tide (ref [7] eq.7.26) --------------------------------*/
static void tide_pole(gtime_t tut, const double *pos, const double *erpv,
                      double *denu)
{
    double xp_bar,yp_bar,m1,m2,cosl,sinl;
    
    trace(3,"tide_pole: pos=%.3f %.3f\n",pos[0]*R2D,pos[1]*R2D);
    
    /* iers mean pole (mas) */
    iers_mean_pole(tut,&xp_bar,&yp_bar);
    
    /* ref [7] eq.7.24 */
    m1= erpv[0]/AS2R-xp_bar*1E-3; /* (as) */
    m2=-erpv[1]/AS2R+yp_bar*1E-3;
    
    /* sin(2*theta) = sin(2*phi), cos(2*theta)=-cos(2*phi) */
    cosl=cos(pos[1]);
    sinl=sin(pos[1]);
    denu[0]=  9E-3*sin(pos[0])    *(m1*sinl-m2*cosl); /* de= Slambda (m) */
    denu[1]= -9E-3*cos(2.0*pos[0])*(m1*cosl+m2*sinl); /* dn=-Stheta  (m) */
    denu[2]=-33E-3*sin(2.0*pos[0])*(m1*cosl+m2*sinl); /* du= Sr      (m) */
    
    trace(5,"tide_pole : denu=%.3f %.3f %.3f\n",denu[0],denu[1],denu[2]);
}
/* tidal displacement ----------------------------------------------------------
* displacements by earth tides
* args   : gtime_t tutc     I   time in utc
*          double *rr       I   site position (ecef) (m)
*          int    opt       I   options (or of the followings)
*                                 1: solid earth tide
*                                 2: ocean tide loading
*                                 4: pole tide
*                                 8: elimate permanent deformation
*          double *erp      I   earth rotation parameters (NULL: not used)
*          double *odisp    I   ocean loading parameters  (NULL: not used)
*                                 odisp[0+i*6]: consituent i amplitude radial(m)
*                                 odisp[1+i*6]: consituent i amplitude west  (m)
*                                 odisp[2+i*6]: consituent i amplitude south (m)
*                                 odisp[3+i*6]: consituent i phase radial  (deg)
*                                 odisp[4+i*6]: consituent i phase west    (deg)
*                                 odisp[5+i*6]: consituent i phase south   (deg)
*                                (i=0:M2,1:S2,2:N2,3:K2,4:K1,5:O1,6:P1,7:Q1,
*                                   8:Mf,9:Mm,10:Ssa)
*          double *dr       O   displacement by earth tides (ecef) (m)
* return : none
* notes  : see ref [1], [2] chap 7
*          see ref [4] 5.2.1, 5.2.2, 5.2.3
*          ver.2.4.0 does not use ocean loading and pole tide corrections
*-----------------------------------------------------------------------------*/
extern void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp, const double *odisp,
                     double *dr) {
  char tstr[40];
  trace(3, "tidedisp: tutc=%s\n", time2str(tutc, tstr, 0));

  double erpv[5] = {0};
  if (erp) {
    geterp(erp, utc2gpst(tutc), erpv);
  }
  gtime_t tut = timeadd(tutc, erpv[2]);

  dr[0] = dr[1] = dr[2] = 0.0;

  if (norm(rr, 3) <= 0.0) return;

  double pos[2];
  pos[0] = asin(rr[2] / norm(rr, 3));
  pos[1] = atan2(rr[1], rr[0]);
  double E[9];
  xyz2enu(pos, E);

  if (opt & 1) {  // Solid earth tides.
    // Sun and moon position in ECEF.
    double rs[3], rm[3];
    sunmoonpos(tutc, erpv, rs, rm, NULL);
    double drt[3];
    dehanttideinel(tutc, (double *)rr, rs, rm, drt);
    for (int i = 0; i < 3; i++) dr[i] += drt[i];
    trace(5, "tidedisp solid: dr=%.3f %.3f %.3f\n", drt[0], drt[1], drt[2]);
  }
  if ((opt & 2) && odisp) {  // Ocean tide loading.
    double denu[3];
    tide_oload(tut, odisp, denu);
    double drt[3];
    matmul("TN", 3, 1, 3, E, denu, drt);
    for (int i = 0; i < 3; i++) dr[i] += drt[i];
  }
  if ((opt & 4) && erp) {  // Pole tide.
    double denu[3];
    tide_pole(tut, pos, erpv, denu);
    double drt[3];
    matmul("TN", 3, 1, 3, E, denu, drt);
    for (int i = 0; i < 3; i++) dr[i] += drt[i];
  }
  trace(5, "tidedisp: dr=%.3f %.3f %.3f\n", dr[0], dr[1], dr[2]);
}
