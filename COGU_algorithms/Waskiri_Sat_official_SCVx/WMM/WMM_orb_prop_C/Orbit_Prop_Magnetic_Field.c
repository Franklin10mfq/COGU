/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: Orbit_Prop_Magnetic_Field.c
 *
 * Code generated for Simulink model 'Orbit_Prop_Magnetic_Field'.
 *
 * Model version                  : 1.1040
 * Simulink Coder version         : 25.2 (R2025b) 28-Jul-2025
 * C/C++ source code generated on : Thu Apr  2 21:45:48 2026
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: ARM Compatible->ARM Cortex-M
 * Code generation objectives:
 *    1. ROM efficiency
 *    2. RAM efficiency
 * Validation result: Not run
 */

#include "Orbit_Prop_Magnetic_Field.h"
#include "Orbit_Prop_Magnetic_Field_types.h"
#include "rtwtypes.h"
#include "rt_nonfinite.h"
#include <math.h>
#include "Orbit_Prop_Magnetic_Field_private.h"
#include <string.h>
#include "rt_assert.h"
#include "zero_crossing_types.h"
#include "rt_defines.h"
#include <float.h>

/* Block signals (default storage) */
rt2B rt2B_n;

/* Block states (default storage) */
rt2DW rt2DW_l;

/* Previous zero-crossings (trigger) states */
rt2PrevZCX rt2PrevZCX_p;

/* External inputs (root inport signals with default storage) */
rt2ExtU rt2U;

/* External outputs (root outports fed by signals with default storage) */
rt2ExtY rt2Y;

/* Real-time model */
static rt2RT_MODEL rt2M_;
rt2RT_MODEL *const rt2M = &rt2M_;

/* Forward declaration for local functions */
static real_T rt2norm(const real_T rt2x[3]);
static real_T rt2sum_i(const real_T rt2x[34]);
static real_T rt2sum(const real_T rt2x[64]);
static real_T rt2sum_iq(const real_T rt2x[20]);
static real_T rt2sum_iqe(const real_T rt2x[7]);
static real_T rt2sum_iqej(const real_T rt2x[3]);
static real_T rt2sum_iqejq(const real_T rt2x[5]);
static real_T rt2sum_iqejqxh(const real_T rt2x[10]);
static real_T rt2sum_iqejqx(const real_T rt2x[40]);
static real_T rt2sum_iqejqxho(const real_T rt2x[6]);
static real_T rt2sum_iqejqxhoh(const real_T rt2x[63]);
static void rt2accel(const rt2captured_var_1 *rt2b_const_GM_earth, const
                     rt2captured_var_1 *rt2b_H_base, const rt2captured_var_1
                     *rt2b_H_scale, const rt2captured_var_1 *rt2b_rho_o, const
                     rt2captured_var_1 *rt2b_atmospheric_drag_terms,
                     rt2captured_var *rt2orig2sun, const rt2captured_var_1
                     *rt2JD, const rt2captured_var_1
                     *rt2b_sun_radiation_pressure_terms, const rt2captured_var_1
                     *rt2b_const_GM_moon, const rt2captured_var_1
                     *rt2b_const_GM_sun, const rt2captured_var
                     *rt2accel_thruster, const real_T rt2r[3], const real_T
                     rt2v[3], real_T rt2ddr[3]);
real_T rt_atan2d_snf(real_T rt2u0, real_T rt2u1)
{
  real_T rt2y;
  if (rtIsNaN(rt2u0) || rtIsNaN(rt2u1)) {
    rt2y = (rtNaN);
  } else if (rtIsInf(rt2u0) && rtIsInf(rt2u1)) {
    int32_T rt2tmp;
    int32_T rt2tmp_0;
    if (rt2u0 > 0.0) {
      rt2tmp = 1;
    } else {
      rt2tmp = -1;
    }

    if (rt2u1 > 0.0) {
      rt2tmp_0 = 1;
    } else {
      rt2tmp_0 = -1;
    }

    rt2y = atan2(rt2tmp, rt2tmp_0);
  } else if (rt2u1 == 0.0) {
    if (rt2u0 > 0.0) {
      rt2y = RT_PI / 2.0;
    } else if (rt2u0 < 0.0) {
      rt2y = -(RT_PI / 2.0);
    } else {
      rt2y = 0.0;
    }
  } else {
    rt2y = atan2(rt2u0, rt2u1);
  }

  return rt2y;
}

real_T rt_modd_snf(real_T rt2u0, real_T rt2u1)
{
  real_T rt2y;
  rt2y = rt2u0;
  if (rt2u1 == 0.0) {
    if (rt2u0 == 0.0) {
      rt2y = rt2u1;
    }
  } else if (rtIsNaN(rt2u0) || rtIsNaN(rt2u1) || rtIsInf(rt2u0)) {
    rt2y = (rtNaN);
  } else if (rt2u0 == 0.0) {
    rt2y = 0.0 / rt2u1;
  } else if (rtIsInf(rt2u1)) {
    if ((rt2u1 < 0.0) != (rt2u0 < 0.0)) {
      rt2y = rt2u1;
    }
  } else {
    boolean_T rt2yEq;
    rt2y = fmod(rt2u0, rt2u1);
    rt2yEq = (rt2y == 0.0);
    if ((!rt2yEq) && (rt2u1 > floor(rt2u1))) {
      real_T rt2q;
      rt2q = fabs(rt2u0 / rt2u1);
      rt2yEq = !(fabs(rt2q - floor(rt2q + 0.5)) > DBL_EPSILON * rt2q);
    }

    if (rt2yEq) {
      rt2y = rt2u1 * 0.0;
    } else if ((rt2u0 < 0.0) != (rt2u1 < 0.0)) {
      rt2y += rt2u1;
    }
  }

  return rt2y;
}

real_T rt_remd_snf(real_T rt2u0, real_T rt2u1)
{
  real_T rt2y;
  if (rtIsNaN(rt2u0) || rtIsNaN(rt2u1) || rtIsInf(rt2u0)) {
    rt2y = (rtNaN);
  } else if (rtIsInf(rt2u1)) {
    rt2y = rt2u0;
  } else if ((rt2u1 != 0.0) && (rt2u1 != trunc(rt2u1))) {
    real_T rt2q;
    rt2q = fabs(rt2u0 / rt2u1);
    if (!(fabs(rt2q - floor(rt2q + 0.5)) > DBL_EPSILON * rt2q)) {
      rt2y = 0.0 * rt2u0;
    } else {
      rt2y = fmod(rt2u0, rt2u1);
    }
  } else {
    rt2y = fmod(rt2u0, rt2u1);
  }

  return rt2y;
}

/* Function for MATLAB Function: '<S1>/HPOP_RK4' */
static real_T rt2norm(const real_T rt2x[3])
{
  real_T rt2absxk;
  real_T rt2scale;
  real_T rt2t;
  real_T rt2y;
  rt2scale = 3.3121686421112381E-170;
  rt2absxk = fabs(rt2x[0]);
  if (rt2absxk > 3.3121686421112381E-170) {
    rt2y = 1.0;
    rt2scale = rt2absxk;
  } else {
    rt2t = rt2absxk / 3.3121686421112381E-170;
    rt2y = rt2t * rt2t;
  }

  rt2absxk = fabs(rt2x[1]);
  if (rt2absxk > rt2scale) {
    rt2t = rt2scale / rt2absxk;
    rt2y = rt2y * rt2t * rt2t + 1.0;
    rt2scale = rt2absxk;
  } else {
    rt2t = rt2absxk / rt2scale;
    rt2y += rt2t * rt2t;
  }

  rt2absxk = fabs(rt2x[2]);
  if (rt2absxk > rt2scale) {
    rt2t = rt2scale / rt2absxk;
    rt2y = rt2y * rt2t * rt2t + 1.0;
    rt2scale = rt2absxk;
  } else {
    rt2t = rt2absxk / rt2scale;
    rt2y += rt2t * rt2t;
  }

  rt2y = rt2scale * sqrt(rt2y);
  if (rtIsNaN(rt2y)) {
    int32_T rt2k;
    rt2k = 0;
    int32_T exitg1;
    do {
      exitg1 = 0;
      if (rt2k < 3) {
        if (rtIsNaN(rt2x[rt2k])) {
          exitg1 = 1;
        } else {
          rt2k++;
        }
      } else {
        rt2y = (rtInf);
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return rt2y;
}

real_T rt_powd_snf(real_T rt2u0, real_T rt2u1)
{
  real_T rt2y;
  if (rtIsNaN(rt2u0) || rtIsNaN(rt2u1)) {
    rt2y = (rtNaN);
  } else {
    real_T rt2tmp;
    real_T rt2tmp_0;
    rt2tmp = fabs(rt2u0);
    rt2tmp_0 = fabs(rt2u1);
    if (rtIsInf(rt2u1)) {
      if (rt2tmp == 1.0) {
        rt2y = 1.0;
      } else if (rt2tmp > 1.0) {
        if (rt2u1 > 0.0) {
          rt2y = (rtInf);
        } else {
          rt2y = 0.0;
        }
      } else if (rt2u1 > 0.0) {
        rt2y = 0.0;
      } else {
        rt2y = (rtInf);
      }
    } else if (rt2tmp_0 == 0.0) {
      rt2y = 1.0;
    } else if (rt2tmp_0 == 1.0) {
      if (rt2u1 > 0.0) {
        rt2y = rt2u0;
      } else {
        rt2y = 1.0 / rt2u0;
      }
    } else if (rt2u1 == 2.0) {
      rt2y = rt2u0 * rt2u0;
    } else if ((rt2u1 == 0.5) && (rt2u0 >= 0.0)) {
      rt2y = sqrt(rt2u0);
    } else if ((rt2u0 < 0.0) && (rt2u1 > floor(rt2u1))) {
      rt2y = (rtNaN);
    } else {
      rt2y = pow(rt2u0, rt2u1);
    }
  }

  return rt2y;
}

/* Function for MATLAB Function: '<S1>/HPOP_RK4' */
static real_T rt2sum_i(const real_T rt2x[34])
{
  real_T rt2y;
  int32_T rt2k;
  rt2y = rt2x[0];
  for (rt2k = 0; rt2k < 33; rt2k++) {
    rt2y += rt2x[rt2k + 1];
  }

  return rt2y;
}

/* Function for MATLAB Function: '<S1>/HPOP_RK4' */
static real_T rt2sum(const real_T rt2x[64])
{
  real_T rt2y;
  int32_T rt2k;
  rt2y = rt2x[0];
  for (rt2k = 0; rt2k < 63; rt2k++) {
    rt2y += rt2x[rt2k + 1];
  }

  return rt2y;
}

/* Function for MATLAB Function: '<S1>/HPOP_RK4' */
static real_T rt2sum_iq(const real_T rt2x[20])
{
  real_T rt2y;
  int32_T rt2k;
  rt2y = rt2x[0];
  for (rt2k = 0; rt2k < 19; rt2k++) {
    rt2y += rt2x[rt2k + 1];
  }

  return rt2y;
}

/* Function for MATLAB Function: '<S1>/HPOP_RK4' */
static real_T rt2sum_iqe(const real_T rt2x[7])
{
  real_T rt2y;
  int32_T rt2k;
  rt2y = rt2x[0];
  for (rt2k = 0; rt2k < 6; rt2k++) {
    rt2y += rt2x[rt2k + 1];
  }

  return rt2y;
}

/* Function for MATLAB Function: '<S1>/HPOP_RK4' */
static real_T rt2sum_iqej(const real_T rt2x[3])
{
  return (rt2x[0] + rt2x[1]) + rt2x[2];
}

/* Function for MATLAB Function: '<S1>/HPOP_RK4' */
static real_T rt2sum_iqejq(const real_T rt2x[5])
{
  return (((rt2x[0] + rt2x[1]) + rt2x[2]) + rt2x[3]) + rt2x[4];
}

/* Function for MATLAB Function: '<S1>/HPOP_RK4' */
static real_T rt2sum_iqejqxh(const real_T rt2x[10])
{
  real_T rt2y;
  int32_T rt2k;
  rt2y = rt2x[0];
  for (rt2k = 0; rt2k < 9; rt2k++) {
    rt2y += rt2x[rt2k + 1];
  }

  return rt2y;
}

/* Function for MATLAB Function: '<S1>/HPOP_RK4' */
static real_T rt2sum_iqejqx(const real_T rt2x[40])
{
  real_T rt2y;
  int32_T rt2k;
  rt2y = rt2x[0];
  for (rt2k = 0; rt2k < 39; rt2k++) {
    rt2y += rt2x[rt2k + 1];
  }

  return rt2y;
}

/* Function for MATLAB Function: '<S1>/HPOP_RK4' */
static real_T rt2sum_iqejqxho(const real_T rt2x[6])
{
  real_T rt2y;
  int32_T rt2k;
  rt2y = rt2x[0];
  for (rt2k = 0; rt2k < 5; rt2k++) {
    rt2y += rt2x[rt2k + 1];
  }

  return rt2y;
}

/* Function for MATLAB Function: '<S1>/HPOP_RK4' */
static real_T rt2sum_iqejqxhoh(const real_T rt2x[63])
{
  real_T rt2y;
  int32_T rt2k;
  rt2y = rt2x[0];
  for (rt2k = 0; rt2k < 62; rt2k++) {
    rt2y += rt2x[rt2k + 1];
  }

  return rt2y;
}

/* Function for MATLAB Function: '<S1>/HPOP_RK4' */
static void rt2accel(const rt2captured_var_1 *rt2b_const_GM_earth, const
                     rt2captured_var_1 *rt2b_H_base, const rt2captured_var_1
                     *rt2b_H_scale, const rt2captured_var_1 *rt2b_rho_o, const
                     rt2captured_var_1 *rt2b_atmospheric_drag_terms,
                     rt2captured_var *rt2orig2sun, const rt2captured_var_1
                     *rt2JD, const rt2captured_var_1
                     *rt2b_sun_radiation_pressure_terms, const rt2captured_var_1
                     *rt2b_const_GM_moon, const rt2captured_var_1
                     *rt2b_const_GM_sun, const rt2captured_var
                     *rt2accel_thruster, const real_T rt2r[3], const real_T
                     rt2v[3], real_T rt2ddr[3])
{
  real_T rt2b_x[64];
  real_T rt2i[64];
  real_T rt2bb[63];
  real_T rt2cb[63];
  real_T rt2eb_a[63];
  real_T rt2g_x[40];
  real_T rt2w[40];
  real_T rt2c_x[34];
  real_T rt2j[34];
  real_T rt2d_x[20];
  real_T rt2l[20];
  real_T rt2h_x[10];
  real_T rt2y[10];
  real_T rt2e_x[7];
  real_T rt2m[7];
  real_T rt2ab[6];
  real_T rt2i_x[6];
  real_T rt2p[5];
  real_T rt2Ehp[3];
  real_T rt2c_c_tmp[3];
  real_T rt2earth_heliocentric_position[3];
  real_T rt2D_sun;
  real_T rt2J_ephem_cen;
  real_T rt2M_moon;
  real_T rt2T;
  real_T rt2a;
  real_T rt2a_tmp;
  real_T rt2c_a;
  real_T rt2ddr_tmp;
  real_T rt2ddr_tmp_tmp;
  real_T rt2ddr_tmp_tmp_0;
  real_T rt2earth_heliocentric_position_latitude;
  real_T rt2earth_heliocentric_position_longitude;
  real_T rt2earth_heliocentric_position_longitude_tmp;
  real_T rt2earth_heliocentric_position_longitude_tmp_0;
  real_T rt2earth_heliocentric_position_longitude_tmp_1;
  real_T rt2i_a;
  real_T rt2i_a_tmp;
  real_T rt2i_a_tmp_0;
  real_T rt2j_a;
  real_T rt2julian_idx_3;
  real_T rt2julian_idx_4;
  real_T rt2n_tmp;
  real_T rt2sun_geocentric_position_latitude;
  real_T rt2u_idx_3;
  int32_T rt2i_0;
  int32_T rt2k;
  static const real_T rt2u_a[64] = { 0.0, 6283.07585, 12566.1517, 5753.3849,
    3.5231, 77713.7715, 7860.4194, 3930.2097, 11506.7698, 529.691, 1577.3435,
    5884.927, 26.298, 398.149, 5223.694, 5507.553, 18849.228, 775.523, 0.067,
    11790.629, 796.298, 10977.079, 5486.778, 2544.314, 5573.143, 6069.777,
    213.299, 2942.463, 20.775, 0.98, 4694.003, 15720.839, 7.114, 2146.17, 155.42,
    161000.69, 6275.96, 71430.7, 17260.15, 12036.46, 5088.63, 3154.69, 801.82,
    9437.76, 8827.39, 7084.9, 6286.6, 14143.5, 6279.55, 12139.55, 1748.02,
    5856.48, 1194.45, 8429.24, 19651.05, 10447.39, 10213.29, 1059.38, 2352.87,
    6812.77, 17789.85, 83996.85, 1349.87, 4690.48 };

  static const real_T rt2d[64] = { 0.0, 4.6692568, 4.6261, 2.7441, 2.8289,
    3.6277, 4.4181, 6.1352, 0.7425, 2.0371, 1.1096, 5.233, 2.045, 3.508, 1.179,
    2.533, 4.583, 4.205, 2.92, 5.849, 1.899, 0.315, 0.345, 4.806, 1.869, 2.4458,
    0.833, 3.411, 1.083, 0.645, 0.636, 0.976, 4.267, 6.21, 0.68, 5.98, 1.3, 3.67,
    1.81, 3.04, 1.76, 3.5, 4.68, 0.83, 3.98, 1.82, 2.78, 4.39, 3.47, 0.19, 1.33,
    0.28, 0.49, 5.37, 2.4, 6.17, 6.04, 2.57, 1.71, 1.78, 0.59, 0.44, 2.74, 3.16
  };

  static const real_T rt2v_a[34] = { 0.0, 6283.07585, 12566.1517, 3.523, 26.298,
    1577.344, 18849.23, 529.69, 398.15, 5507.55, 5223.69, 155.42, 796.3, 775.52,
    7.11, 0.98, 5486.78, 213.3, 6275.96, 2544.31, 2146.17, 10977.08, 1748.02,
    5088.63, 1194.45, 4694.0, 553.57, 3286.6, 1349.87, 242.73, 951.72, 2352.87,
    9437.76, 4690.48 };

  static const real_T rt2e[34] = { 0.0, 2.678235, 2.6351, 1.59, 5.796, 2.966,
    2.59, 1.14, 1.87, 4.41, 2.89, 2.17, 0.4, 0.47, 2.65, 5.34, 1.85, 4.97, 2.99,
    0.03, 1.43, 1.21, 2.83, 3.26, 5.27, 2.08, 0.77, 1.3, 4.24, 2.7, 5.64, 5.3,
    2.65, 4.67 };

  static const real_T rt2w_a[20] = { 0.0, 6283.0758, 12566.152, 3.52, 26.3,
    155.42, 18849.23, 77713.77, 775.52, 1577.34, 7.11, 5573.14, 796.3, 5507.55,
    242.73, 529.69, 398.15, 553.57, 5223.69, 0.98 };

  static const real_T rt2f[20] = { 0.0, 1.0721, 0.867, 0.05, 5.19, 3.68, 0.76,
    2.06, 0.83, 4.66, 1.03, 3.44, 5.14, 6.05, 1.19, 6.12, 0.31, 2.28, 4.38, 3.75
  };

  static const real_T rt2x_a[7] = { 6283.076, 0.0, 12566.15, 155.42, 3.52,
    18849.23, 242.73 };

  static const real_T rt2g[7] = { 5.844, 0.0, 5.49, 5.2, 4.72, 5.3, 5.97 };

  static const int32_T rt2i_1[64] = { 175347046, 3341656, 34894, 3497, 3418,
    3136, 2676, 2343, 1324, 1273, 1199, 990, 902, 857, 780, 753, 505, 492, 357,
    317, 284, 271, 243, 206, 205, 202, 156, 132, 126, 115, 103, 102, 102, 99, 98,
    86, 85, 85, 80, 79, 71, 74, 74, 70, 62, 61, 57, 56, 56, 52, 52, 51, 49, 41,
    41, 39, 37, 37, 36, 36, 33, 30, 30, 25 };

  static const real_T rt2j_0[34] = { 6.28331966747E+11, 206059.0, 4303.0, 425.0,
    119.0, 109.0, 93.0, 72.0, 68.0, 67.0, 59.0, 56.0, 45.0, 36.0, 29.0, 21.0,
    19.0, 19.0, 17.0, 16.0, 16.0, 15.0, 12.0, 12.0, 12.0, 12.0, 11.0, 10.0, 10.0,
    9.0, 9.0, 8.0, 6.0, 6.0 };

  static const uint16_T rt2l_0[20] = { 52919U, 8720U, 309U, 27U, 16U, 16U, 10U,
    9U, 7U, 5U, 4U, 4U, 3U, 3U, 3U, 3U, 3U, 3U, 2U, 2U };

  static const int16_T rt2m_0[7] = { 289, 35, 17, 3, 1, 1, 1 };

  static const int16_T rt2p_0[5] = { 280, 102, 80, 44, 32 };

  static const real_T rt2ab_a[5] = { 84334.662, 5507.553, 5223.69, 2352.87,
    1577.34 };

  static const real_T rt2o[5] = { 3.199, 5.422, 3.88, 3.7, 4.0 };

  static const real_T rt2bb_a[40] = { 0.0, 6283.07585, 12566.1517, 77713.7715,
    5753.3849, 7860.4194, 11506.77, 3930.21, 5884.927, 5507.553, 5223.694,
    5573.143, 11790.629, 1577.344, 10977.079, 18849.228, 5486.778, 6069.78,
    15720.84, 161000.69, 17260.15, 529.69, 83996.85, 71430.7, 2544.31, 775.52,
    9437.76, 6275.96, 4694.0, 8827.39, 19651.05, 12139.55, 12036.46, 2942.46,
    7084.9, 5088.63, 398.15, 6286.6, 6279.55, 10447.39 };

  static const real_T rt2q[40] = { 0.0, 3.0984635, 3.05525, 5.1985, 1.1739,
    2.8469, 5.453, 4.564, 3.661, 0.964, 5.9, 0.299, 4.273, 5.847, 5.022, 3.012,
    5.055, 0.89, 5.69, 1.27, 0.27, 0.92, 2.01, 5.24, 3.25, 2.58, 5.54, 6.01,
    5.36, 2.39, 0.83, 4.9, 1.67, 1.84, 0.24, 0.18, 1.78, 1.21, 1.9, 4.59 };

  static const real_T rt2cb_a[10] = { 6283.07585, 12566.1517, 0.0, 18849.23,
    5507.55, 5223.69, 1577.34, 10977.08, 6275.96, 5486.78 };

  static const real_T rt2s[10] = { 1.10749, 1.0644, 3.142, 1.02, 2.84, 1.32,
    1.42, 5.91, 1.42, 0.27 };

  static const real_T rt2db_a[6] = { 6283.0758, 12566.152, 0.0, 77713.77,
    5573.14, 18849.0 };

  static const real_T rt2t[6] = { 5.7846, 5.579, 3.14, 3.63, 1.87, 5.47 };

  static const int32_T rt2w_0[40] = { 100013989, 1670700, 13956, 3084, 1628,
    1576, 925, 542, 472, 346, 329, 307, 243, 212, 186, 175, 110, 98, 86, 86, 85,
    63, 57, 56, 49, 47, 45, 43, 39, 38, 37, 37, 36, 35, 33, 32, 32, 28, 28, 26 };

  static const int32_T rt2y_0[10] = { 103019, 1721, 702, 32, 31, 25, 18, 10, 9,
    9 };

  static const int16_T rt2ab_0[6] = { 4359, 124, 12, 9, 6, 3 };

  static const int8_T rt2eb_a_0[315] = { 0, -2, 0, 0, 0, 0, -2, 0, 0, -2, -2, -2,
    0, 2, 0, 2, 0, 0, -2, 0, 2, 0, 0, -2, 0, -2, 0, 0, 2, -2, 0, -2, 0, 0, 2, 2,
    0, -2, 0, 2, 2, -2, -2, 2, 2, 0, -2, -2, 0, -2, -2, 0, -1, -2, 1, 0, 0, -1,
    0, 0, 2, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 1, 0, -1, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0, 0, 0,
    -1, -1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, -1, 1, -1, -1, 0, -1, 0, 0, 0, 0, 0, 1,
    0, 0, 1, 0, 1, 0, -1, 0, 1, -1, -1, 1, 2, -2, 0, 2, 2, 1, 0, 0, -1, 0, -1, 0,
    0, 1, 0, 2, -1, 1, 0, 1, 0, 0, 1, 2, 1, -2, 0, 1, 0, 0, 2, 2, 0, 1, 1, 0, 0,
    1, -2, 1, 1, 1, -1, 3, 0, 0, 2, 2, 0, 0, 0, 2, 2, 2, 2, 0, 2, 2, 0, 0, 2, 0,
    2, 0, 2, 2, 2, 0, 2, 2, 2, 2, 0, 0, 2, 0, 0, 0, -2, 2, 2, 2, 0, 2, 2, 0, 2,
    2, 0, 0, 0, 2, 0, 2, 0, 2, -2, 0, 0, 0, 2, 2, 0, 0, 2, 2, 2, 2, 1, 2, 2, 2,
    0, 0, 2, 1, 2, 2, 0, 1, 2, 0, 1, 2, 1, 1, 0, 1, 2, 2, 0, 2, 0, 0, 1, 0, 1, 2,
    1, 1, 1, 0, 1, 2, 2, 0, 2, 1, 0, 2, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
    2, 0, 0, 2, 2, 2, 2 };

  static const real_T rt2fb_a[63] = { 8.9, -3.1, -0.5, 0.5, -0.1, 0.0, -0.6, 0.0,
    -0.1, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const int32_T rt2bb_0[63] = { 92025, 5736, 977, -895, 54, -7, 224, 200,
    129, -95, 0, -70, -53, 0, -33, 26, 32, 27, 0, -24, 16, 13, 0, -12, 0, 0, -10,
    0, -8, 7, 9, 7, 6, 0, 5, 3, -3, 0, 3, 3, 0, -3, -3, 3, 3, 0, 3, 3, 3, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  static const real_T rt2gb_a[63] = { -174.2, -1.6, -0.2, 0.2, -3.4, 0.1, 1.2,
    -0.4, 0.0, -0.5, 0.0, 0.1, 0.0, 0.0, 0.1, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const int32_T rt2cb_0[63] = { -171996, -13187, -2274, 2062, 1426, 712,
    -517, -386, -301, 217, -158, 129, 123, 63, 63, -59, -58, -51, 48, 46, -38,
    -31, 29, 29, 26, -22, 21, 17, 16, -16, -15, -13, -12, 11, -10, -8, 7, -7, -7,
    -7, 6, 6, 6, -6, -6, 5, -5, -5, -5, 4, 4, 4, -4, -4, -4, 3, -3, -3, -3, -3,
    -3, -3, -3 };

  rt2a_tmp = rt2norm(rt2r);
  rt2a = -rt2b_const_GM_earth->contents / rt_powd_snf(rt2a_tmp, 3.0);
  rt2u_idx_3 = rt2a_tmp;
  rt2i_a_tmp_0 = 6.378135E+6 / rt2a_tmp;
  rt2j_a = rt2r[2] / rt2a_tmp;
  rt2i_a_tmp = rt2b_const_GM_earth->contents / (rt2a_tmp * rt2a_tmp);
  rt2i_a = rt2i_a_tmp * -0.0016239 * (rt2i_a_tmp_0 * rt2i_a_tmp_0);
  rt2c_a = rt2i_a_tmp * 1.26E-6 * rt_powd_snf(rt2i_a_tmp_0, 3.0);
  rt2i_a_tmp_0 = rt2i_a_tmp * 1.0062500000000001E-6 * rt_powd_snf(rt2i_a_tmp_0,
    5.0);
  rt2a_tmp = exp(-((rt2a_tmp - 6.3781E+6) - rt2b_H_base->contents) /
                 rt2b_H_scale->contents) * rt2b_rho_o->contents * (-0.5 *
    rt2b_atmospheric_drag_terms->contents);
  rt2i_a_tmp = rt2norm(rt2v);
  rt2J_ephem_cen = ((rt2JD->contents + 0.00037249999999999995) - 2.451545E+6) /
    36525.0;
  rt2julian_idx_3 = rt2J_ephem_cen;
  rt2julian_idx_4 = rt2J_ephem_cen / 10.0;
  for (rt2k = 0; rt2k < 64; rt2k++) {
    rt2b_x[rt2k] = cos(rt2u_a[rt2k] * rt2julian_idx_4 + rt2d[rt2k]);
  }

  for (rt2k = 0; rt2k < 34; rt2k++) {
    rt2c_x[rt2k] = cos(rt2v_a[rt2k] * rt2julian_idx_4 + rt2e[rt2k]);
  }

  for (rt2k = 0; rt2k < 20; rt2k++) {
    rt2d_x[rt2k] = cos(rt2w_a[rt2k] * rt2julian_idx_4 + rt2f[rt2k]);
  }

  for (rt2k = 0; rt2k < 7; rt2k++) {
    rt2e_x[rt2k] = cos(rt2x_a[rt2k] * rt2julian_idx_4 + rt2g[rt2k]);
  }

  rt2J_ephem_cen = rt_powd_snf(rt2julian_idx_4, 3.0);
  rt2sun_geocentric_position_latitude = rt_powd_snf(rt2julian_idx_4, 4.0);
  for (rt2k = 0; rt2k < 64; rt2k++) {
    rt2i[rt2k] = (real_T)rt2i_1[rt2k] * rt2b_x[rt2k];
  }

  for (rt2k = 0; rt2k < 34; rt2k++) {
    rt2j[rt2k] = rt2j_0[rt2k] * rt2c_x[rt2k];
  }

  for (rt2k = 0; rt2k < 20; rt2k++) {
    rt2l[rt2k] = (real_T)rt2l_0[rt2k] * rt2d_x[rt2k];
  }

  for (rt2k = 0; rt2k < 7; rt2k++) {
    rt2m[rt2k] = (real_T)rt2m_0[rt2k] * rt2e_x[rt2k];
  }

  rt2earth_heliocentric_position_longitude = 0.0 * rt2julian_idx_4;
  rt2Ehp[0] = cos(rt2earth_heliocentric_position_longitude + 3.142) * 114.0;
  rt2n_tmp = 6283.08 * rt2julian_idx_4;
  rt2Ehp[1] = cos(rt2n_tmp + 4.13) * 8.0;
  rt2T = 12566.15 * rt2julian_idx_4;
  rt2Ehp[2] = cos(rt2T + 3.84);
  rt2earth_heliocentric_position_longitude_tmp = rt2julian_idx_4 *
    rt2julian_idx_4;
  rt2earth_heliocentric_position_longitude = (((((rt2sum_i(rt2j) *
    rt2julian_idx_4 + rt2sum(rt2i)) +
    rt2earth_heliocentric_position_longitude_tmp * rt2sum_iq(rt2l)) + rt2sum_iqe
    (rt2m) * rt2J_ephem_cen) + rt2sum_iqej(rt2Ehp) *
    rt2sun_geocentric_position_latitude) + cos
    (rt2earth_heliocentric_position_longitude + 3.14) * rt_powd_snf
    (rt2julian_idx_4, 5.0)) / 1.0E+8 * 180.0 / 3.1415926535897931;
  rt2earth_heliocentric_position_longitude -= floor
    (rt2earth_heliocentric_position_longitude / 360.0) * 360.0;
  if (rt2earth_heliocentric_position_longitude < 0.0) {
    rt2earth_heliocentric_position_longitude += 360.0;
  }

  for (rt2k = 0; rt2k < 5; rt2k++) {
    rt2p[rt2k] = cos(rt2ab_a[rt2k] * rt2julian_idx_4 + rt2o[rt2k]) * (real_T)
      rt2p_0[rt2k];
  }

  rt2earth_heliocentric_position_latitude = ((cos(5507.55 * rt2julian_idx_4 +
    3.9) * 9.0 + cos(5223.69 * rt2julian_idx_4 + 1.73) * 6.0) * rt2julian_idx_4
    + rt2sum_iqejq(rt2p)) / 1.0E+8 * 180.0 / 3.1415926535897931;
  rt2earth_heliocentric_position_latitude -= floor
    (rt2earth_heliocentric_position_latitude / 360.0) * 360.0;
  if (rt2earth_heliocentric_position_latitude < 0.0) {
    rt2earth_heliocentric_position_latitude += 360.0;
  }

  for (rt2k = 0; rt2k < 40; rt2k++) {
    rt2g_x[rt2k] = cos(rt2bb_a[rt2k] * rt2julian_idx_4 + rt2q[rt2k]);
  }

  for (rt2k = 0; rt2k < 10; rt2k++) {
    rt2h_x[rt2k] = cos(rt2cb_a[rt2k] * rt2julian_idx_4 + rt2s[rt2k]);
  }

  for (rt2k = 0; rt2k < 6; rt2k++) {
    rt2i_x[rt2k] = cos(rt2db_a[rt2k] * rt2julian_idx_4 + rt2t[rt2k]);
  }

  for (rt2k = 0; rt2k < 40; rt2k++) {
    rt2w[rt2k] = (real_T)rt2w_0[rt2k] * rt2g_x[rt2k];
  }

  for (rt2k = 0; rt2k < 10; rt2k++) {
    rt2y[rt2k] = (real_T)rt2y_0[rt2k] * rt2h_x[rt2k];
  }

  for (rt2k = 0; rt2k < 6; rt2k++) {
    rt2ab[rt2k] = (real_T)rt2ab_0[rt2k] * rt2i_x[rt2k];
  }

  rt2earth_heliocentric_position[0] = (((cos(6283.076 * rt2julian_idx_4 + 4.273)
    * 145.0 + cos(rt2T + 3.92) * 7.0) * rt2J_ephem_cen + ((rt2sum_iqejqxh(rt2y) *
    rt2julian_idx_4 + rt2sum_iqejqx(rt2w)) +
    rt2earth_heliocentric_position_longitude_tmp * rt2sum_iqejqxho(rt2ab))) +
    cos(rt2n_tmp + 2.56) * 4.0 * rt2sun_geocentric_position_latitude) / 1.0E+8;
  rt2J_ephem_cen = (rt2earth_heliocentric_position_longitude + 180.0) - floor
    ((rt2earth_heliocentric_position_longitude + 180.0) / 360.0) * 360.0;
  if (rt2J_ephem_cen < 0.0) {
    rt2J_ephem_cen += 360.0;
  }

  rt2sun_geocentric_position_latitude = -rt2earth_heliocentric_position_latitude
    - floor(-rt2earth_heliocentric_position_latitude / 360.0) * 360.0;
  if (rt2sun_geocentric_position_latitude < 0.0) {
    rt2sun_geocentric_position_latitude += 360.0;
  }

  rt2earth_heliocentric_position_latitude = rt2J_ephem_cen;
  rt2J_ephem_cen = rt_powd_snf(rt2julian_idx_3, 3.0);
  rt2n_tmp = rt2julian_idx_3 * rt2julian_idx_3;
  rt2p[0] = ((rt2n_tmp * -0.0019142 + 5.2777689814961416E-6 * rt2J_ephem_cen) +
             445267.11148 * rt2julian_idx_3) + 297.85036;
  rt2p[1] = ((rt2n_tmp * -0.0001603 + -3.3333333333333333E-6 * rt2J_ephem_cen) +
             35999.05034 * rt2julian_idx_3) + 357.52772;
  rt2p[2] = ((rt2n_tmp * 0.0086972 + 1.7777777777777777E-5 * rt2J_ephem_cen) +
             477198.867398 * rt2julian_idx_3) + 134.96298;
  rt2p[3] = ((rt2n_tmp * -0.0036825 + 3.0555810187307116E-6 * rt2J_ephem_cen) +
             483202.017538 * rt2julian_idx_3) + 93.27191;
  rt2p[4] = ((rt2n_tmp * 0.0020708 + 2.2222222222222221E-6 * rt2J_ephem_cen) +
             -1934.136261 * rt2julian_idx_3) + 125.04452;
  memset(&rt2eb_a[0], 0, 63U * sizeof(real_T));
  for (rt2k = 0; rt2k < 5; rt2k++) {
    rt2n_tmp = rt2p[rt2k];
    for (rt2i_0 = 0; rt2i_0 < 63; rt2i_0++) {
      rt2eb_a[rt2i_0] += (real_T)rt2eb_a_0[63 * rt2k + rt2i_0] * rt2n_tmp;
    }
  }

  rt2J_ephem_cen = rt2julian_idx_4 / 10.0;
  for (rt2k = 0; rt2k < 63; rt2k++) {
    rt2earth_heliocentric_position_longitude = rt2eb_a[rt2k] *
      3.1415926535897931 / 180.0;
    rt2bb[rt2k] = (rt2fb_a[rt2k] * rt2julian_idx_3 + (real_T)rt2bb_0[rt2k]) *
      cos(rt2earth_heliocentric_position_longitude);
    rt2cb[rt2k] = (rt2gb_a[rt2k] * rt2julian_idx_3 + (real_T)rt2cb_0[rt2k]) *
      sin(rt2earth_heliocentric_position_longitude);
  }

  rt2earth_heliocentric_position_longitude = (((((((((((2.45 * rt_powd_snf
    (rt2J_ephem_cen, 10.0) + 5.79 * rt_powd_snf(rt2J_ephem_cen, 9.0)) + 27.87 *
    rt_powd_snf(rt2J_ephem_cen, 8.0)) + 7.12 * rt_powd_snf(rt2J_ephem_cen, 7.0))
    + -39.05 * rt_powd_snf(rt2J_ephem_cen, 6.0)) + -249.67 * rt_powd_snf
    (rt2J_ephem_cen, 5.0)) + -51.38 * rt_powd_snf(rt2J_ephem_cen, 4.0)) +
    1999.25 * rt_powd_snf(rt2J_ephem_cen, 3.0)) + rt2J_ephem_cen *
    rt2J_ephem_cen * -1.55) + -4680.93 * rt2J_ephem_cen) + 84381.448) / 3600.0 +
    rt2sum_iqejqxhoh(rt2bb) / 3.6E+7) * 3.1415926535897931 / 180.0;
  rt2J_ephem_cen = cos(rt2earth_heliocentric_position_longitude);
  rt2earth_heliocentric_position_longitude = sin
    (rt2earth_heliocentric_position_longitude);
  rt2julian_idx_3 = ((rt2sum_iqejqxhoh(rt2cb) / 3.6E+7 +
                      rt2earth_heliocentric_position_latitude) + -20.4898 /
                     (3600.0 * rt2earth_heliocentric_position[0])) *
    3.1415926535897931 / 180.0;
  rt2earth_heliocentric_position_latitude = sin(rt2julian_idx_3);
  rt2julian_idx_4 = rt2sun_geocentric_position_latitude * 3.1415926535897931 /
    180.0;
  rt2sun_geocentric_position_latitude = rt_atan2d_snf
    (rt2earth_heliocentric_position_latitude * rt2J_ephem_cen - tan
     (rt2julian_idx_4) * rt2earth_heliocentric_position_longitude, cos
     (rt2julian_idx_3)) * 180.0 / 3.1415926535897931;
  rt2sun_geocentric_position_latitude -= floor
    (rt2sun_geocentric_position_latitude / 360.0) * 360.0;
  if (rt2sun_geocentric_position_latitude < 0.0) {
    rt2sun_geocentric_position_latitude += 360.0;
  }

  rt2n_tmp = asin(cos(rt2julian_idx_4) *
                  rt2earth_heliocentric_position_longitude *
                  rt2earth_heliocentric_position_latitude + sin(rt2julian_idx_4)
                  * rt2J_ephem_cen);
  rt2Ehp[2] = rt2sun_geocentric_position_latitude;
  rt2J_ephem_cen = rt2earth_heliocentric_position[0] * 1.49597870691E+11;
  rt2sun_geocentric_position_latitude = rt2n_tmp * 180.0 / 3.1415926535897931 *
    3.1415926535897931 / 180.0;
  rt2earth_heliocentric_position_longitude = rt2Ehp[2] * 3.1415926535897931 /
    180.0;
  rt2n_tmp = rt2J_ephem_cen * cos(rt2sun_geocentric_position_latitude);
  rt2orig2sun->contents[0] = rt2n_tmp * cos
    (rt2earth_heliocentric_position_longitude);
  rt2orig2sun->contents[1] = rt2n_tmp * sin
    (rt2earth_heliocentric_position_longitude);
  rt2orig2sun->contents[2] = rt2J_ephem_cen * sin
    (rt2sun_geocentric_position_latitude);
  rt2sun_geocentric_position_latitude = sqrt((rt2r[0] * rt2r[0] + rt2r[1] *
    rt2r[1]) + rt2r[2] * rt2r[2]);
  rt2earth_heliocentric_position_longitude = 6.371E+6 /
    rt2sun_geocentric_position_latitude;
  if (rt2earth_heliocentric_position_longitude > 1.0) {
    rt2earth_heliocentric_position_longitude = 1.0;
  }

  rt2julian_idx_3 = rt2orig2sun->contents[0] - rt2r[0];
  rt2earth_heliocentric_position[0] = rt2julian_idx_3;
  rt2julian_idx_4 = rt2orig2sun->contents[1] - rt2r[1];
  rt2earth_heliocentric_position[1] = rt2julian_idx_4;
  rt2n_tmp = rt2orig2sun->contents[2] - rt2r[2];
  rt2earth_heliocentric_position[2] = rt2n_tmp;
  rt2J_ephem_cen = rt2norm(rt2earth_heliocentric_position);
  rt2earth_heliocentric_position_latitude = 1.495978707E+11 / rt2J_ephem_cen;
  rt2sun_geocentric_position_latitude = -(real_T)(acos(-((rt2r[0] *
    rt2orig2sun->contents[0] + rt2r[1] * rt2orig2sun->contents[1]) + rt2r[2] *
    rt2orig2sun->contents[2]) / (sqrt((rt2orig2sun->contents[0] *
    rt2orig2sun->contents[0] + rt2orig2sun->contents[1] * rt2orig2sun->contents
    [1]) + rt2orig2sun->contents[2] * rt2orig2sun->contents[2]) *
    rt2sun_geocentric_position_latitude)) > asin
    (rt2earth_heliocentric_position_longitude)) *
    rt2b_sun_radiation_pressure_terms->contents *
    (rt2earth_heliocentric_position_latitude *
     rt2earth_heliocentric_position_latitude);
  rt2T = ((rt2JD->contents - 2.451545E+6) + 0.00037249999999999995) / 36525.0;
  rt2earth_heliocentric_position_longitude = rt_powd_snf(rt2T, 3.0);
  rt2D_sun = rt2T * rt2T;
  rt2M_moon = ((8328.6914228849528 * rt2T + 2.355548394) + rt2D_sun *
               0.000151795) + 3.103E-7 *
    rt2earth_heliocentric_position_longitude;
  rt2earth_heliocentric_position_latitude = ((8433.4661583190045 * rt2T +
    1.627901934) - rt2D_sun * 6.4272E-5) + 5.34E-8 *
    rt2earth_heliocentric_position_longitude;
  rt2earth_heliocentric_position_longitude_tmp = (((7771.3771461739689 * rt2T +
    5.198469514) - rt2D_sun * 3.34086E-5) + 9.22E-8 *
    rt2earth_heliocentric_position_longitude) * 2.0;
  rt2earth_heliocentric_position_longitude_tmp_0 = rt2M_moon -
    rt2earth_heliocentric_position_longitude_tmp;
  rt2earth_heliocentric_position_longitude_tmp_1 = 2.0 * rt2M_moon;
  rt2earth_heliocentric_position_longitude = (((((((481267.8813 * rt2T + 218.32)
    + 6.29 * sin(rt2M_moon)) - sin
    (rt2earth_heliocentric_position_longitude_tmp_0) * 1.27) + sin
    (rt2earth_heliocentric_position_longitude_tmp) * 0.66) + sin
    (rt2earth_heliocentric_position_longitude_tmp_1) * 0.21) - sin
    (((628.30195601077912 * rt2T + 6.2400359) - rt2D_sun * 2.7974E-5) - 5.81E-8 *
     rt2earth_heliocentric_position_longitude) * 0.19) - sin(2.0 *
    rt2earth_heliocentric_position_latitude) * 0.11) * 3.1415926535897931 /
    180.0;
  rt2earth_heliocentric_position_latitude = (((sin(rt2M_moon +
    rt2earth_heliocentric_position_latitude) * 0.28 + 5.13 * sin
    (rt2earth_heliocentric_position_latitude)) - sin
    (rt2earth_heliocentric_position_latitude - rt2M_moon) * 0.28) - sin
    (rt2earth_heliocentric_position_latitude -
     rt2earth_heliocentric_position_longitude_tmp) * 0.17) * 3.1415926535897931 /
    180.0;
  rt2T = (23.439291 - 0.0130042 * rt2T) * 3.1415926535897931 / 180.0;
  rt2M_moon = 1.0 / sin(((((0.0518 * cos(rt2M_moon) + 0.9508) + cos
    (rt2earth_heliocentric_position_longitude_tmp_0) * 0.0095) + cos
    (rt2earth_heliocentric_position_longitude_tmp) * 0.0078) + cos
    (rt2earth_heliocentric_position_longitude_tmp_1) * 0.0028) *
                        3.1415926535897931 / 180.0) * 6.371008E+6;
  rt2D_sun = cos(rt2earth_heliocentric_position_latitude);
  rt2earth_heliocentric_position_longitude_tmp = sin(rt2T);
  rt2earth_heliocentric_position_longitude_tmp_0 = sin
    (rt2earth_heliocentric_position_longitude);
  rt2T = cos(rt2T);
  rt2earth_heliocentric_position_latitude = sin
    (rt2earth_heliocentric_position_latitude);
  rt2Ehp[0] = rt2D_sun * cos(rt2earth_heliocentric_position_longitude) *
    rt2M_moon;
  rt2Ehp[1] = (rt2T * rt2D_sun * rt2earth_heliocentric_position_longitude_tmp_0
               - rt2earth_heliocentric_position_longitude_tmp *
               rt2earth_heliocentric_position_latitude) * rt2M_moon;
  rt2Ehp[2] = (rt2earth_heliocentric_position_longitude_tmp * rt2D_sun *
               rt2earth_heliocentric_position_longitude_tmp_0 + rt2T *
               rt2earth_heliocentric_position_latitude) * rt2M_moon;
  rt2earth_heliocentric_position_latitude = rt_powd_snf(rt2norm(rt2Ehp), 3.0);
  rt2earth_heliocentric_position[0] = rt2Ehp[0] - rt2r[0];
  rt2c_c_tmp[0] = rt2julian_idx_3;
  rt2earth_heliocentric_position[1] = rt2Ehp[1] - rt2r[1];
  rt2c_c_tmp[1] = rt2julian_idx_4;
  rt2earth_heliocentric_position[2] = rt2Ehp[2] - rt2r[2];
  rt2c_c_tmp[2] = rt2n_tmp;
  rt2earth_heliocentric_position_longitude = rt_powd_snf(rt2norm
    (rt2earth_heliocentric_position), 3.0);
  rt2T = rt_powd_snf(rt2norm(rt2c_c_tmp), 3.0);
  rt2M_moon = rt_powd_snf(rt2norm(rt2orig2sun->contents), 3.0);
  rt2earth_heliocentric_position_longitude_tmp = rt_powd_snf(rt2j_a, 4.0);
  rt2D_sun = rt2r[0] / rt2u_idx_3;
  rt2earth_heliocentric_position_longitude_tmp_0 = rt2j_a * rt2j_a;
  rt2earth_heliocentric_position_longitude_tmp_1 = (7.0 * rt_powd_snf(rt2j_a,
    3.0) - rt2j_a * 3.0) * 5.0;
  rt2ddr_tmp_tmp = rt2earth_heliocentric_position_longitude_tmp_0 * 5.0;
  rt2ddr_tmp_tmp_0 = 63.0 * rt2earth_heliocentric_position_longitude_tmp;
  rt2ddr_tmp = (3.0 - rt2earth_heliocentric_position_longitude_tmp_0 * 42.0) +
    rt2ddr_tmp_tmp_0;
  rt2ddr[0] = (((((((1.0 - rt2ddr_tmp_tmp) * rt2D_sun * rt2i_a +
                    rt2earth_heliocentric_position_longitude_tmp_1 * rt2D_sun *
                    rt2c_a) + rt2ddr_tmp * rt2D_sun * rt2i_a_tmp_0) + (rt2a *
    rt2r[0] + rt2accel_thruster->contents[0])) + rt2a_tmp * rt2v[0] * rt2i_a_tmp)
                + rt2julian_idx_3 / rt2J_ephem_cen *
                rt2sun_geocentric_position_latitude) +
               (rt2earth_heliocentric_position[0] /
                rt2earth_heliocentric_position_longitude - rt2Ehp[0] /
                rt2earth_heliocentric_position_latitude) *
               rt2b_const_GM_moon->contents) + (rt2julian_idx_3 / rt2T -
    rt2orig2sun->contents[0] / rt2M_moon) * rt2b_const_GM_sun->contents;
  rt2D_sun = rt2r[1] / rt2u_idx_3;
  rt2ddr[1] = (((((((1.0 - rt2ddr_tmp_tmp) * rt2D_sun * rt2i_a +
                    rt2earth_heliocentric_position_longitude_tmp_1 * rt2D_sun *
                    rt2c_a) + rt2ddr_tmp * rt2D_sun * rt2i_a_tmp_0) + (rt2a *
    rt2r[1] + rt2accel_thruster->contents[1])) + rt2a_tmp * rt2v[1] * rt2i_a_tmp)
                + rt2julian_idx_4 / rt2J_ephem_cen *
                rt2sun_geocentric_position_latitude) +
               (rt2earth_heliocentric_position[1] /
                rt2earth_heliocentric_position_longitude - rt2Ehp[1] /
                rt2earth_heliocentric_position_latitude) *
               rt2b_const_GM_moon->contents) + (rt2julian_idx_4 / rt2T -
    rt2orig2sun->contents[1] / rt2M_moon) * rt2b_const_GM_sun->contents;
  rt2ddr[2] = ((((((((rt2earth_heliocentric_position_longitude_tmp_0 * 10.0 -
                      11.666666666666666 *
                      rt2earth_heliocentric_position_longitude_tmp) - 1.0) * 3.0
                    * rt2c_a + (3.0 - rt2ddr_tmp_tmp) * rt2j_a * rt2i_a) +
                   -((15.0 - rt2earth_heliocentric_position_longitude_tmp_0 *
                      70.0) + rt2ddr_tmp_tmp_0) * rt2j_a * rt2i_a_tmp_0) + (rt2a
    * rt2r[2] + rt2accel_thruster->contents[2])) + rt2a_tmp * rt2v[2] *
                 rt2i_a_tmp) + rt2n_tmp / rt2J_ephem_cen *
                rt2sun_geocentric_position_latitude) +
               (rt2earth_heliocentric_position[2] /
                rt2earth_heliocentric_position_longitude - rt2Ehp[2] /
                rt2earth_heliocentric_position_latitude) *
               rt2b_const_GM_moon->contents) + (rt2n_tmp / rt2T -
    rt2orig2sun->contents[2] / rt2M_moon) * rt2b_const_GM_sun->contents;
}

/* Model step function */
void Orbit_Prop_Magnetic_Field_step(void)
{
  rt2captured_var rt2accel_thruster;
  rt2captured_var rt2orig2sun;
  rt2captured_var_1 rt2JD;
  rt2captured_var_1 rt2b_H_base;
  rt2captured_var_1 rt2b_H_scale;
  rt2captured_var_1 rt2b_atmospheric_drag_terms;
  rt2captured_var_1 rt2b_const_GM_earth;
  rt2captured_var_1 rt2b_const_GM_moon;
  rt2captured_var_1 rt2b_const_GM_sun;
  rt2captured_var_1 rt2b_rho_o;
  rt2captured_var_1 rt2b_sun_radiation_pressure_terms;
  real_T rt2b_tc_old[169];
  real_T rt2b_TmpSignalConversionAtppnInport1[13];
  real_T rt2b_Assignment[11];
  real_T rt2b_Assignment1[11];
  real_T rt2b_R_ECEF2ECI[9];
  real_T rt2b_Transpose[9];
  real_T rt2rt2b_Transpose[9];
  real_T rt2a2[3];
  real_T rt2a2_tmp[3];
  real_T rt2a3[3];
  real_T rt2a3_tmp[3];
  real_T rt2a4[3];
  real_T rt2a4_tmp[3];
  real_T rt2b_MatrixMultiply2[3];
  real_T rt2tmp[3];
  real_T rt2E4;
  real_T rt2E6;
  real_T rt2E7;
  real_T rt2b_Abs1;
  real_T rt2b_ECEF_Position_to_LLA1_o2;
  real_T rt2b_Selector1_f;
  real_T rt2b_Sum_n;
  real_T rt2b_x1;
  real_T rt2c;
  real_T rt2q;
  real_T rt2rt2Merge;
  real_T rt2rt2Merge_c;
  real_T rt2rt2Merge_c_tmp;
  real_T rt2rt2UnitDelay2_DSTATE_idx_2;
  real_T rt2rt2UnitDelay2_DSTATE_idx_3;
  real_T rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0;
  real_T rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1;
  real_T rt2rt2b_Sum1_c_idx_3;
  real_T rt2rt2b_sincos_o1_idx_1;
  real_T rt2rt2b_sincos_o1_j_idx_0;
  real_T rt2rt2b_sincos_o1_n_idx_1;
  real_T rt2rt2b_sincos_o2_idx_1;
  int32_T rt2i;
  int32_T rt2rt2Sum1_tmp;
  int32_T rt2rt2b_Product_tmp;
  int32_T rt2s48_iter;
  int32_T rt2s89_iter;
  int8_T rt2rt2b_R_ECEF2ECI_tmp[3];
  uint8_T rt2ForIterator_IterationMarker[6];
  boolean_T rt2b_LogicalOperator;
  boolean_T rt2b_max_relop_h;
  static const int8_T rt2c_0[3] = { 1, 0, 0 };

  static const real_T rt2b_b[9] = { 0.99999995201534, -0.000309789149845, 0.0,
    0.000309789149845, 0.99999995201534, 0.0, 0.0, 0.0, 1.0 };

  static const real_T rt2c_b[9] = { 0.99999992651175, 0.0, 0.000383375135609,
    0.0, 1.0, 0.0, -0.000383375135609, 0.0, 0.99999992651175 };

  static const real_T rt2d_b[9] = { 1.0, 0.0, 0.0, 0.0, 0.999999999999749,
    -7.087976E-7, 0.0, 7.087976E-7, 0.999999999999749 };

  /* MATLAB Function: '<S5>/MATLAB Function' incorporates:
   *  Inport: '<Root>/JD'
   */
  rt2c = ((rt2U.rt2JD - 2.451545E+6) * 360.98564736629 + 280.46061837) *
    0.017453292519943295;
  if (rtIsNaN(rt2c)) {
    rt2q = (rtNaN);
  } else if (rtIsInf(rt2c)) {
    rt2q = (rtNaN);
  } else {
    rt2q = fabs(rt2c / 6.2831853071795862);
    if (fabs(rt2q - floor(rt2q + 0.5)) > 2.2204460492503131E-16 * rt2q) {
      rt2q = fmod(rt2c, 6.2831853071795862);
    } else {
      rt2q = 0.0;
    }

    if (rt2q == 0.0) {
      rt2q = 0.0;
    } else if (rt2q < 0.0) {
      rt2q += 6.2831853071795862;
    }
  }

  rt2c = cos(rt2q);
  rt2q = sin(rt2q);
  rt2b_R_ECEF2ECI[0] = rt2c;
  rt2b_R_ECEF2ECI[3] = -rt2q;
  rt2b_R_ECEF2ECI[6] = 0.0;
  rt2b_R_ECEF2ECI[1] = rt2q;
  rt2b_R_ECEF2ECI[4] = rt2c;
  rt2b_R_ECEF2ECI[7] = 0.0;
  rt2rt2b_R_ECEF2ECI_tmp[0] = 0;
  rt2b_R_ECEF2ECI[2] = 0.0;
  rt2rt2b_R_ECEF2ECI_tmp[1] = 0;
  rt2b_R_ECEF2ECI[5] = 0.0;
  rt2rt2b_R_ECEF2ECI_tmp[2] = 1;
  rt2b_R_ECEF2ECI[8] = 1.0;

  /* End of MATLAB Function: '<S5>/MATLAB Function' */

  /* Math: '<S1>/Transpose3' incorporates:
   *  Math: '<S4>/Transpose'
   */
  for (rt2s89_iter = 0; rt2s89_iter < 3; rt2s89_iter++) {
    rt2b_Transpose[3 * rt2s89_iter] = rt2b_R_ECEF2ECI[rt2s89_iter];
    rt2b_Transpose[3 * rt2s89_iter + 1] = rt2b_R_ECEF2ECI[rt2s89_iter + 3];
    rt2b_Transpose[3 * rt2s89_iter + 2] = rt2b_R_ECEF2ECI[rt2s89_iter + 6];
  }

  /* End of Math: '<S1>/Transpose3' */

  /* Delay: '<S1>/Delay2' incorporates:
   *  Inport: '<Root>/r_flag'
   *  Inport: '<Root>/r_pred_init'
   */
  if ((rt_ZCFcn(ANY_ZERO_CROSSING,&rt2PrevZCX_p.rt2Delay2_Reset_ZCE,
                (rt2U.rt2r_flag)) != NO_ZCEVENT) || (rt2U.rt2r_flag != 0.0)) {
    rt2DW_l.rt2icLoad = true;
  }

  if (rt2DW_l.rt2icLoad) {
    rt2DW_l.rt2Delay2_DSTATE[0] = rt2U.rt2r_pred_init[0];
    rt2DW_l.rt2Delay2_DSTATE[1] = rt2U.rt2r_pred_init[1];
    rt2DW_l.rt2Delay2_DSTATE[2] = rt2U.rt2r_pred_init[2];
  }

  /* Product: '<S1>/Matrix Multiply2' incorporates:
   *  Delay: '<S1>/Delay2'
   *  Math: '<S4>/Transpose'
   */
  rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 = 0.0;
  rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 = 0.0;
  rt2rt2b_sincos_o1_n_idx_1 = 0.0;
  for (rt2s89_iter = 0; rt2s89_iter < 3; rt2s89_iter++) {
    rt2b_Sum_n = rt2DW_l.rt2Delay2_DSTATE[rt2s89_iter];
    rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 += rt2b_Transpose[3 * rt2s89_iter] *
      rt2b_Sum_n;
    rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 += rt2b_Transpose[3 * rt2s89_iter + 1]
      * rt2b_Sum_n;
    rt2rt2b_sincos_o1_n_idx_1 += rt2b_Transpose[3 * rt2s89_iter + 2] *
      rt2b_Sum_n;
  }

  /* ECEF2LLA: '<S1>/ECEF_Position_to_LLA1' incorporates:
   *  Product: '<S1>/Matrix Multiply2'
   */
  rt2c = sqrt(rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 *
              rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 +
              rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 *
              rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1);
  rt2b_Abs1 = rt_atan2d_snf(rt2rt2b_sincos_o1_n_idx_1, 0.99664718933525254 *
    rt2c);
  rt2q = sin(rt2b_Abs1);
  rt2rt2Merge = cos(rt2b_Abs1);
  rt2q = rt_atan2d_snf(42841.311513313573 * rt2q * rt2q * rt2q +
                       rt2rt2b_sincos_o1_n_idx_1, rt2c - 42697.672707179969 *
                       rt2rt2Merge * rt2rt2Merge * rt2rt2Merge);
  rt2rt2Merge = rt_atan2d_snf(0.99664718933525254 * sin(rt2q), cos(rt2q));
  rt2s89_iter = 0;
  while ((rt2b_Abs1 != rt2rt2Merge) && (rt2s89_iter < 5)) {
    rt2b_Abs1 = rt2rt2Merge;
    rt2q = sin(rt2rt2Merge);
    rt2rt2Merge = cos(rt2rt2Merge);
    rt2q = rt_atan2d_snf(42841.311513313573 * rt2q * rt2q * rt2q +
                         rt2rt2b_sincos_o1_n_idx_1, rt2c - 42697.672707179969 *
                         rt2rt2Merge * rt2rt2Merge * rt2rt2Merge);
    rt2rt2Merge = rt_atan2d_snf(0.99664718933525254 * sin(rt2q), cos(rt2q));
    rt2s89_iter++;
  }

  rt2b_ECEF_Position_to_LLA1_o2 = fabs(rt2q);
  rt2b_Abs1 = rt2q;
  rt2rt2Merge = rt_atan2d_snf(rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1,
    rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0);
  if (rt2b_ECEF_Position_to_LLA1_o2 > 3.1415926535897931) {
    rt2b_Abs1 = rt_modd_snf(rt2q + 3.1415926535897931, 6.2831853071795862) -
      3.1415926535897931;
    rt2b_ECEF_Position_to_LLA1_o2 = fabs(rt2b_Abs1);
  }

  if (rt2b_ECEF_Position_to_LLA1_o2 > 1.5707963267948966) {
    rt2rt2Merge += 3.1415926535897931;
    if (rtIsNaN(rt2b_Abs1)) {
      rt2b_Sum_n = (rtNaN);
    } else if (rt2b_Abs1 < 0.0) {
      rt2b_Sum_n = -1.0;
    } else {
      rt2b_Sum_n = (rt2b_Abs1 > 0.0);
    }

    rt2b_Abs1 = (1.5707963267948966 - (rt2b_ECEF_Position_to_LLA1_o2 -
      1.5707963267948966)) * rt2b_Sum_n;
  }

  if (fabs(rt2rt2Merge) > 3.1415926535897931) {
    rt2rt2Merge = rt_remd_snf(rt2rt2Merge, 6.2831853071795862);
    rt2rt2Merge -= trunc(rt2rt2Merge / 3.1415926535897931) * 6.2831853071795862;
  }

  rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 = rt2b_Abs1 * 180.0 /
    3.1415926535897931;
  rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 = rt2rt2Merge * 180.0 /
    3.1415926535897931;
  rt2b_Abs1 = sin(rt2q);
  rt2rt2Merge_c = 6.378137E+6 / sqrt(1.0 - rt2b_Abs1 * rt2b_Abs1 *
    0.0066943799901413165);
  rt2b_ECEF_Position_to_LLA1_o2 = ((rt2rt2Merge_c * 0.0066943799901413165 *
    rt2b_Abs1 + rt2rt2b_sincos_o1_n_idx_1) * rt2b_Abs1 + rt2c * cos(rt2q)) -
    rt2rt2Merge_c;

  /* End of ECEF2LLA: '<S1>/ECEF_Position_to_LLA1' */

  /* Switch: '<S35>/Switch' incorporates:
   *  Abs: '<S35>/Abs'
   *  Bias: '<S35>/Bias'
   *  Bias: '<S35>/Bias1'
   *  Constant: '<S35>/Constant2'
   *  Constant: '<S38>/Constant'
   *  Math: '<S35>/Math Function1'
   *  RelationalOperator: '<S38>/Compare'
   */
  if (fabs(rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0) > 180.0) {
    rt2c = rt_modd_snf(rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 + 180.0, 360.0) -
      180.0;
  } else {
    rt2c = rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0;
  }

  /* End of Switch: '<S35>/Switch' */

  /* Abs: '<S31>/Abs1' */
  rt2b_Abs1 = fabs(rt2c);

  /* RelationalOperator: '<S33>/Compare' incorporates:
   *  Constant: '<S33>/Constant'
   */
  rt2b_max_relop_h = (rt2b_Abs1 > 90.0);

  /* Switch: '<S25>/Switch1' incorporates:
   *  Constant: '<S25>/Constant'
   *  Constant: '<S25>/Constant1'
   */
  if (rt2b_max_relop_h) {
    rt2s89_iter = 180;
  } else {
    rt2s89_iter = 0;
  }

  /* Sum: '<S25>/Sum' incorporates:
   *  Switch: '<S25>/Switch1'
   */
  rt2q = (real_T)rt2s89_iter + rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1;

  /* Switch: '<S32>/Switch' incorporates:
   *  Abs: '<S32>/Abs'
   *  Bias: '<S32>/Bias'
   *  Bias: '<S32>/Bias1'
   *  Constant: '<S32>/Constant2'
   *  Constant: '<S39>/Constant'
   *  Math: '<S32>/Math Function1'
   *  RelationalOperator: '<S39>/Compare'
   */
  if (fabs(rt2q) > 180.0) {
    rt2q = rt_modd_snf(rt2q + 180.0, 360.0) - 180.0;
  }

  /* End of Switch: '<S32>/Switch' */

  /* If: '<S34>/If' incorporates:
   *  Abs: '<S34>/Abs'
   */
  if (rt2b_Abs1 > 0.0) {
    /* Outputs for IfAction SubSystem: '<S34>/If Action Subsystem' incorporates:
     *  ActionPort: '<S36>/Action Port'
     */
    /* Merge: '<S34>/Merge' incorporates:
     *  Product: '<S36>/Divide'
     */
    rt2rt2Merge = rt2c / rt2b_Abs1;

    /* End of Outputs for SubSystem: '<S34>/If Action Subsystem' */
  } else {
    /* Outputs for IfAction SubSystem: '<S34>/If Action Subsystem1' incorporates:
     *  ActionPort: '<S37>/Action Port'
     */
    /* Merge: '<S34>/Merge' incorporates:
     *  Constant: '<S37>/Constant'
     *  SignalConversion generated from: '<S37>/zero'
     */
    rt2rt2Merge = 0.0;

    /* End of Outputs for SubSystem: '<S34>/If Action Subsystem1' */
  }

  /* End of If: '<S34>/If' */

  /* Switch: '<S31>/Switch' incorporates:
   *  Bias: '<S31>/Bias'
   *  Bias: '<S31>/Bias1'
   *  Gain: '<S31>/Gain'
   *  Product: '<S31>/Divide1'
   */
  if (rt2b_max_relop_h) {
    rt2c = (-(rt2b_Abs1 - 90.0) + 90.0) * rt2rt2Merge;
  }

  /* End of Switch: '<S31>/Switch' */

  /* UnitConversion: '<S90>/Unit Conversion' */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  rt2rt2Merge = 0.017453292519943295 * rt2q;
  rt2rt2b_sincos_o1_n_idx_1 = 0.017453292519943295 * rt2c;

  /* Trigonometry: '<S43>/sincos' */
  rt2b_Abs1 = sin(rt2rt2Merge);
  rt2rt2Merge = cos(rt2rt2Merge);
  rt2rt2b_sincos_o1_idx_1 = sin(rt2rt2b_sincos_o1_n_idx_1);
  rt2rt2b_sincos_o2_idx_1 = cos(rt2rt2b_sincos_o1_n_idx_1);

  /* Outputs for Enabled SubSystem: '<S29>/Convert from geodetic to  spherical coordinates ' incorporates:
   *  EnablePort: '<S42>/Enable'
   */
  /* RelationalOperator: '<S45>/Relational Operator' incorporates:
   *  Memory: '<S45>/olon'
   */
  if (rt2q != rt2DW_l.rt2olon_PreviousInput) {
    /* Outputs for Iterator SubSystem: '<S42>/For Iterator Subsystem' incorporates:
     *  ForIterator: '<S89>/For Iterator'
     */
    for (rt2s89_iter = 0; rt2s89_iter < 11; rt2s89_iter++) {
      /* Switch: '<S89>/cp[m-1] sp[m-1]' incorporates:
       *  SignalConversion generated from: '<S42>/cp[2]'
       *  SignalConversion generated from: '<S42>/sp[2]'
       *  UnitDelay: '<S89>/Unit Delay1'
       */
      if (rt2s89_iter > 0) {
        rt2rt2b_sincos_o1_n_idx_1 = rt2DW_l.rt2UnitDelay1_DSTATE[0];
        rt2rt2b_sincos_o1_j_idx_0 = rt2DW_l.rt2UnitDelay1_DSTATE[1];
      } else {
        rt2rt2b_sincos_o1_n_idx_1 = rt2rt2Merge;
        rt2rt2b_sincos_o1_j_idx_0 = rt2b_Abs1;
      }

      /* End of Switch: '<S89>/cp[m-1] sp[m-1]' */

      /* Sum: '<S89>/Sum2' incorporates:
       *  Product: '<S89>/Product1'
       *  Product: '<S89>/Product2'
       *  SignalConversion generated from: '<S42>/cp[2]'
       *  SignalConversion generated from: '<S42>/sp[2]'
       */
      rt2b_Sum_n = rt2rt2b_sincos_o1_n_idx_1 * rt2b_Abs1 +
        rt2rt2b_sincos_o1_j_idx_0 * rt2rt2Merge;

      /* Sum: '<S89>/Sum1' incorporates:
       *  Product: '<S89>/Product3'
       *  Product: '<S89>/Product8'
       *  SignalConversion generated from: '<S42>/cp[2]'
       *  SignalConversion generated from: '<S42>/sp[2]'
       */
      rt2rt2b_sincos_o1_n_idx_1 = rt2rt2b_sincos_o1_n_idx_1 * rt2rt2Merge -
        rt2rt2b_sincos_o1_j_idx_0 * rt2b_Abs1;

      /* Assignment: '<S89>/Assignment' incorporates:
       *  Assignment: '<S89>/Assignment1'
       *  Constant: '<S89>/Constant'
       *  Constant: '<S89>/Constant1'
       */
      if (rt2s89_iter == 0) {
        for (rt2i = 0; rt2i < 11; rt2i++) {
          rt2rt2b_sincos_o1_j_idx_0 = rt2ConstP_d.rt2pooled13[rt2i];
          rt2b_Assignment[rt2i] = rt2rt2b_sincos_o1_j_idx_0;
          rt2b_Assignment1[rt2i] = rt2rt2b_sincos_o1_j_idx_0;
        }
      }

      rt2b_Assignment[rt2s89_iter] = rt2b_Sum_n;

      /* End of Assignment: '<S89>/Assignment' */

      /* Assignment: '<S89>/Assignment1' */
      rt2b_Assignment1[rt2s89_iter] = rt2rt2b_sincos_o1_n_idx_1;

      /* Update for UnitDelay: '<S89>/Unit Delay1' */
      rt2DW_l.rt2UnitDelay1_DSTATE[0] = rt2rt2b_sincos_o1_n_idx_1;
      rt2DW_l.rt2UnitDelay1_DSTATE[1] = rt2b_Sum_n;
    }

    /* End of Outputs for SubSystem: '<S42>/For Iterator Subsystem' */

    /* SignalConversion generated from: '<S42>/cp[13]' incorporates:
     *  Constant: '<S42>/cp[1]'
     *  SignalConversion generated from: '<S42>/cp[2]'
     */
    rt2B_n.rt2OutportBufferForcp13[0] = 1.0;
    rt2B_n.rt2OutportBufferForcp13[1] = rt2rt2Merge;

    /* SignalConversion generated from: '<S42>/sp[13]' incorporates:
     *  Constant: '<S42>/sp[1]'
     *  SignalConversion generated from: '<S42>/sp[2]'
     */
    rt2B_n.rt2OutportBufferForsp13[0] = 0.0;
    rt2B_n.rt2OutportBufferForsp13[1] = rt2b_Abs1;

    /* SignalConversion generated from: '<S42>/cp[13]' */
    memcpy(&rt2B_n.rt2OutportBufferForcp13[2], &rt2b_Assignment1[0], 11U *
           sizeof(real_T));

    /* SignalConversion generated from: '<S42>/sp[13]' */
    memcpy(&rt2B_n.rt2OutportBufferForsp13[2], &rt2b_Assignment[0], 11U * sizeof
           (real_T));
  }

  /* End of RelationalOperator: '<S45>/Relational Operator' */
  /* End of Outputs for SubSystem: '<S29>/Convert from geodetic to  spherical coordinates ' */

  /* Math: '<S4>/Transpose' */
  for (rt2s89_iter = 0; rt2s89_iter < 3; rt2s89_iter++) {
    rt2rt2b_Transpose[3 * rt2s89_iter] = rt2b_Transpose[rt2s89_iter];
    rt2rt2b_Transpose[3 * rt2s89_iter + 1] = rt2b_Transpose[rt2s89_iter + 3];
    rt2rt2b_Transpose[3 * rt2s89_iter + 2] = rt2b_Transpose[rt2s89_iter + 6];
  }

  memcpy(&rt2b_Transpose[0], &rt2rt2b_Transpose[0], 9U * sizeof(real_T));

  /* End of Math: '<S4>/Transpose' */

  /* Outport: '<Root>/R_ECEF2ECI' */
  memcpy(&rt2Y.rt2R_ECEF2ECI[0], &rt2b_R_ECEF2ECI[0], 9U * sizeof(real_T));

  /* UnitConversion: '<S18>/Unit Conversion' */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  rt2b_Abs1 = 0.017453292519943295 * rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0;

  /* Trigonometry: '<S6>/sincos' */
  rt2rt2b_sincos_o1_j_idx_0 = cos(rt2b_Abs1);
  rt2rt2Merge = sin(rt2b_Abs1);

  /* UnitConversion: '<S18>/Unit Conversion' */
  rt2b_Abs1 = 0.017453292519943295 * rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1;

  /* Trigonometry: '<S6>/sincos' */
  rt2rt2b_sincos_o1_n_idx_1 = cos(rt2b_Abs1);
  rt2b_Abs1 = sin(rt2b_Abs1);

  /* UnaryMinus: '<S9>/Unary Minus' incorporates:
   *  Product: '<S9>/u(1)*u(4)'
   *  Trigonometry: '<S6>/sincos'
   */
  rt2b_R_ECEF2ECI[0] = -(rt2rt2Merge * rt2rt2b_sincos_o1_n_idx_1);

  /* UnaryMinus: '<S12>/Unary Minus' incorporates:
   *  Trigonometry: '<S6>/sincos'
   */
  rt2b_R_ECEF2ECI[1] = -rt2b_Abs1;

  /* UnaryMinus: '<S15>/Unary Minus' incorporates:
   *  Product: '<S15>/u(3)*u(4)'
   */
  rt2b_R_ECEF2ECI[2] = -(rt2rt2b_sincos_o1_j_idx_0 * rt2rt2b_sincos_o1_n_idx_1);

  /* UnaryMinus: '<S10>/Unary Minus' incorporates:
   *  Product: '<S10>/u(1)*u(2)'
   *  Trigonometry: '<S6>/sincos'
   */
  rt2b_R_ECEF2ECI[3] = -(rt2rt2Merge * rt2b_Abs1);

  /* SignalConversion generated from: '<S19>/Vector Concatenate' */
  rt2b_R_ECEF2ECI[4] = rt2rt2b_sincos_o1_n_idx_1;

  /* UnaryMinus: '<S16>/Unary Minus' incorporates:
   *  Product: '<S16>/u(2)*u(3)'
   *  Trigonometry: '<S6>/sincos'
   */
  rt2b_R_ECEF2ECI[5] = -(rt2rt2b_sincos_o1_j_idx_0 * rt2b_Abs1);

  /* SignalConversion generated from: '<S19>/Vector Concatenate' */
  rt2b_R_ECEF2ECI[6] = rt2rt2b_sincos_o1_j_idx_0;

  /* SignalConversion generated from: '<S19>/Vector Concatenate' incorporates:
   *  Constant: '<S14>/Constant'
   */
  rt2b_R_ECEF2ECI[7] = 0.0;

  /* UnaryMinus: '<S17>/Unary Minus' incorporates:
   *  Trigonometry: '<S6>/sincos'
   */
  rt2b_R_ECEF2ECI[8] = -rt2rt2Merge;

  /* MATLAB Function: '<S4>/MATLAB Function1' incorporates:
   *  Inport: '<Root>/JD'
   */
  rt2b_Abs1 = (rt2U.rt2JD - 2.451545E+6) / 365.25 + 2000.0;

  /* RelationalOperator: '<S46>/Relational Operator' incorporates:
   *  Memory: '<S46>/otime'
   */
  rt2b_max_relop_h = (rt2b_Abs1 != rt2DW_l.rt2otime_PreviousInput);

  /* Gain: '<S8>/Gain' */
  rt2rt2Merge = 0.001 * rt2b_ECEF_Position_to_LLA1_o2;

  /* Logic: '<S44>/Logical Operator' incorporates:
   *  Memory: '<S44>/oalt'
   *  Memory: '<S44>/olat'
   *  RelationalOperator: '<S44>/Relational Operator'
   *  RelationalOperator: '<S44>/Relational Operator1'
   */
  rt2b_LogicalOperator = ((rt2c != rt2DW_l.rt2olat_PreviousInput) ||
    (rt2rt2Merge != rt2DW_l.rt2oalt_PreviousInput));

  /* Product: '<S43>/Product' */
  rt2rt2b_sincos_o1_n_idx_1 = rt2rt2b_sincos_o1_idx_1 * rt2rt2b_sincos_o1_idx_1;

  /* Product: '<S43>/Product1' */
  rt2rt2b_sincos_o1_j_idx_0 = rt2rt2b_sincos_o2_idx_1 * rt2rt2b_sincos_o2_idx_1;

  /* Outputs for Enabled SubSystem: '<S29>/Convert from geodetic to  spherical coordinates' incorporates:
   *  EnablePort: '<S41>/Enable'
   */
  if (rt2b_LogicalOperator) {
    /* Sqrt: '<S83>/sqrt' incorporates:
     *  Product: '<S83>/Product10'
     *  Product: '<S83>/Product9'
     *  Sum: '<S83>/Sum7'
     */
    rt2rt2Merge_c = sqrt(rt2rt2b_sincos_o1_j_idx_0 * rt2ConstB.rt2a2 +
                         rt2ConstB.rt2b2 * rt2rt2b_sincos_o1_n_idx_1);

    /* Sqrt: '<S84>/sqrt' incorporates:
     *  Product: '<S84>/Product'
     *  Sum: '<S84>/Sum'
     */
    rt2b_Selector1_f = sqrt(rt2ConstB.rt2a2 - rt2rt2b_sincos_o1_n_idx_1 *
      rt2ConstB.rt2c2);

    /* Product: '<S41>/Product1' */
    rt2b_x1 = rt2b_Selector1_f * rt2rt2Merge;

    /* Sqrt: '<S41>/sqrt' incorporates:
     *  Gain: '<S86>/Gain'
     *  Product: '<S86>/Product1'
     *  Product: '<S86>/Product6'
     *  Product: '<S86>/Product7'
     *  Product: '<S86>/Product8'
     *  Sum: '<S86>/Sum5'
     *  Sum: '<S86>/Sum6'
     */
    rt2B_n.rt2sqrt = sqrt((rt2ConstB.rt2a4 - rt2ConstB.rt2c4 *
      rt2rt2b_sincos_o1_n_idx_1) / (rt2b_Selector1_f * rt2b_Selector1_f) + (2.0 *
      rt2b_x1 + rt2rt2Merge * rt2rt2Merge));

    /* Product: '<S81>/Product11' incorporates:
     *  Sum: '<S81>/Sum8'
     */
    rt2B_n.rt2Product11 = (rt2rt2Merge + rt2rt2Merge_c) / rt2B_n.rt2sqrt;

    /* Product: '<S87>/Product12' incorporates:
     *  Product: '<S87>/Product1'
     */
    rt2B_n.rt2Product12 = rt2ConstB.rt2c2_d / (rt2B_n.rt2sqrt * rt2rt2Merge_c) *
      rt2rt2b_sincos_o2_idx_1 * rt2rt2b_sincos_o1_idx_1;

    /* Sum: '<S85>/Sum2' */
    rt2rt2Merge_c = rt2ConstB.rt2a2 + rt2b_x1;

    /* Sum: '<S85>/Sum1' */
    rt2b_x1 += rt2ConstB.rt2b2;

    /* Product: '<S82>/Product4' incorporates:
     *  Product: '<S82>/Product3'
     *  Product: '<S85>/Product1'
     *  Product: '<S85>/Product2'
     *  Sqrt: '<S82>/sqrt'
     *  Sum: '<S82>/Sum3'
     */
    rt2B_n.rt2Product4 = rt2rt2b_sincos_o1_idx_1 / sqrt(rt2rt2Merge_c *
      rt2rt2Merge_c / (rt2b_x1 * rt2b_x1) * rt2rt2b_sincos_o1_j_idx_0 +
      rt2rt2b_sincos_o1_n_idx_1);

    /* Sqrt: '<S88>/sqrt' incorporates:
     *  Constant: '<S88>/Constant'
     *  Product: '<S88>/Product5'
     *  Sum: '<S88>/Sum4'
     */
    rt2B_n.rt2sqrt_g = sqrt(1.0 - rt2B_n.rt2Product4 * rt2B_n.rt2Product4);
  }

  /* End of Outputs for SubSystem: '<S29>/Convert from geodetic to  spherical coordinates' */

  /* Product: '<S29>/aor' incorporates:
   *  Constant: '<S29>/re'
   */
  rt2rt2b_sincos_o1_n_idx_1 = 6371.2 / rt2B_n.rt2sqrt;

  /* Outputs for Iterator SubSystem: '<S29>/Compute magnetic vector in spherical coordinates' incorporates:
   *  ForIterator: '<S40>/For Iterator'
   */
  /* InitializeConditions for UnitDelay: '<S40>/Unit Delay2' */
  rt2rt2b_sincos_o1_idx_1 = 0.0;
  rt2rt2b_sincos_o2_idx_1 = 0.0;
  rt2rt2UnitDelay2_DSTATE_idx_2 = 0.0;
  rt2rt2UnitDelay2_DSTATE_idx_3 = 0.0;

  /* InitializeConditions for UnitDelay: '<S40>/Unit Delay' */
  rt2rt2Merge_c = 0.0;
  for (rt2s89_iter = 0; rt2s89_iter < 12; rt2s89_iter++) {
    /* Switch: '<S40>/ar(n)' incorporates:
     *  Product: '<S29>/ar'
     */
    if (rt2s89_iter <= 0) {
      rt2rt2Merge_c = rt2rt2b_sincos_o1_n_idx_1 * rt2rt2b_sincos_o1_n_idx_1;
    }

    /* End of Switch: '<S40>/ar(n)' */

    /* Product: '<S40>/Product8' */
    rt2rt2b_sincos_o1_j_idx_0 = rt2rt2Merge_c * rt2rt2b_sincos_o1_n_idx_1;

    /* Outputs for Iterator SubSystem: '<S40>/For Iterator Subsystem' incorporates:
     *  ForIterator: '<S48>/For Iterator'
     */
    /* InitializeConditions for UnitDelay: '<S49>/Unit Delay3' */
    rt2b_x1 = 0.0;

    /* InitializeConditions for UnitDelay: '<S49>/Unit Delay2' */
    rt2E4 = 0.0;

    /* InitializeConditions for UnitDelay: '<S49>/Unit Delay1' */
    rt2E6 = 0.0;

    /* InitializeConditions for UnitDelay: '<S49>/Unit Delay4' */
    rt2E7 = 0.0;
    for (rt2i = 0; rt2i < 6; rt2i++) {
      rt2ForIterator_IterationMarker[rt2i] = 1U;
    }

    /* Sum: '<S40>/Sum' incorporates:
     *  Constant: '<S40>/Constant'
     *  Constant: '<S48>/Constant'
     *  Constant: '<S55>/Constant'
     *  Constant: '<S55>/Constant1'
     *  Logic: '<S55>/Logical Operator'
     *  RelationalOperator: '<S55>/Relational Operator'
     *  RelationalOperator: '<S55>/Relational Operator1'
     *  Sum: '<S48>/Sum1'
     */
    for (rt2s48_iter = 1; rt2s48_iter <= rt2s89_iter + 2; rt2s48_iter++) {
      /* Outputs for Enabled SubSystem: '<S48>/Time adjust the gauss coefficients' incorporates:
       *  EnablePort: '<S51>/Enable'
       */
      if (rt2b_max_relop_h) {
        /* Assignment: '<S51>/Assignment' */
        if (rt2ForIterator_IterationMarker[4] < 2) {
          rt2ForIterator_IterationMarker[4] = 2U;

          /* Assignment: '<S51>/Assignment' incorporates:
           *  UnitDelay: '<S51>/Unit Delay'
           */
          memcpy(&rt2B_n.rt2Assignment[0], &rt2DW_l.rt2UnitDelay_DSTATE_h[0],
                 169U * sizeof(real_T));
        }

        /* Outputs for Atomic SubSystem: '<S51>/If Action Subsystem' */
        /* Selector: '<S77>/c[m][n]' incorporates:
         *  Assignment: '<S51>/Assignment'
         *  Constant: '<S51>/cd[maxdef][maxdef]'
         *  Selector: '<S77>/cd[m][n]'
         */
        rt2i = ((rt2s89_iter + 1) * 13 + rt2s48_iter) - 1;

        /* Assignment: '<S51>/Assignment' incorporates:
         *  Constant: '<S29>/epoch'
         *  Constant: '<S51>/c[maxdef][maxdef]'
         *  Constant: '<S51>/cd[maxdef][maxdef]'
         *  Product: '<S77>/Product'
         *  Selector: '<S77>/c[m][n]'
         *  Selector: '<S77>/cd[m][n]'
         *  Sum: '<S29>/Sum'
         *  Sum: '<S77>/Sum'
         */
        rt2B_n.rt2Assignment[rt2i] = (rt2b_Abs1 - 2025.0) *
          rt2ConstP_d.rt2cdmaxdefmaxdef_Value[rt2i] +
          rt2ConstP_d.rt2cmaxdefmaxdef_Value[rt2i];

        /* End of Outputs for SubSystem: '<S51>/If Action Subsystem' */

        /* Switch: '<S78>/tc_old' incorporates:
         *  UnitDelay: '<S78>/Unit Delay'
         */
        for (rt2i = 0; rt2i < 169; rt2i++) {
          if (rt2s89_iter > 0) {
            rt2b_tc_old[rt2i] = rt2DW_l.rt2UnitDelay_DSTATE_hg[rt2i];
          } else {
            rt2b_tc_old[rt2i] = 0.0;
          }
        }

        /* End of Switch: '<S78>/tc_old' */

        /* If: '<S78>/If' incorporates:
         *  Constant: '<S48>/Constant'
         *  Sum: '<S48>/Sum1'
         */
        if (rt2s48_iter - 1LL != 0LL) {
          /* Outputs for IfAction SubSystem: '<S78>/If Action Subsystem1' incorporates:
           *  ActionPort: '<S79>/Action Port'
           */
          /* Assignment: '<S79>/Assignment2' */
          if (rt2ForIterator_IterationMarker[5] < 2) {
            rt2ForIterator_IterationMarker[5] = 2U;

            /* Assignment: '<S79>/Assignment2' incorporates:
             *  Switch: '<S78>/tc_old'
             */
            memcpy(&rt2B_n.rt2Assignment2[0], &rt2b_tc_old[0], 169U * sizeof
                   (real_T));
          }

          /* Selector: '<S79>/c[m][n]' incorporates:
           *  Assignment: '<S79>/Assignment2'
           *  Constant: '<S51>/cd[maxdef][maxdef]'
           *  Selector: '<S79>/cd[m][n]'
           */
          rt2i = ((rt2s48_iter - 2) * 13 + rt2s89_iter) + 1;

          /* Assignment: '<S79>/Assignment2' incorporates:
           *  Constant: '<S29>/epoch'
           *  Constant: '<S51>/c[maxdef][maxdef]'
           *  Constant: '<S51>/cd[maxdef][maxdef]'
           *  Product: '<S79>/Product'
           *  Selector: '<S79>/c[m][n]'
           *  Selector: '<S79>/cd[m][n]'
           *  Sum: '<S29>/Sum'
           *  Sum: '<S79>/Sum'
           */
          rt2B_n.rt2Assignment2[rt2i] = (rt2b_Abs1 - 2025.0) *
            rt2ConstP_d.rt2cdmaxdefmaxdef_Value[rt2i] +
            rt2ConstP_d.rt2cmaxdefmaxdef_Value[rt2i];

          /* Gain: '<S79>/Gain' incorporates:
           *  Assignment: '<S79>/Assignment2'
           *  Merge: '<S78>/Merge'
           */
          memcpy(&rt2b_tc_old[0], &rt2B_n.rt2Assignment2[0], 169U * sizeof
                 (real_T));

          /* End of Outputs for SubSystem: '<S78>/If Action Subsystem1' */
        }

        /* End of If: '<S78>/If' */
        for (rt2i = 0; rt2i < 169; rt2i++) {
          /* Sum: '<S51>/Sum2' incorporates:
           *  Assignment: '<S51>/Assignment'
           *  Merge: '<S78>/Merge'
           */
          rt2b_Selector1_f = rt2B_n.rt2Assignment[rt2i];
          rt2b_Sum_n = rt2b_tc_old[rt2i];

          /* Sum: '<S51>/Sum2' incorporates:
           *  Assignment: '<S51>/Assignment'
           *  Merge: '<S78>/Merge'
           */
          rt2B_n.rt2Sum2_b[rt2i] = rt2b_Selector1_f + rt2b_Sum_n;

          /* Update for UnitDelay: '<S51>/Unit Delay' incorporates:
           *  Assignment: '<S51>/Assignment'
           */
          rt2DW_l.rt2UnitDelay_DSTATE_h[rt2i] = rt2b_Selector1_f;

          /* Update for UnitDelay: '<S78>/Unit Delay' incorporates:
           *  Merge: '<S78>/Merge'
           */
          rt2DW_l.rt2UnitDelay_DSTATE_hg[rt2i] = rt2b_Sum_n;
        }
      }

      /* End of Outputs for SubSystem: '<S48>/Time adjust the gauss coefficients' */

      /* If: '<S54>/If' incorporates:
       *  Constant: '<S48>/Constant'
       *  Sum: '<S48>/Sum1'
       */
      if (rt2s48_iter - 1LL == 0LL) {
        /* Outputs for IfAction SubSystem: '<S54>/If Action Subsystem' incorporates:
         *  ActionPort: '<S60>/Action Port'
         */
        /* Product: '<S60>/Product' incorporates:
         *  Constant: '<S60>/Constant'
         *  Selector: '<S60>/Selector'
         *  Sum: '<S51>/Sum2'
         *  Sum: '<S60>/Sum'
         */
        rt2b_Selector1_f = rt2B_n.rt2Sum2_b[(rt2s89_iter + 1) * 13];

        /* Merge: '<S54>/Merge' incorporates:
         *  Constant: '<S54>/Constant1'
         *  Gain: '<S60>/Gain1'
         *  Product: '<S60>/Product'
         *  Selector: '<S54>/cp[m+1]'
         *  Selector: '<S60>/Selector'
         *  Sum: '<S54>/Sum4'
         */
        rt2rt2Merge_c = rt2B_n.rt2OutportBufferForcp13[rt2s48_iter - 1] *
          rt2b_Selector1_f;

        /* Merge: '<S54>/Merge1' incorporates:
         *  Constant: '<S54>/Constant1'
         *  Gain: '<S60>/Gain2'
         *  Product: '<S60>/Product'
         *  Selector: '<S54>/sp[m+1]'
         *  Sum: '<S54>/Sum4'
         */
        rt2rt2b_Sum1_c_idx_3 = rt2B_n.rt2OutportBufferForsp13[rt2s48_iter - 1] *
          rt2b_Selector1_f;

        /* End of Outputs for SubSystem: '<S54>/If Action Subsystem' */
      } else {
        /* Outputs for IfAction SubSystem: '<S54>/If Action Subsystem1' incorporates:
         *  ActionPort: '<S61>/Action Port'
         */
        /* Product: '<S61>/Product' incorporates:
         *  Constant: '<S62>/Constant'
         *  Selector: '<S61>/Selector'
         *  Sum: '<S51>/Sum2'
         *  Sum: '<S62>/Sum'
         */
        rt2b_Selector1_f = rt2B_n.rt2Sum2_b[((rt2s89_iter + 1) * 13 +
          rt2s48_iter) - 1];

        /* Product: '<S61>/Product1' incorporates:
         *  Constant: '<S54>/Constant1'
         *  Product: '<S61>/Product'
         *  Selector: '<S54>/sp[m+1]'
         *  Selector: '<S61>/Selector1'
         *  Sum: '<S51>/Sum2'
         *  Sum: '<S54>/Sum4'
         *  Sum: '<S63>/Sum'
         */
        rt2b_Sum_n = rt2B_n.rt2OutportBufferForsp13[rt2s48_iter - 1];
        rt2rt2b_Sum1_c_idx_3 = rt2B_n.rt2Sum2_b[((rt2s48_iter - 2) * 13 +
          rt2s89_iter) + 1];

        /* Product: '<S61>/Product' incorporates:
         *  Constant: '<S54>/Constant1'
         *  Product: '<S61>/Product1'
         *  Selector: '<S54>/cp[m+1]'
         *  Sum: '<S54>/Sum4'
         */
        rt2rt2Merge_c_tmp = rt2B_n.rt2OutportBufferForcp13[rt2s48_iter - 1];

        /* Merge: '<S54>/Merge' incorporates:
         *  Product: '<S61>/Product'
         *  Product: '<S61>/Product1'
         *  Selector: '<S54>/cp[m+1]'
         *  Selector: '<S54>/sp[m+1]'
         *  Selector: '<S61>/Selector'
         *  Selector: '<S61>/Selector1'
         *  Sum: '<S61>/Sum'
         */
        rt2rt2Merge_c = rt2b_Selector1_f * rt2rt2Merge_c_tmp +
          rt2rt2b_Sum1_c_idx_3 * rt2b_Sum_n;

        /* Merge: '<S54>/Merge1' incorporates:
         *  Product: '<S61>/Product'
         *  Product: '<S61>/Product1'
         *  Sum: '<S61>/Sum1'
         */
        rt2rt2b_Sum1_c_idx_3 = rt2b_Selector1_f * rt2b_Sum_n -
          rt2rt2b_Sum1_c_idx_3 * rt2rt2Merge_c_tmp;

        /* End of Outputs for SubSystem: '<S54>/If Action Subsystem1' */
      }

      /* End of If: '<S54>/If' */

      /* Outputs for Enabled SubSystem: '<S48>/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations' incorporates:
       *  EnablePort: '<S50>/Enable'
       */
      if (rt2b_LogicalOperator) {
        /* If: '<S50>/if n == m elseif (n==1&m==0) elseif (n>1&m~=n)' incorporates:
         *  Constant: '<S48>/Constant'
         *  Sum: '<S48>/Sum1'
         */
        if (rt2s89_iter + 1 == rt2s48_iter - 1) {
          /* Outputs for IfAction SubSystem: '<S50>/If Action Subsystem' incorporates:
           *  ActionPort: '<S64>/Action Port'
           */
          /* UnitDelay: '<S50>/Unit Delay1' incorporates:
           *  Constant: '<S68>/Constant'
           *  Gain: '<S68>/Gain'
           *  Selector: '<S64>/Selector'
           *  Sum: '<S68>/Sum1'
           *  Sum: '<S68>/Sum2'
           */
          rt2b_Selector1_f = rt2DW_l.rt2UnitDelay1_DSTATE_a[(int32_T)
            ((rt2s48_iter - 2LL) * 13LL + (rt2s89_iter + 1)) - 1];

          /* Merge: '<S50>/Merge1' incorporates:
           *  Product: '<S64>/Product'
           *  Product: '<S64>/Product2'
           *  Selector: '<S64>/Selector'
           *  Selector: '<S64>/Selector1'
           *  Sum: '<S64>/Sum'
           *  UnitDelay: '<S50>/Unit Delay'
           *  UnitDelay: '<S50>/Unit Delay1'
           */
          rt2B_n.rt2Merge1_c = rt2DW_l.rt2UnitDelay_DSTATE_m[(13 * rt2s89_iter +
            rt2s48_iter) - 2] * rt2B_n.rt2sqrt_g + rt2b_Selector1_f *
            rt2B_n.rt2Product4;

          /* Merge: '<S50>/Merge' incorporates:
           *  Product: '<S64>/Product1'
           */
          rt2B_n.rt2Merge_o = rt2b_Selector1_f * rt2B_n.rt2sqrt_g;

          /* End of Outputs for SubSystem: '<S50>/If Action Subsystem' */
        } else if ((rt2s89_iter == 0) && (rt2s48_iter - 1LL == 0LL)) {
          /* Outputs for IfAction SubSystem: '<S50>/If Action Subsystem1' incorporates:
           *  ActionPort: '<S65>/Action Port'
           */
          /* Merge: '<S50>/Merge1' incorporates:
           *  Product: '<S65>/Product'
           *  Product: '<S65>/Product2'
           *  Selector: '<S65>/Selector'
           *  Selector: '<S65>/Selector1'
           *  Sum: '<S65>/Sum'
           *  Sum: '<S71>/Sum'
           *  UnitDelay: '<S50>/Unit Delay'
           *  UnitDelay: '<S50>/Unit Delay1'
           */
          rt2B_n.rt2Merge1_c = rt2B_n.rt2Product4 *
            rt2DW_l.rt2UnitDelay_DSTATE_m[0] - rt2B_n.rt2sqrt_g *
            rt2DW_l.rt2UnitDelay1_DSTATE_a[0];

          /* Merge: '<S50>/Merge' incorporates:
           *  Product: '<S65>/Product3'
           *  Selector: '<S65>/Selector'
           *  UnitDelay: '<S50>/Unit Delay1'
           */
          rt2B_n.rt2Merge_o = rt2B_n.rt2Product4 *
            rt2DW_l.rt2UnitDelay1_DSTATE_a[0];

          /* End of Outputs for SubSystem: '<S50>/If Action Subsystem1' */
        } else if ((rt2s89_iter > 0) && (rt2s89_iter + 1 != rt2s48_iter - 1)) {
          /* Outputs for IfAction SubSystem: '<S50>/If Action Subsystem2' incorporates:
           *  ActionPort: '<S66>/Action Port'
           */
          /* Gain: '<S72>/Gain' */
          rt2i = (int32_T)((rt2s48_iter - 1LL) * 13LL);

          /* Switch: '<S66>/Switch1' incorporates:
           *  Constant: '<S66>/Constant1'
           *  Constant: '<S72>/Constant1'
           *  Constant: '<S76>/Constant1'
           *  RelationalOperator: '<S76>/Relational Operator'
           *  Selector: '<S66>/Selector'
           *  Sum: '<S72>/Sum1'
           *  Sum: '<S72>/Sum2'
           *  Sum: '<S76>/Sum2'
           *  UnitDelay: '<S50>/Unit Delay1'
           */
          if ((rt2s89_iter + 1) - 2LL >= rt2s48_iter - 1LL) {
            rt2b_Sum_n = rt2DW_l.rt2UnitDelay1_DSTATE_a[(int32_T)(((int64_T)
              (rt2s89_iter + 1) + rt2i) - 1LL) - 1];
          } else {
            rt2b_Sum_n = 0.0;
          }

          /* UnitDelay: '<S50>/Unit Delay1' incorporates:
           *  Selector: '<S66>/Selector'
           *  Sum: '<S72>/Sum1'
           */
          rt2b_Selector1_f = rt2DW_l.rt2UnitDelay1_DSTATE_a[rt2s89_iter + rt2i];

          /* Selector: '<S66>/Selector2' incorporates:
           *  Constant: '<S66>/k[13][13]'
           */
          rt2rt2Merge_c_tmp = rt2ConstP_d.rt2pooled8[((rt2s89_iter + 1) * 13 +
            rt2s48_iter) - 1];

          /* Merge: '<S50>/Merge' incorporates:
           *  Constant: '<S66>/k[13][13]'
           *  Product: '<S66>/Product2'
           *  Product: '<S66>/Product3'
           *  Selector: '<S66>/Selector'
           *  Selector: '<S66>/Selector2'
           *  Sum: '<S66>/Sum1'
           *  Switch: '<S66>/Switch1'
           *  UnitDelay: '<S50>/Unit Delay1'
           */
          rt2B_n.rt2Merge_o = rt2b_Selector1_f * rt2B_n.rt2Product4 -
            rt2rt2Merge_c_tmp * rt2b_Sum_n;

          /* Switch: '<S66>/Switch' incorporates:
           *  Constant: '<S66>/Constant'
           *  Constant: '<S75>/Constant1'
           *  RelationalOperator: '<S75>/Relational Operator'
           *  Selector: '<S66>/Selector1'
           *  Sum: '<S75>/Sum2'
           *  UnitDelay: '<S50>/Unit Delay'
           */
          if ((rt2s89_iter + 1) - 2LL >= rt2s48_iter - 1LL) {
            rt2b_Sum_n = rt2DW_l.rt2UnitDelay_DSTATE_m[((rt2s89_iter - 1) * 13 +
              rt2s48_iter) - 1];
          } else {
            rt2b_Sum_n = 0.0;
          }

          /* Merge: '<S50>/Merge1' incorporates:
           *  Product: '<S66>/Product'
           *  Product: '<S66>/Product1'
           *  Product: '<S66>/Product4'
           *  Selector: '<S66>/Selector1'
           *  Sum: '<S66>/Sum'
           *  Switch: '<S66>/Switch'
           *  UnitDelay: '<S50>/Unit Delay'
           */
          rt2B_n.rt2Merge1_c = (rt2DW_l.rt2UnitDelay_DSTATE_m[(13 * rt2s89_iter
            + rt2s48_iter) - 1] * rt2B_n.rt2Product4 - rt2b_Selector1_f *
                                rt2B_n.rt2sqrt_g) - rt2rt2Merge_c_tmp *
            rt2b_Sum_n;

          /* End of Outputs for SubSystem: '<S50>/If Action Subsystem2' */
        }

        /* End of If: '<S50>/if n == m elseif (n==1&m==0) elseif (n>1&m~=n)' */

        /* Assignment: '<S50>/Assignment_snorm' incorporates:
         *  Constant: '<S50>/Constant'
         *  Constant: '<S67>/Constant'
         *  Gain: '<S67>/Gain'
         *  Sum: '<S48>/Sum1'
         *  Sum: '<S50>/Sum'
         *  Sum: '<S67>/Sum1'
         *  Sum: '<S67>/Sum2'
         *  UnitDelay: '<S50>/Unit Delay1'
         */
        if (rt2ForIterator_IterationMarker[3] < 2) {
          rt2ForIterator_IterationMarker[3] = 2U;
          memcpy(&rt2B_n.rt2Assignment_snorm[0],
                 &rt2DW_l.rt2UnitDelay1_DSTATE_a[0], 169U * sizeof(real_T));
        }

        rt2B_n.rt2Assignment_snorm[(int32_T)(((rt2s89_iter + 1) + 13LL *
          (rt2s48_iter - 1LL)) + 1LL) - 1] = rt2B_n.rt2Merge_o;

        /* End of Assignment: '<S50>/Assignment_snorm' */

        /* Assignment: '<S50>/Assignment' incorporates:
         *  UnitDelay: '<S50>/Unit Delay'
         */
        if (rt2ForIterator_IterationMarker[2] < 2) {
          rt2ForIterator_IterationMarker[2] = 2U;
          memcpy(&rt2B_n.rt2Assignment_m[0], &rt2DW_l.rt2UnitDelay_DSTATE_m[0],
                 169U * sizeof(real_T));
        }

        rt2B_n.rt2Assignment_m[(rt2s48_iter + 13 * (rt2s89_iter + 1)) - 1] =
          rt2B_n.rt2Merge1_c;

        /* End of Assignment: '<S50>/Assignment' */

        /* Update for UnitDelay: '<S50>/Unit Delay' incorporates:
         *  Assignment: '<S50>/Assignment'
         */
        memcpy(&rt2DW_l.rt2UnitDelay_DSTATE_m[0], &rt2B_n.rt2Assignment_m[0],
               169U * sizeof(real_T));

        /* Update for UnitDelay: '<S50>/Unit Delay1' */
        memcpy(&rt2DW_l.rt2UnitDelay1_DSTATE_a[0], &rt2B_n.rt2Assignment_snorm[0],
               169U * sizeof(real_T));
      }

      /* End of Outputs for SubSystem: '<S48>/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations' */

      /* Product: '<S49>/par' incorporates:
       *  Constant: '<S48>/Constant'
       *  Constant: '<S53>/Constant'
       *  Gain: '<S53>/Gain'
       *  Selector: '<S49>/snorm[n+m*13]'
       *  Sum: '<S48>/Sum1'
       *  Sum: '<S53>/Sum1'
       */
      rt2b_Selector1_f = rt2B_n.rt2Assignment_snorm[(int32_T)(((rt2s48_iter -
        1LL) * 13LL + (rt2s89_iter + 1)) + 1LL) - 1] * rt2rt2b_sincos_o1_j_idx_0;

      /* Sum: '<S49>/Sum2' incorporates:
       *  Constant: '<S49>/fm'
       *  Product: '<S49>/Product1'
       *  Selector: '<S49>/fm[m]'
       *  UnitDelay: '<S49>/Unit Delay3'
       */
      rt2B_n.rt2Sum2 = rt2ConstP_d.rt2fm_Value[rt2s48_iter - 1] *
        rt2b_Selector1_f * rt2rt2b_Sum1_c_idx_3 + rt2b_x1;

      /* Sum: '<S49>/Sum3' incorporates:
       *  Constant: '<S49>/fn'
       *  Product: '<S49>/Product2'
       *  Selector: '<S49>/fn[m]'
       *  UnitDelay: '<S49>/Unit Delay2'
       */
      rt2B_n.rt2Sum3 = rt2ConstP_d.rt2fn_Value[rt2s89_iter + 1] *
        rt2b_Selector1_f * rt2rt2Merge_c + rt2E4;

      /* Outputs for Enabled SubSystem: '<S49>/Special case - North//South Geographic Pole' incorporates:
       *  EnablePort: '<S52>/Enable'
       */
      /* Outputs for IfAction SubSystem: '<S52>/If Action Subsystem2' incorporates:
       *  ActionPort: '<S57>/Action Port'
       */
      /* If: '<S52>/n ==1' incorporates:
       *  Selector: '<S49>/dp[n][m]'
       *  Selector: '<S57>/Selector2'
       */
      rt2rt2Sum1_tmp = ((rt2s89_iter + 1) * 13 + rt2s48_iter) - 1;

      /* End of Outputs for SubSystem: '<S52>/If Action Subsystem2' */
      /* End of Outputs for SubSystem: '<S49>/Special case - North//South Geographic Pole' */

      /* Sum: '<S49>/Sum1' incorporates:
       *  Assignment: '<S50>/Assignment'
       *  Product: '<S49>/Product'
       *  Selector: '<S49>/dp[n][m]'
       *  UnitDelay: '<S49>/Unit Delay1'
       */
      rt2B_n.rt2Sum1 = rt2E6 - rt2B_n.rt2Assignment_m[rt2rt2Sum1_tmp] *
        rt2rt2Merge_c * rt2rt2b_sincos_o1_j_idx_0;

      /* Outputs for Enabled SubSystem: '<S49>/Special case - North//South Geographic Pole' incorporates:
       *  EnablePort: '<S52>/Enable'
       */
      if ((rt2B_n.rt2sqrt_g == 0.0) && (rt2s48_iter - 1LL == 1LL)) {
        rt2DW_l.rt2SpecialcaseNorthSouthGeographicPole_MODE = true;

        /* If: '<S52>/n ==1' incorporates:
         *  Assignment: '<S57>/Assignment2'
         */
        if (rt2s89_iter == 0) {
          /* Outputs for IfAction SubSystem: '<S52>/If Action Subsystem1' incorporates:
           *  ActionPort: '<S56>/Action Port'
           */
          /* Assignment: '<S56>/Assignment2' incorporates:
           *  Selector: '<S56>/pp[n-1]'
           *  UnitDelay: '<S52>/Unit Delay1'
           */
          if (rt2ForIterator_IterationMarker[0] < 2) {
            rt2ForIterator_IterationMarker[0] = 2U;
            memcpy(&rt2B_n.rt2Assignment2_b[0], &rt2DW_l.rt2UnitDelay1_DSTATE_b
                   [0], 13U * sizeof(real_T));
          }

          rt2B_n.rt2Assignment2_b[1] = rt2DW_l.rt2UnitDelay1_DSTATE_b[0];

          /* End of Assignment: '<S56>/Assignment2' */
          /* End of Outputs for SubSystem: '<S52>/If Action Subsystem1' */
        } else {
          /* Outputs for IfAction SubSystem: '<S52>/If Action Subsystem2' incorporates:
           *  ActionPort: '<S57>/Action Port'
           */
          if (rt2ForIterator_IterationMarker[1] < 2) {
            /* Assignment: '<S57>/Assignment2' */
            rt2ForIterator_IterationMarker[1] = 2U;

            /* Assignment: '<S57>/Assignment2' incorporates:
             *  UnitDelay: '<S52>/Unit Delay1'
             */
            memcpy(&rt2B_n.rt2Assignment2_p[0], &rt2DW_l.rt2UnitDelay1_DSTATE_b
                   [0], 13U * sizeof(real_T));
          }

          /* Assignment: '<S57>/Assignment2' incorporates:
           *  Constant: '<S57>/k[13][13]'
           *  Product: '<S57>/Product1'
           *  Product: '<S57>/Product2'
           *  Selector: '<S57>/Selector2'
           *  Selector: '<S57>/pp[n-2] pp[n-1]'
           *  Sum: '<S57>/Sum1'
           *  UnitDelay: '<S52>/Unit Delay1'
           */
          rt2B_n.rt2Assignment2_p[rt2s89_iter + 1] = rt2B_n.rt2Product4 *
            rt2DW_l.rt2UnitDelay1_DSTATE_b[rt2s89_iter] -
            rt2DW_l.rt2UnitDelay1_DSTATE_b[rt2s89_iter - 1] *
            rt2ConstP_d.rt2pooled8[rt2rt2Sum1_tmp];

          /* End of Outputs for SubSystem: '<S52>/If Action Subsystem2' */
        }

        /* SignalConversion generated from: '<S52>/pp[n]' incorporates:
         *  UnitDelay: '<S52>/Unit Delay1'
         */
        rt2b_TmpSignalConversionAtppnInport1[0] =
          rt2DW_l.rt2UnitDelay1_DSTATE_b[0];
        rt2b_TmpSignalConversionAtppnInport1[1] = rt2B_n.rt2Assignment2_b[1];
        memcpy(&rt2b_TmpSignalConversionAtppnInport1[2],
               &rt2B_n.rt2Assignment2_p[2], 11U * sizeof(real_T));

        /* Product: '<S52>/Product2' incorporates:
         *  Selector: '<S52>/pp[n]'
         */
        rt2B_n.rt2Product2 = rt2b_TmpSignalConversionAtppnInport1[rt2s89_iter +
          1] * rt2rt2b_sincos_o1_j_idx_0 * rt2rt2b_Sum1_c_idx_3;

        /* Update for UnitDelay: '<S52>/Unit Delay1' */
        memcpy(&rt2DW_l.rt2UnitDelay1_DSTATE_b[0],
               &rt2b_TmpSignalConversionAtppnInport1[0], 13U * sizeof(real_T));
      } else if (rt2DW_l.rt2SpecialcaseNorthSouthGeographicPole_MODE) {
        /* Disable for Product: '<S52>/Product2' incorporates:
         *  Outport: '<S52>/bpp'
         */
        rt2B_n.rt2Product2 = 0.0;
        rt2DW_l.rt2SpecialcaseNorthSouthGeographicPole_MODE = false;
      }

      /* End of Outputs for SubSystem: '<S49>/Special case - North//South Geographic Pole' */

      /* Sum: '<S49>/Sum5' incorporates:
       *  Constant: '<S48>/Constant'
       *  Constant: '<S55>/Constant'
       *  Constant: '<S55>/Constant1'
       *  Logic: '<S55>/Logical Operator'
       *  RelationalOperator: '<S55>/Relational Operator'
       *  RelationalOperator: '<S55>/Relational Operator1'
       *  Sum: '<S48>/Sum1'
       *  UnitDelay: '<S49>/Unit Delay4'
       */
      rt2B_n.rt2Sum5 = rt2E7 + rt2B_n.rt2Product2;

      /* Update for UnitDelay: '<S49>/Unit Delay3' */
      rt2b_x1 = rt2B_n.rt2Sum2;

      /* Update for UnitDelay: '<S49>/Unit Delay2' */
      rt2E4 = rt2B_n.rt2Sum3;

      /* Update for UnitDelay: '<S49>/Unit Delay1' */
      rt2E6 = rt2B_n.rt2Sum1;

      /* Update for UnitDelay: '<S49>/Unit Delay4' */
      rt2E7 = rt2B_n.rt2Sum5;
    }

    /* End of Sum: '<S40>/Sum' */
    /* End of Outputs for SubSystem: '<S40>/For Iterator Subsystem' */

    /* Sum: '<S40>/Sum1' incorporates:
     *  UnitDelay: '<S40>/Unit Delay2'
     */
    rt2E4 = rt2rt2b_sincos_o1_idx_1 + rt2B_n.rt2Sum1;
    rt2E6 = rt2rt2b_sincos_o2_idx_1 + rt2B_n.rt2Sum2;
    rt2E7 = rt2rt2UnitDelay2_DSTATE_idx_2 + rt2B_n.rt2Sum3;
    rt2rt2b_Sum1_c_idx_3 = rt2rt2UnitDelay2_DSTATE_idx_3 + rt2B_n.rt2Sum5;

    /* Update for UnitDelay: '<S40>/Unit Delay2' */
    rt2rt2b_sincos_o1_idx_1 = rt2E4;
    rt2rt2b_sincos_o2_idx_1 = rt2E6;
    rt2rt2UnitDelay2_DSTATE_idx_2 = rt2E7;
    rt2rt2UnitDelay2_DSTATE_idx_3 = rt2rt2b_Sum1_c_idx_3;

    /* Update for UnitDelay: '<S40>/Unit Delay' */
    rt2rt2Merge_c = rt2rt2b_sincos_o1_j_idx_0;
  }

  /* End of Outputs for SubSystem: '<S29>/Compute magnetic vector in spherical coordinates' */

  /* Switch: '<S92>/Switch' incorporates:
   *  Product: '<S92>/Product'
   */
  if (rt2B_n.rt2sqrt_g != 0.0) {
    rt2rt2b_sincos_o1_n_idx_1 = rt2E6 / rt2B_n.rt2sqrt_g;
  } else {
    rt2rt2b_sincos_o1_n_idx_1 = rt2rt2b_Sum1_c_idx_3;
  }

  /* End of Switch: '<S92>/Switch' */

  /* Sum: '<S91>/Sum1' incorporates:
   *  Product: '<S91>/Product1'
   *  Product: '<S91>/Product4'
   */
  rt2rt2b_sincos_o1_j_idx_0 = (0.0 - rt2B_n.rt2Product11 * rt2E4) - rt2E7 *
    rt2B_n.rt2Product12;

  /* Sum: '<S93>/Sum1' incorporates:
   *  Product: '<S93>/Product1'
   *  Product: '<S93>/Product4'
   */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  rt2b_Selector1_f = rt2B_n.rt2Product12 * rt2E4 - rt2E7 * rt2B_n.rt2Product11;

  /* Sum: '<S94>/Sum' incorporates:
   *  Product: '<S94>/Product'
   *  Product: '<S94>/Product1'
   */
  rt2rt2Merge_c = rt2rt2b_sincos_o1_n_idx_1 * rt2rt2b_sincos_o1_n_idx_1 +
    rt2rt2b_sincos_o1_j_idx_0 * rt2rt2b_sincos_o1_j_idx_0;

  /* UnitConversion: '<S30>/Unit Conversion' incorporates:
   *  Sqrt: '<S94>/sqrt1'
   *  Trigonometry: '<S94>/Trigonometric Function'
   *  Trigonometry: '<S94>/Trigonometric Function1'
   *  UnitConversion: '<S95>/Unit Conversion'
   *  UnitConversion: '<S96>/Unit Conversion'
   */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  rt2rt2b_sincos_o1_j_idx_0 = 57.295779513082323 * rt_atan2d_snf
    (rt2rt2b_sincos_o1_n_idx_1, rt2rt2b_sincos_o1_j_idx_0) *
    0.017453292519943295;
  rt2rt2b_sincos_o1_n_idx_1 = 57.295779513082323 * rt_atan2d_snf
    (rt2b_Selector1_f, sqrt(rt2rt2Merge_c)) * 0.017453292519943295;

  /* Sqrt: '<S94>/sqrt' incorporates:
   *  Product: '<S94>/Product2'
   *  Sum: '<S94>/Sum1'
   */
  rt2rt2Merge_c = sqrt(rt2b_Selector1_f * rt2b_Selector1_f + rt2rt2Merge_c);

  /* Product: '<S23>/h1' incorporates:
   *  Trigonometry: '<S23>/sincos'
   */
  rt2b_Selector1_f = cos(rt2rt2b_sincos_o1_n_idx_1) * rt2rt2Merge_c;

  /* Product: '<S4>/Matrix Multiply' incorporates:
   *  Product: '<S23>/x1'
   *  Product: '<S23>/y1'
   *  Product: '<S23>/z1'
   *  Trigonometry: '<S23>/sincos'
   */
  rt2b_x1 = cos(rt2rt2b_sincos_o1_j_idx_0) * rt2b_Selector1_f;
  rt2rt2b_sincos_o1_j_idx_0 = sin(rt2rt2b_sincos_o1_j_idx_0) * rt2b_Selector1_f;
  rt2rt2b_sincos_o1_n_idx_1 = sin(rt2rt2b_sincos_o1_n_idx_1) * rt2rt2Merge_c;

  /* Product: '<S4>/Matrix Multiply1' incorporates:
   *  Outport: '<Root>/ECI_B_ECI2SC'
   */
  rt2b_Selector1_f = 0.0;
  rt2E4 = 0.0;
  rt2E6 = 0.0;
  for (rt2s89_iter = 0; rt2s89_iter < 3; rt2s89_iter++) {
    /* Product: '<S4>/Matrix Multiply' incorporates:
     *  Concatenate: '<S19>/Vector Concatenate'
     *  Math: '<S4>/Transpose1'
     *  Product: '<S4>/Matrix Multiply1'
     */
    rt2i = 3 * rt2s89_iter + 1;
    rt2s48_iter = 3 * rt2s89_iter + 2;
    rt2b_MatrixMultiply2[rt2s89_iter] = (rt2b_R_ECEF2ECI[3 * rt2s89_iter] *
      rt2b_x1 + rt2b_R_ECEF2ECI[rt2i] * rt2rt2b_sincos_o1_j_idx_0) +
      rt2b_R_ECEF2ECI[rt2s48_iter] * rt2rt2b_sincos_o1_n_idx_1;

    /* Product: '<S4>/Matrix Multiply1' incorporates:
     *  Math: '<S4>/Transpose'
     *  Outport: '<Root>/ECI_B_ECI2SC'
     */
    rt2b_Sum_n = rt2b_MatrixMultiply2[rt2s89_iter];
    rt2b_Selector1_f += rt2b_Transpose[3 * rt2s89_iter] * rt2b_Sum_n;
    rt2E4 += rt2b_Transpose[rt2i] * rt2b_Sum_n;
    rt2E6 += rt2b_Transpose[rt2s48_iter] * rt2b_Sum_n;
  }

  /* Product: '<S4>/Matrix Multiply1' incorporates:
   *  Outport: '<Root>/ECI_B_ECI2SC'
   */
  rt2Y.rt2ECI_B_ECI2SC[2] = rt2E6;
  rt2Y.rt2ECI_B_ECI2SC[1] = rt2E4;
  rt2Y.rt2ECI_B_ECI2SC[0] = rt2b_Selector1_f;

  /* Delay: '<S1>/Delay3' incorporates:
   *  Inport: '<Root>/v_flag'
   *  Inport: '<Root>/v_pred_init'
   */
  if ((rt_ZCFcn(ANY_ZERO_CROSSING,&rt2PrevZCX_p.rt2Delay3_Reset_ZCE,
                (rt2U.rt2v_flag)) != NO_ZCEVENT) || (rt2U.rt2v_flag != 0.0)) {
    rt2DW_l.rt2icLoad_g = true;
  }

  if (rt2DW_l.rt2icLoad_g) {
    rt2DW_l.rt2Delay3_DSTATE[0] = rt2U.rt2v_pred_init[0];
    rt2DW_l.rt2Delay3_DSTATE[1] = rt2U.rt2v_pred_init[1];
    rt2DW_l.rt2Delay3_DSTATE[2] = rt2U.rt2v_pred_init[2];
  }

  /* MATLAB Function: '<S1>/HPOP_RK4' incorporates:
   *  Inport: '<Root>/JD'
   */
  rt2JD.contents = rt2U.rt2JD;

  /* Outport: '<Root>/v_pred' incorporates:
   *  Delay: '<S1>/Delay3'
   */
  rt2Y.rt2v_pred[0] = rt2DW_l.rt2Delay3_DSTATE[0];

  /* MATLAB Function: '<S1>/HPOP_RK4' incorporates:
   *  Inport: '<Root>/accel_thruster'
   */
  rt2accel_thruster.contents[0] = rt2U.rt2accel_thruster[0];

  /* Outport: '<Root>/v_pred' incorporates:
   *  Delay: '<S1>/Delay3'
   */
  rt2Y.rt2v_pred[1] = rt2DW_l.rt2Delay3_DSTATE[1];

  /* MATLAB Function: '<S1>/HPOP_RK4' incorporates:
   *  Inport: '<Root>/accel_thruster'
   */
  rt2accel_thruster.contents[1] = rt2U.rt2accel_thruster[1];

  /* Outport: '<Root>/v_pred' incorporates:
   *  Delay: '<S1>/Delay3'
   */
  rt2Y.rt2v_pred[2] = rt2DW_l.rt2Delay3_DSTATE[2];

  /* MATLAB Function: '<S1>/HPOP_RK4' incorporates:
   *  Delay: '<S1>/Delay2'
   *  Delay: '<S1>/Delay3'
   *  Inport: '<Root>/accel_thruster'
   *  Inport: '<Root>/dt'
   */
  rt2accel_thruster.contents[2] = rt2U.rt2accel_thruster[2];
  rt2b_atmospheric_drag_terms.contents = 0.027500000000000004;
  rt2b_sun_radiation_pressure_terms.contents = 4.9508333333333345E-8;
  rt2b_const_GM_moon.contents = 4.902801056E+12;
  rt2b_const_GM_earth.contents = 3.986004418E+14;
  rt2b_const_GM_sun.contents = 1.3271741784E+20;
  rt2b_H_base.contents = 500000.0;
  rt2b_H_scale.contents = 63828.0;
  rt2b_rho_o.contents = 6.967E-13;
  rt2accel(&rt2b_const_GM_earth, &rt2b_H_base, &rt2b_H_scale, &rt2b_rho_o,
           &rt2b_atmospheric_drag_terms, &rt2orig2sun, &rt2JD,
           &rt2b_sun_radiation_pressure_terms, &rt2b_const_GM_moon,
           &rt2b_const_GM_sun, &rt2accel_thruster, rt2DW_l.rt2Delay2_DSTATE,
           rt2DW_l.rt2Delay3_DSTATE, rt2b_MatrixMultiply2);
  rt2b_Selector1_f = 0.5 * rt2U.rt2dt;
  rt2a2_tmp[0] = rt2b_Selector1_f * rt2b_MatrixMultiply2[0] +
    rt2DW_l.rt2Delay3_DSTATE[0];
  rt2tmp[0] = rt2b_Selector1_f * rt2DW_l.rt2Delay3_DSTATE[0] +
    rt2DW_l.rt2Delay2_DSTATE[0];
  rt2a2_tmp[1] = rt2b_Selector1_f * rt2b_MatrixMultiply2[1] +
    rt2DW_l.rt2Delay3_DSTATE[1];
  rt2tmp[1] = rt2b_Selector1_f * rt2DW_l.rt2Delay3_DSTATE[1] +
    rt2DW_l.rt2Delay2_DSTATE[1];
  rt2a2_tmp[2] = rt2b_Selector1_f * rt2b_MatrixMultiply2[2] +
    rt2DW_l.rt2Delay3_DSTATE[2];
  rt2tmp[2] = rt2b_Selector1_f * rt2DW_l.rt2Delay3_DSTATE[2] +
    rt2DW_l.rt2Delay2_DSTATE[2];
  rt2accel(&rt2b_const_GM_earth, &rt2b_H_base, &rt2b_H_scale, &rt2b_rho_o,
           &rt2b_atmospheric_drag_terms, &rt2orig2sun, &rt2JD,
           &rt2b_sun_radiation_pressure_terms, &rt2b_const_GM_moon,
           &rt2b_const_GM_sun, &rt2accel_thruster, rt2tmp, rt2a2_tmp, rt2a2);
  rt2a3_tmp[0] = rt2b_Selector1_f * rt2a2[0] + rt2DW_l.rt2Delay3_DSTATE[0];
  rt2tmp[0] = rt2b_Selector1_f * rt2a2_tmp[0] + rt2DW_l.rt2Delay2_DSTATE[0];
  rt2a3_tmp[1] = rt2b_Selector1_f * rt2a2[1] + rt2DW_l.rt2Delay3_DSTATE[1];
  rt2tmp[1] = rt2b_Selector1_f * rt2a2_tmp[1] + rt2DW_l.rt2Delay2_DSTATE[1];
  rt2a3_tmp[2] = rt2b_Selector1_f * rt2a2[2] + rt2DW_l.rt2Delay3_DSTATE[2];
  rt2tmp[2] = rt2b_Selector1_f * rt2a2_tmp[2] + rt2DW_l.rt2Delay2_DSTATE[2];
  rt2accel(&rt2b_const_GM_earth, &rt2b_H_base, &rt2b_H_scale, &rt2b_rho_o,
           &rt2b_atmospheric_drag_terms, &rt2orig2sun, &rt2JD,
           &rt2b_sun_radiation_pressure_terms, &rt2b_const_GM_moon,
           &rt2b_const_GM_sun, &rt2accel_thruster, rt2tmp, rt2a3_tmp, rt2a3);
  rt2a4_tmp[0] = rt2U.rt2dt * rt2a3[0] + rt2DW_l.rt2Delay3_DSTATE[0];
  rt2tmp[0] = rt2U.rt2dt * rt2a3_tmp[0] + rt2DW_l.rt2Delay2_DSTATE[0];
  rt2a4_tmp[1] = rt2U.rt2dt * rt2a3[1] + rt2DW_l.rt2Delay3_DSTATE[1];
  rt2tmp[1] = rt2U.rt2dt * rt2a3_tmp[1] + rt2DW_l.rt2Delay2_DSTATE[1];
  rt2a4_tmp[2] = rt2U.rt2dt * rt2a3[2] + rt2DW_l.rt2Delay3_DSTATE[2];
  rt2tmp[2] = rt2U.rt2dt * rt2a3_tmp[2] + rt2DW_l.rt2Delay2_DSTATE[2];
  rt2accel(&rt2b_const_GM_earth, &rt2b_H_base, &rt2b_H_scale, &rt2b_rho_o,
           &rt2b_atmospheric_drag_terms, &rt2orig2sun, &rt2JD,
           &rt2b_sun_radiation_pressure_terms, &rt2b_const_GM_moon,
           &rt2b_const_GM_sun, &rt2accel_thruster, rt2tmp, rt2a4_tmp, rt2a4);
  rt2b_Selector1_f = rt2U.rt2dt / 6.0;

  /* Outport: '<Root>/orig2sun' incorporates:
   *  MATLAB Function: '<S1>/HPOP_RK4'
   */
  rt2Y.rt2orig2sun[0] = rt2orig2sun.contents[0];

  /* Outport: '<Root>/r_pred' incorporates:
   *  Delay: '<S1>/Delay2'
   */
  rt2Y.rt2r_pred[0] = rt2DW_l.rt2Delay2_DSTATE[0];

  /* Outport: '<Root>/orig2sun' incorporates:
   *  MATLAB Function: '<S1>/HPOP_RK4'
   */
  rt2Y.rt2orig2sun[1] = rt2orig2sun.contents[1];

  /* Outport: '<Root>/r_pred' incorporates:
   *  Delay: '<S1>/Delay2'
   */
  rt2Y.rt2r_pred[1] = rt2DW_l.rt2Delay2_DSTATE[1];

  /* Outport: '<Root>/orig2sun' incorporates:
   *  MATLAB Function: '<S1>/HPOP_RK4'
   */
  rt2Y.rt2orig2sun[2] = rt2orig2sun.contents[2];

  /* Outport: '<Root>/r_pred' incorporates:
   *  Delay: '<S1>/Delay2'
   */
  rt2Y.rt2r_pred[2] = rt2DW_l.rt2Delay2_DSTATE[2];

  /* Outport: '<Root>/LLA_deg_ADCS' */
  rt2Y.rt2LLA_deg_ADCS[0] = rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0;
  rt2Y.rt2LLA_deg_ADCS[1] = rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1;
  rt2Y.rt2LLA_deg_ADCS[2] = rt2b_ECEF_Position_to_LLA1_o2;

  /* Assertion: '<S24>/Assertion' incorporates:
   *  Constant: '<S24>/max_val'
   *  Constant: '<S24>/min_val'
   *  Constant: '<S29>/epoch'
   *  Logic: '<S24>/conjunction'
   *  RelationalOperator: '<S24>/max_relop'
   *  RelationalOperator: '<S24>/min_relop'
   *  Sum: '<S29>/Sum'
   */
  utAssert((rt2b_Abs1 - 2025.0 >= 0.0) && (rt2b_Abs1 - 2025.0 <= 5.0));

  /* Assertion: '<S20>/Assertion' incorporates:
   *  Constant: '<S20>/max_val'
   *  Constant: '<S20>/min_val'
   *  Logic: '<S20>/conjunction'
   *  RelationalOperator: '<S20>/max_relop'
   *  RelationalOperator: '<S20>/min_relop'
   */
  utAssert((rt2b_ECEF_Position_to_LLA1_o2 >= -1000.0) &&
           (rt2b_ECEF_Position_to_LLA1_o2 <= 850000.0));

  /* Assertion: '<S21>/Assertion' incorporates:
   *  Constant: '<S21>/max_val'
   *  Constant: '<S21>/min_val'
   *  Logic: '<S21>/conjunction'
   *  RelationalOperator: '<S21>/max_relop'
   *  RelationalOperator: '<S21>/min_relop'
   */
  utAssert((rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 >= -90.0) &&
           (rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 <= 90.0));

  /* Assertion: '<S22>/Assertion' incorporates:
   *  Constant: '<S22>/max_val'
   *  Constant: '<S22>/min_val'
   *  Logic: '<S22>/conjunction'
   *  RelationalOperator: '<S22>/max_relop'
   *  RelationalOperator: '<S22>/min_relop'
   */
  utAssert((rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 >= -180.0) &&
           (rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 <= 180.0));

  /* MATLAB Function: '<S1>/MATLAB Function3' incorporates:
   *  Inport: '<Root>/JD'
   */
  rt2b_Sum_n = ((rt2U.rt2JD - 2.451545E+6) + 0.00037249999999999995) / 36525.0;
  rt2rt2b_sincos_o1_n_idx_1 = (125.045 - ((rt2U.rt2JD - 2.451545E+6) +
    0.00037249999999999995) * 0.0529921) * 0.017453292519943295;
  rt2rt2b_sincos_o1_j_idx_0 = (250.089 - ((rt2U.rt2JD - 2.451545E+6) +
    0.00037249999999999995) * 0.1059842) * 0.017453292519943295;
  rt2b_x1 = (((rt2U.rt2JD - 2.451545E+6) + 0.00037249999999999995) * 13.0120009
             + 260.008) * 0.017453292519943295;
  rt2E4 = (((rt2U.rt2JD - 2.451545E+6) + 0.00037249999999999995) * 13.3407154 +
           176.625) * 0.017453292519943295;
  rt2E6 = (((rt2U.rt2JD - 2.451545E+6) + 0.00037249999999999995) * 26.4057084 +
           311.589) * 0.017453292519943295;
  rt2E7 = (((rt2U.rt2JD - 2.451545E+6) + 0.00037249999999999995) * 13.064993 +
           134.963) * 0.017453292519943295;
  rt2rt2b_Sum1_c_idx_3 = (15.134 - ((rt2U.rt2JD - 2.451545E+6) +
    0.00037249999999999995) * 0.1589763) * 0.017453292519943295;
  rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 = (((rt2U.rt2JD - 2.451545E+6) +
    0.00037249999999999995) * 12.9590088 + 25.053) * 0.017453292519943295;
  rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 = sin(rt2rt2b_sincos_o1_n_idx_1);
  rt2rt2b_sincos_o1_idx_1 = sin(rt2rt2b_sincos_o1_j_idx_0);
  rt2rt2b_sincos_o2_idx_1 = sin(rt2b_x1);
  rt2rt2UnitDelay2_DSTATE_idx_2 = sin(rt2E4);
  rt2rt2UnitDelay2_DSTATE_idx_3 = sin(rt2E6);
  rt2rt2Merge_c_tmp = sin(rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0);
  rt2rt2Merge_c = 0.0052 * sin(rt2rt2b_Sum1_c_idx_3);
  rt2b_ECEF_Position_to_LLA1_o2 = (((((((((((((((((rt2U.rt2JD - 2.451545E+6) +
    0.00037249999999999995) * 13.17635815 + 38.3213) - ((rt2U.rt2JD -
    2.451545E+6) + 0.00037249999999999995) * ((rt2U.rt2JD - 2.451545E+6) +
    0.00037249999999999995) * 1.4E-12) + 3.561 *
    rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1) + 0.1208 * rt2rt2b_sincos_o1_idx_1)
    - 0.0642 * rt2rt2b_sincos_o2_idx_1) + 0.0158 * rt2rt2UnitDelay2_DSTATE_idx_2)
    + sin((((rt2U.rt2JD - 2.451545E+6) + 0.00037249999999999995) * 0.9856003 +
           357.529) * 0.017453292519943295) * 0.0252) - 0.0066 *
    rt2rt2UnitDelay2_DSTATE_idx_3) - 0.0047 * sin(rt2E7)) - sin((((rt2U.rt2JD -
    2.451545E+6) + 0.00037249999999999995) * 0.3287146 + 276.617) *
    0.017453292519943295) * 0.0046) + sin((((rt2U.rt2JD - 2.451545E+6) +
    0.00037249999999999995) * 1.7484877 + 34.226) * 0.017453292519943295) *
    0.0028) + rt2rt2Merge_c) + sin((((rt2U.rt2JD - 2.451545E+6) +
    0.00037249999999999995) * 0.0036096 + 119.743) * 0.017453292519943295) *
    0.004) + sin((((rt2U.rt2JD - 2.451545E+6) + 0.00037249999999999995) *
                  0.1643573 + 239.961) * 0.017453292519943295) * 0.0019) -
    0.0044 * rt2rt2Merge_c_tmp) * 0.017453292519943295;
  rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 = ((((((((0.0031 * rt2b_Sum_n +
    269.9949) - 3.8787 * rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1) - 0.1204 *
    rt2rt2b_sincos_o1_idx_1) + 0.07 * rt2rt2b_sincos_o2_idx_1) - 0.0172 *
    rt2rt2UnitDelay2_DSTATE_idx_2) + 0.0072 * rt2rt2UnitDelay2_DSTATE_idx_3) -
    rt2rt2Merge_c) + 0.0043 * rt2rt2Merge_c_tmp) * 0.017453292519943295 +
    1.5707963267948966;
  rt2b_x1 = 1.5707963267948966 - (((((((((0.013 * rt2b_Sum_n + 66.5392) + 1.5419
    * cos(rt2rt2b_sincos_o1_n_idx_1)) + 0.0239 * cos(rt2rt2b_sincos_o1_j_idx_0))
    - 0.0278 * cos(rt2b_x1)) + 0.0068 * cos(rt2E4)) - 0.0029 * cos(rt2E6)) +
    0.0009 * cos(rt2E7)) + 0.0008 * cos(rt2rt2b_Sum1_c_idx_3)) - 0.0009 * cos
    (rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0)) * 0.017453292519943295;
  rt2b_Sum_n = sin(rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1);
  rt2rt2b_sincos_o1_n_idx_1 = cos(rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1);
  rt2rt2b_sincos_o1_j_idx_0 = sin(rt2b_x1);
  rt2b_x1 = cos(rt2b_x1);
  rt2E4 = sin(rt2b_ECEF_Position_to_LLA1_o2);
  rt2b_ECEF_Position_to_LLA1_o2 = cos(rt2b_ECEF_Position_to_LLA1_o2);
  rt2b_Transpose[0] = rt2rt2b_sincos_o1_n_idx_1;
  rt2b_Transpose[3] = -rt2b_Sum_n;
  rt2b_Transpose[6] = 0.0;
  rt2b_Transpose[1] = rt2b_Sum_n;
  rt2b_Transpose[4] = rt2rt2b_sincos_o1_n_idx_1;
  rt2b_Transpose[7] = 0.0;
  rt2rt2b_Transpose[1] = 0.0;
  rt2rt2b_Transpose[4] = rt2b_x1;
  rt2rt2b_Transpose[7] = -rt2rt2b_sincos_o1_j_idx_0;
  rt2rt2b_Transpose[2] = 0.0;
  rt2rt2b_Transpose[5] = rt2rt2b_sincos_o1_j_idx_0;
  rt2rt2b_Transpose[8] = rt2b_x1;
  for (rt2s89_iter = 0; rt2s89_iter < 3; rt2s89_iter++) {
    rt2rt2Sum1_tmp = 3 * rt2s89_iter + 2;
    rt2b_Transpose[rt2rt2Sum1_tmp] = rt2rt2b_R_ECEF2ECI_tmp[rt2s89_iter];
    rt2rt2b_Transpose[3 * rt2s89_iter] = rt2c_0[rt2s89_iter];
    rt2b_R_ECEF2ECI[3 * rt2s89_iter] = 0.0;
    rt2b_R_ECEF2ECI[3 * rt2s89_iter + 1] = 0.0;
    rt2b_R_ECEF2ECI[rt2rt2Sum1_tmp] = 0.0;
  }

  for (rt2s89_iter = 0; rt2s89_iter < 3; rt2s89_iter++) {
    rt2rt2b_sincos_o1_n_idx_1 = rt2b_R_ECEF2ECI[3 * rt2s89_iter];
    rt2rt2Sum1_tmp = 3 * rt2s89_iter + 1;
    rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 = rt2b_R_ECEF2ECI[rt2rt2Sum1_tmp];
    rt2rt2b_Product_tmp = 3 * rt2s89_iter + 2;
    rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 = rt2b_R_ECEF2ECI[rt2rt2b_Product_tmp];
    for (rt2s48_iter = 0; rt2s48_iter < 3; rt2s48_iter++) {
      rt2b_Sum_n = rt2rt2b_Transpose[3 * rt2s89_iter + rt2s48_iter];
      rt2rt2b_sincos_o1_n_idx_1 += rt2b_Transpose[3 * rt2s48_iter] * rt2b_Sum_n;
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 += rt2b_Transpose[3 * rt2s48_iter +
        1] * rt2b_Sum_n;
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 += rt2b_Transpose[3 * rt2s48_iter +
        2] * rt2b_Sum_n;
    }

    rt2b_R_ECEF2ECI[rt2rt2b_Product_tmp] =
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1;
    rt2b_R_ECEF2ECI[rt2rt2Sum1_tmp] = rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0;
    rt2b_R_ECEF2ECI[3 * rt2s89_iter] = rt2rt2b_sincos_o1_n_idx_1;
  }

  rt2rt2b_Transpose[0] = rt2b_ECEF_Position_to_LLA1_o2;
  rt2rt2b_Transpose[3] = -rt2E4;
  rt2rt2b_Transpose[6] = 0.0;
  rt2rt2b_Transpose[1] = rt2E4;
  rt2rt2b_Transpose[4] = rt2b_ECEF_Position_to_LLA1_o2;
  rt2rt2b_Transpose[7] = 0.0;
  for (rt2s89_iter = 0; rt2s89_iter < 3; rt2s89_iter++) {
    rt2i = 3 * rt2s89_iter + 2;
    rt2rt2b_Transpose[rt2i] = rt2rt2b_R_ECEF2ECI_tmp[rt2s89_iter];
    rt2rt2b_sincos_o1_n_idx_1 = 0.0;
    rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 = 0.0;
    rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 = 0.0;
    for (rt2s48_iter = 0; rt2s48_iter < 3; rt2s48_iter++) {
      rt2b_Sum_n = rt2rt2b_Transpose[3 * rt2s89_iter + rt2s48_iter];
      rt2rt2b_sincos_o1_n_idx_1 += rt2b_R_ECEF2ECI[3 * rt2s48_iter] * rt2b_Sum_n;
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 += rt2b_R_ECEF2ECI[3 * rt2s48_iter
        + 1] * rt2b_Sum_n;
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 += rt2b_R_ECEF2ECI[3 * rt2s48_iter
        + 2] * rt2b_Sum_n;
    }

    rt2b_Transpose[rt2i] = rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1;
    rt2b_Transpose[3 * rt2s89_iter + 1] = rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0;
    rt2b_Transpose[3 * rt2s89_iter] = rt2rt2b_sincos_o1_n_idx_1;
  }

  for (rt2s89_iter = 0; rt2s89_iter < 3; rt2s89_iter++) {
    rt2rt2b_sincos_o1_n_idx_1 = 0.0;
    rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 = 0.0;
    rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 = 0.0;
    for (rt2s48_iter = 0; rt2s48_iter < 3; rt2s48_iter++) {
      rt2i = 3 * rt2s89_iter + rt2s48_iter;
      rt2b_Sum_n = rt2b_b[rt2i];
      rt2rt2b_sincos_o1_n_idx_1 += rt2b_Transpose[3 * rt2s48_iter] * rt2b_Sum_n;
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 += rt2b_Transpose[3 * rt2s48_iter +
        1] * rt2b_Sum_n;
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 += rt2b_Transpose[3 * rt2s48_iter +
        2] * rt2b_Sum_n;
      rt2rt2b_Transpose[rt2i] = 0.0;
    }

    rt2b_R_ECEF2ECI[3 * rt2s89_iter + 2] =
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1;
    rt2b_R_ECEF2ECI[3 * rt2s89_iter + 1] =
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0;
    rt2b_R_ECEF2ECI[3 * rt2s89_iter] = rt2rt2b_sincos_o1_n_idx_1;
  }

  for (rt2s89_iter = 0; rt2s89_iter < 3; rt2s89_iter++) {
    rt2rt2b_sincos_o1_n_idx_1 = rt2rt2b_Transpose[3 * rt2s89_iter];
    rt2rt2Sum1_tmp = 3 * rt2s89_iter + 1;
    rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 = rt2rt2b_Transpose[rt2rt2Sum1_tmp];
    rt2rt2b_Product_tmp = 3 * rt2s89_iter + 2;
    rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 =
      rt2rt2b_Transpose[rt2rt2b_Product_tmp];
    for (rt2s48_iter = 0; rt2s48_iter < 3; rt2s48_iter++) {
      rt2i = 3 * rt2s89_iter + rt2s48_iter;
      rt2b_Sum_n = rt2c_b[rt2i];
      rt2rt2b_sincos_o1_n_idx_1 += rt2b_R_ECEF2ECI[3 * rt2s48_iter] * rt2b_Sum_n;
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 += rt2b_R_ECEF2ECI[3 * rt2s48_iter
        + 1] * rt2b_Sum_n;
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 += rt2b_R_ECEF2ECI[3 * rt2s48_iter
        + 2] * rt2b_Sum_n;
      rt2b_Transpose[rt2i] = 0.0;
    }

    rt2rt2b_Transpose[rt2rt2b_Product_tmp] =
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1;
    rt2rt2b_Transpose[rt2rt2Sum1_tmp] = rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0;
    rt2rt2b_Transpose[3 * rt2s89_iter] = rt2rt2b_sincos_o1_n_idx_1;
  }

  for (rt2s89_iter = 0; rt2s89_iter < 3; rt2s89_iter++) {
    rt2b_Sum_n = rt2b_Transpose[3 * rt2s89_iter];
    rt2s48_iter = 3 * rt2s89_iter + 1;
    rt2b_ECEF_Position_to_LLA1_o2 = rt2b_Transpose[rt2s48_iter];
    rt2i = 3 * rt2s89_iter + 2;
    rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 = rt2b_Transpose[rt2i];
    for (rt2rt2Sum1_tmp = 0; rt2rt2Sum1_tmp < 3; rt2rt2Sum1_tmp++) {
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1 = rt2d_b[3 * rt2s89_iter +
        rt2rt2Sum1_tmp];
      rt2b_Sum_n += rt2rt2b_Transpose[3 * rt2rt2Sum1_tmp] *
        rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1;
      rt2b_ECEF_Position_to_LLA1_o2 += rt2rt2b_Transpose[3 * rt2rt2Sum1_tmp + 1]
        * rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1;
      rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0 += rt2rt2b_Transpose[3 *
        rt2rt2Sum1_tmp + 2] * rt2rt2b_ECEF_Position_to_LLA1_o1_idx_1;
    }

    rt2b_Transpose[rt2i] = rt2rt2b_ECEF_Position_to_LLA1_o1_idx_0;
    rt2b_Transpose[rt2s48_iter] = rt2b_ECEF_Position_to_LLA1_o2;
    rt2b_Transpose[3 * rt2s89_iter] = rt2b_Sum_n;
  }

  /* Update for Delay: '<S1>/Delay2' */
  rt2DW_l.rt2icLoad = false;

  /* Update for Memory: '<S45>/olon' */
  rt2DW_l.rt2olon_PreviousInput = rt2q;

  /* Update for Memory: '<S46>/otime' */
  rt2DW_l.rt2otime_PreviousInput = rt2b_Abs1;

  /* Update for Memory: '<S44>/olat' */
  rt2DW_l.rt2olat_PreviousInput = rt2c;

  /* Update for Memory: '<S44>/oalt' */
  rt2DW_l.rt2oalt_PreviousInput = rt2rt2Merge;

  /* Update for Delay: '<S1>/Delay3' */
  rt2DW_l.rt2icLoad_g = false;
  for (rt2s89_iter = 0; rt2s89_iter < 3; rt2s89_iter++) {
    /* Outport: '<Root>/R_MCI_moon' incorporates:
     *  MATLAB Function: '<S1>/MATLAB Function3'
     */
    rt2Y.rt2R_MCI_moon[3 * rt2s89_iter] = rt2b_Transpose[rt2s89_iter];
    rt2Y.rt2R_MCI_moon[3 * rt2s89_iter + 1] = rt2b_Transpose[rt2s89_iter + 3];
    rt2Y.rt2R_MCI_moon[3 * rt2s89_iter + 2] = rt2b_Transpose[rt2s89_iter + 6];

    /* Delay: '<S1>/Delay3' incorporates:
     *  Delay: '<S1>/Delay2'
     */
    rt2c = rt2DW_l.rt2Delay3_DSTATE[rt2s89_iter];

    /* Update for Delay: '<S1>/Delay2' incorporates:
     *  Delay: '<S1>/Delay3'
     *  MATLAB Function: '<S1>/HPOP_RK4'
     */
    rt2DW_l.rt2Delay2_DSTATE[rt2s89_iter] += (((2.0 * rt2a2_tmp[rt2s89_iter] +
      rt2c) + 2.0 * rt2a3_tmp[rt2s89_iter]) + rt2a4_tmp[rt2s89_iter]) *
      rt2b_Selector1_f;

    /* Update for Delay: '<S1>/Delay3' incorporates:
     *  MATLAB Function: '<S1>/HPOP_RK4'
     */
    rt2DW_l.rt2Delay3_DSTATE[rt2s89_iter] = (((2.0 * rt2a2[rt2s89_iter] +
      rt2b_MatrixMultiply2[rt2s89_iter]) + 2.0 * rt2a3[rt2s89_iter]) +
      rt2a4[rt2s89_iter]) * rt2b_Selector1_f + rt2c;
  }
}

/* Model initialize function */
void Orbit_Prop_Magnetic_Field_initialize(void)
{
  /* Registration code */

  /* initialize error status */
  rtmSetErrorStatus(rt2M, (NULL));

  /* block I/O */
  (void) memset(((void *) &rt2B_n), 0,
                sizeof(rt2B));

  /* states (dwork) */
  (void) memset((void *)&rt2DW_l, 0,
                sizeof(rt2DW));

  /* external inputs */
  (void)memset(&rt2U, 0, sizeof(rt2ExtU));

  /* external outputs */
  (void)memset(&rt2Y, 0, sizeof(rt2ExtY));

  {
    int32_T rt2i;
    rt2PrevZCX_p.rt2Delay2_Reset_ZCE = UNINITIALIZED_ZCSIG;
    rt2PrevZCX_p.rt2Delay3_Reset_ZCE = UNINITIALIZED_ZCSIG;

    /* InitializeConditions for Delay: '<S1>/Delay2' */
    rt2DW_l.rt2icLoad = true;

    /* InitializeConditions for Memory: '<S45>/olon' */
    rt2DW_l.rt2olon_PreviousInput = -1000.0;

    /* InitializeConditions for Memory: '<S46>/otime' */
    rt2DW_l.rt2otime_PreviousInput = -1000.0;

    /* InitializeConditions for Memory: '<S44>/olat' */
    rt2DW_l.rt2olat_PreviousInput = -1000.0;

    /* InitializeConditions for Memory: '<S44>/oalt' */
    rt2DW_l.rt2oalt_PreviousInput = -1000.0;

    /* InitializeConditions for Delay: '<S1>/Delay3' */
    rt2DW_l.rt2icLoad_g = true;

    /* SystemInitialize for Enabled SubSystem: '<S29>/Convert from geodetic to  spherical coordinates ' */
    /* SystemInitialize for SignalConversion generated from: '<S42>/sp[13]' incorporates:
     *  Outport: '<S42>/sp[13]'
     */
    memcpy(&rt2B_n.rt2OutportBufferForsp13[0], &rt2ConstP_d.rt2sp13_Y0[0], 13U *
           sizeof(real_T));

    /* SystemInitialize for SignalConversion generated from: '<S42>/cp[13]' incorporates:
     *  Outport: '<S42>/cp[13]'
     */
    memcpy(&rt2B_n.rt2OutportBufferForcp13[0], &rt2ConstP_d.rt2cp13_Y0[0], 13U *
           sizeof(real_T));

    /* End of SystemInitialize for SubSystem: '<S29>/Convert from geodetic to  spherical coordinates ' */

    /* SystemInitialize for Enabled SubSystem: '<S29>/Convert from geodetic to  spherical coordinates' */
    /* SystemInitialize for Sqrt: '<S41>/sqrt' incorporates:
     *  Outport: '<S41>/r'
     */
    rt2B_n.rt2sqrt = 6378.137;

    /* SystemInitialize for Product: '<S82>/Product4' incorporates:
     *  Outport: '<S41>/ct'
     */
    rt2B_n.rt2Product4 = 1.0;

    /* End of SystemInitialize for SubSystem: '<S29>/Convert from geodetic to  spherical coordinates' */

    /* SystemInitialize for Iterator SubSystem: '<S29>/Compute magnetic vector in spherical coordinates' */
    /* SystemInitialize for Iterator SubSystem: '<S40>/For Iterator Subsystem' */
    /* SystemInitialize for Enabled SubSystem: '<S48>/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations' */
    /* SystemInitialize for Merge: '<S50>/Merge' */
    rt2B_n.rt2Merge_o = 0.0;

    /* SystemInitialize for Merge: '<S50>/Merge1' */
    rt2B_n.rt2Merge1_c = 0.0;

    /* InitializeConditions for UnitDelay: '<S50>/Unit Delay1' */
    memcpy(&rt2DW_l.rt2UnitDelay1_DSTATE_a[0], &rt2ConstP_d.rt2pooled10[0], 169U
           * sizeof(real_T));

    /* SystemInitialize for Assignment: '<S50>/Assignment_snorm' incorporates:
     *  Outport: '<S50>/snorm[169]'
     */
    memcpy(&rt2B_n.rt2Assignment_snorm[0], &rt2ConstP_d.rt2pooled10[0], 169U *
           sizeof(real_T));

    /* End of SystemInitialize for SubSystem: '<S48>/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations' */

    /* SystemInitialize for Enabled SubSystem: '<S49>/Special case - North//South Geographic Pole' */
    /* SystemInitialize for IfAction SubSystem: '<S52>/If Action Subsystem2' */
    /* SystemInitialize for IfAction SubSystem: '<S52>/If Action Subsystem1' */
    for (rt2i = 0; rt2i < 13; rt2i++) {
      /* InitializeConditions for UnitDelay: '<S52>/Unit Delay1' */
      rt2DW_l.rt2UnitDelay1_DSTATE_b[rt2i] = 1.0;

      /* SystemInitialize for Assignment: '<S56>/Assignment2' incorporates:
       *  Outport: '<S56>/pp[13]'
       */
      rt2B_n.rt2Assignment2_b[rt2i] = 1.0;

      /* SystemInitialize for Assignment: '<S57>/Assignment2' incorporates:
       *  Outport: '<S57>/pp[13]'
       */
      rt2B_n.rt2Assignment2_p[rt2i] = 1.0;
    }

    /* End of SystemInitialize for SubSystem: '<S52>/If Action Subsystem1' */
    /* End of SystemInitialize for SubSystem: '<S52>/If Action Subsystem2' */
    /* End of SystemInitialize for SubSystem: '<S49>/Special case - North//South Geographic Pole' */
    /* End of SystemInitialize for SubSystem: '<S40>/For Iterator Subsystem' */
    /* End of SystemInitialize for SubSystem: '<S29>/Compute magnetic vector in spherical coordinates' */
  }
}

/* Model terminate function */
void Orbit_Prop_Magnetic_Field_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
