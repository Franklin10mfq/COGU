/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: rt_nonfinite.h
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

#ifndef rt_nonfinite_h_
#define rt_nonfinite_h_
#include "rtwtypes.h"

extern real_T rtInf;
extern real_T rtMinusInf;
extern real_T rtNaN;
extern real32_T rtInfF;
extern real32_T rtMinusInfF;
extern real32_T rtNaNF;
extern boolean_T rtIsInf(real_T value);
extern boolean_T rtIsInfF(real32_T value);
extern boolean_T rtIsNaN(real_T value);
extern boolean_T rtIsNaNF(real32_T value);

#endif                                 /* rt_nonfinite_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
