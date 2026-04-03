/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: Orbit_Prop_Magnetic_Field.h
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

#ifndef Orbit_Prop_Magnetic_Field_h_
#define Orbit_Prop_Magnetic_Field_h_
#ifndef Orbit_Prop_Magnetic_Field_COMMON_INCLUDES_
#define Orbit_Prop_Magnetic_Field_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rt_nonfinite.h"
#include "math.h"
#endif                          /* Orbit_Prop_Magnetic_Field_COMMON_INCLUDES_ */

#include "Orbit_Prop_Magnetic_Field_types.h"
#include "rt_zcfcn.h"
#include "rtGetInf.h"
#include "rtGetNaN.h"
#include <stddef.h>
#include <string.h>
#include "zero_crossing_types.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

/* Block signals (default storage) */
typedef struct {
  real_T rt2OutportBufferForcp13[13];
  real_T rt2OutportBufferForsp13[13];
  real_T rt2sqrt;                      /* '<S41>/sqrt' */
  real_T rt2Product11;                 /* '<S81>/Product11' */
  real_T rt2Product12;                 /* '<S87>/Product12' */
  real_T rt2Product4;                  /* '<S82>/Product4' */
  real_T rt2sqrt_g;                    /* '<S88>/sqrt' */
  real_T rt2Sum2;                      /* '<S49>/Sum2' */
  real_T rt2Sum3;                      /* '<S49>/Sum3' */
  real_T rt2Sum1;                      /* '<S49>/Sum1' */
  real_T rt2Sum5;                      /* '<S49>/Sum5' */
  real_T rt2Assignment[169];           /* '<S51>/Assignment' */
  real_T rt2Sum2_b[169];               /* '<S51>/Sum2' */
  real_T rt2Assignment2[169];          /* '<S79>/Assignment2' */
  real_T rt2Merge_o;                   /* '<S50>/Merge' */
  real_T rt2Assignment_snorm[169];     /* '<S50>/Assignment_snorm' */
  real_T rt2Merge1_c;                  /* '<S50>/Merge1' */
  real_T rt2Assignment_m[169];         /* '<S50>/Assignment' */
  real_T rt2Product2;                  /* '<S52>/Product2' */
  real_T rt2Assignment2_p[13];         /* '<S57>/Assignment2' */
  real_T rt2Assignment2_b[13];         /* '<S56>/Assignment2' */
} rt2B;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T rt2Delay2_DSTATE[3];          /* '<S1>/Delay2' */
  real_T rt2Delay3_DSTATE[3];          /* '<S1>/Delay3' */
  real_T rt2UnitDelay1_DSTATE[2];      /* '<S89>/Unit Delay1' */
  real_T rt2UnitDelay_DSTATE_h[169];   /* '<S51>/Unit Delay' */
  real_T rt2UnitDelay_DSTATE_hg[169];  /* '<S78>/Unit Delay' */
  real_T rt2UnitDelay_DSTATE_m[169];   /* '<S50>/Unit Delay' */
  real_T rt2UnitDelay1_DSTATE_a[169];  /* '<S50>/Unit Delay1' */
  real_T rt2UnitDelay1_DSTATE_b[13];   /* '<S52>/Unit Delay1' */
  real_T rt2olon_PreviousInput;        /* '<S45>/olon' */
  real_T rt2otime_PreviousInput;       /* '<S46>/otime' */
  real_T rt2olat_PreviousInput;        /* '<S44>/olat' */
  real_T rt2oalt_PreviousInput;        /* '<S44>/oalt' */
  boolean_T rt2icLoad;                 /* '<S1>/Delay2' */
  boolean_T rt2icLoad_g;               /* '<S1>/Delay3' */
  boolean_T rt2SpecialcaseNorthSouthGeographicPole_MODE;
                       /* '<S49>/Special case - North//South Geographic Pole' */
} rt2DW;

/* Zero-crossing (trigger) state */
typedef struct {
  ZCSigState rt2Delay2_Reset_ZCE;      /* '<S1>/Delay2' */
  ZCSigState rt2Delay3_Reset_ZCE;      /* '<S1>/Delay3' */
} rt2PrevZCX;

/* Invariant block signals (default storage) */
typedef struct {
  const real_T rt2a2;                  /* '<S41>/a2' */
  const real_T rt2b2;                  /* '<S41>/b2' */
  const real_T rt2c2;                  /* '<S84>/Sum1' */
  const real_T rt2a4;                  /* '<S86>/a4' */
  const real_T rt2b4;                  /* '<S86>/b4' */
  const real_T rt2c4;                  /* '<S86>/Sum9' */
  const real_T rt2c2_d;                /* '<S87>/Sum1' */
} rt2ConstB_f;

/* Constant parameters (default storage) */
typedef struct {
  /* Pooled Parameter (Expression: k)
   * Referenced by:
   *   '<S66>/k[13][13]'
   *   '<S57>/k[13][13]'
   */
  real_T rt2pooled8[169];

  /* Pooled Parameter (Expression: snorm)
   * Referenced by:
   *   '<S50>/snorm[169]'
   *   '<S50>/Unit Delay1'
   */
  real_T rt2pooled10[169];

  /* Expression: c
   * Referenced by: '<S51>/c[maxdef][maxdef]'
   */
  real_T rt2cmaxdefmaxdef_Value[169];

  /* Expression: cd
   * Referenced by: '<S51>/cd[maxdef][maxdef]'
   */
  real_T rt2cdmaxdefmaxdef_Value[169];

  /* Expression: fm
   * Referenced by: '<S49>/fm'
   */
  real_T rt2fm_Value[13];

  /* Expression: fn
   * Referenced by: '<S49>/fn'
   */
  real_T rt2fn_Value[13];

  /* Pooled Parameter (Mixed Expressions)
   * Referenced by:
   *   '<S89>/sp[11]'
   *   '<S89>/cp[11]'
   *   '<S89>/Constant'
   *   '<S89>/Constant1'
   */
  real_T rt2pooled13[11];

  /* Expression: [0 0 (1:(maxdef-1))]
   * Referenced by: '<S42>/sp[13]'
   */
  real_T rt2sp13_Y0[13];

  /* Expression: [1 1 (1:(maxdef-1))]
   * Referenced by: '<S42>/cp[13]'
   */
  real_T rt2cp13_Y0[13];
} rt2ConstP;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T rt2v_flag;                    /* '<Root>/v_flag' */
  real_T rt2v_pred_init[3];            /* '<Root>/v_pred_init' */
  real_T rt2r_flag;                    /* '<Root>/r_flag' */
  real_T rt2r_pred_init[3];            /* '<Root>/r_pred_init' */
  real_T rt2accel_thruster[3];         /* '<Root>/accel_thruster' */
  real_T rt2JD;                        /* '<Root>/JD' */
  real_T rt2dt;                        /* '<Root>/dt' */
} rt2ExtU;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T rt2v_pred[3];                 /* '<Root>/v_pred' */
  real_T rt2r_pred[3];                 /* '<Root>/r_pred' */
  real_T rt2orig2sun[3];               /* '<Root>/orig2sun' */
  real_T rt2R_ECEF2ECI[9];             /* '<Root>/R_ECEF2ECI' */
  real_T rt2ECI_B_ECI2SC[3];           /* '<Root>/ECI_B_ECI2SC' */
  real_T rt2R_MCI_moon[9];             /* '<Root>/R_MCI_moon' */
  real_T rt2LLA_deg_ADCS[3];           /* '<Root>/LLA_deg_ADCS' */
} rt2ExtY;

/* Real-time Model Data Structure */
struct rt2tag_RTM {
  const char_T * volatile errorStatus;
};

/* Block signals (default storage) */
extern rt2B rt2B_n;

/* Block states (default storage) */
extern rt2DW rt2DW_l;

/* Zero-crossing (trigger) state */
extern rt2PrevZCX rt2PrevZCX_p;

/* External inputs (root inport signals with default storage) */
extern rt2ExtU rt2U;

/* External outputs (root outports fed by signals with default storage) */
extern rt2ExtY rt2Y;
extern const rt2ConstB_f rt2ConstB;    /* constant block i/o */

/* Constant parameters (default storage) */
extern const rt2ConstP rt2ConstP_d;

/* Model entry point functions */
extern void Orbit_Prop_Magnetic_Field_initialize(void);
extern void Orbit_Prop_Magnetic_Field_step(void);
extern void Orbit_Prop_Magnetic_Field_terminate(void);

/* Real-time Model object */
extern rt2RT_MODEL *const rt2M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S4>/Scope' : Unused code path elimination
 * Block '<S4>/Scope1' : Unused code path elimination
 * Block '<S4>/Scope2' : Unused code path elimination
 * Block '<S4>/Scope3' : Unused code path elimination
 * Block '<S27>/Unit Conversion' : Unused code path elimination
 * Block '<S28>/Unit Conversion' : Unused code path elimination
 * Block '<S1>/Data Type Conversion2' : Eliminate redundant data type conversion
 * Block '<S19>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S4>/Reshape1' : Reshape block reduction
 * Block '<S24>/maxtype' : Eliminate redundant data type conversion
 * Block '<S24>/mintype' : Eliminate redundant data type conversion
 * Block '<S26>/Unit Conversion' : Eliminated nontunable gain of 1
 * Block '<S8>/Unit Conversion2' : Eliminated nontunable gain of 1
 * Block '<S57>/Reshape' : Reshape block reduction
 * Block '<S64>/Reshape' : Reshape block reduction
 * Block '<S65>/Reshape' : Reshape block reduction
 * Block '<S66>/Reshape' : Reshape block reduction
 * Block '<S66>/Reshape1' : Reshape block reduction
 * Block '<S42>/Gain' : Eliminated nontunable gain of 1
 * Block '<S42>/Gain1' : Eliminated nontunable gain of 1
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Note that this particular code originates from a subsystem build,
 * and has its own system numbers different from the parent model.
 * Refer to the system hierarchy for this subsystem below, and use the
 * MATLAB hilite_system command to trace the generated code back
 * to the parent model.  For example,
 *
 * hilite_system('WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field')    - opens subsystem WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field
 * hilite_system('WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Kp') - opens and selects block Kp
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'WaskiriSat_COGU_orb_prop_WMM'
 * '<S1>'   : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field'
 * '<S2>'   : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/HPOP_RK4'
 * '<S3>'   : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/MATLAB Function3'
 * '<S4>'   : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model'
 * '<S5>'   : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Subsystem2'
 * '<S6>'   : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/Direction Cosine Matrix ECEF to NED'
 * '<S7>'   : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/MATLAB Function1'
 * '<S8>'   : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model'
 * '<S9>'   : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/Direction Cosine Matrix ECEF to NED/A11'
 * '<S10>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/Direction Cosine Matrix ECEF to NED/A12'
 * '<S11>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/Direction Cosine Matrix ECEF to NED/A13'
 * '<S12>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/Direction Cosine Matrix ECEF to NED/A21'
 * '<S13>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/Direction Cosine Matrix ECEF to NED/A22'
 * '<S14>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/Direction Cosine Matrix ECEF to NED/A23'
 * '<S15>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/Direction Cosine Matrix ECEF to NED/A31'
 * '<S16>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/Direction Cosine Matrix ECEF to NED/A32'
 * '<S17>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/Direction Cosine Matrix ECEF to NED/A33'
 * '<S18>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/Direction Cosine Matrix ECEF to NED/Angle Conversion'
 * '<S19>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/Direction Cosine Matrix ECEF to NED/Create Transformation Matrix'
 * '<S20>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/Check Altitude'
 * '<S21>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/Check Latitude'
 * '<S22>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/Check Longitude'
 * '<S23>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/Compute x,y,z, and h components of magnetic field'
 * '<S24>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/Is time within model limits'
 * '<S25>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/LatLong wrap'
 * '<S26>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/Length Conversion'
 * '<S27>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/MagField Conversion'
 * '<S28>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/MagField Conversion1'
 * '<S29>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag'
 * '<S30>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/Compute x,y,z, and h components of magnetic field/Angle Conversion'
 * '<S31>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/LatLong wrap/Latitude Wrap 90'
 * '<S32>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/LatLong wrap/Wrap Longitude'
 * '<S33>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/LatLong wrap/Latitude Wrap 90/Compare To Constant'
 * '<S34>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/LatLong wrap/Latitude Wrap 90/Sign1 Or Zero'
 * '<S35>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/LatLong wrap/Latitude Wrap 90/Wrap Angle 180'
 * '<S36>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/LatLong wrap/Latitude Wrap 90/Sign1 Or Zero/If Action Subsystem'
 * '<S37>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/LatLong wrap/Latitude Wrap 90/Sign1 Or Zero/If Action Subsystem1'
 * '<S38>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/LatLong wrap/Latitude Wrap 90/Wrap Angle 180/Compare To Constant'
 * '<S39>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/LatLong wrap/Wrap Longitude/Compare To Constant'
 * '<S40>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates'
 * '<S41>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Convert from geodetic to  spherical coordinates'
 * '<S42>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Convert from geodetic to  spherical coordinates '
 * '<S43>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Get Cosine and Sine  of Latitude and Longitude'
 * '<S44>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Has altitude or latitude changed'
 * '<S45>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Has longitude changed '
 * '<S46>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Has time changed'
 * '<S47>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity'
 * '<S48>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem'
 * '<S49>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion'
 * '<S50>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations'
 * '<S51>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Time adjust the gauss coefficients'
 * '<S52>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/Special case - North//South Geographic Pole'
 * '<S53>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/calculate  index'
 * '<S54>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/calculate temp values'
 * '<S55>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/special case'
 * '<S56>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/Special case - North//South Geographic Pole/If Action Subsystem1'
 * '<S57>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/Special case - North//South Geographic Pole/If Action Subsystem2'
 * '<S58>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/Special case - North//South Geographic Pole/If Action Subsystem2/calculate  indices'
 * '<S59>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/Special case - North//South Geographic Pole/If Action Subsystem2/calculate  row and column'
 * '<S60>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/calculate temp values/If Action Subsystem'
 * '<S61>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/calculate temp values/If Action Subsystem1'
 * '<S62>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/calculate temp values/If Action Subsystem1/m,n'
 * '<S63>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Accumulate terms of the  spherical harmonic expansion/calculate temp values/If Action Subsystem1/n,m-1'
 * '<S64>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem'
 * '<S65>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem1'
 * '<S66>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem2'
 * '<S67>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/calculate  index'
 * '<S68>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem/calculate  index'
 * '<S69>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem/calculate  row and column'
 * '<S70>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem1/calculate  index'
 * '<S71>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem1/calculate  row and column'
 * '<S72>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem2/calculate  indices'
 * '<S73>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem2/calculate  row and column1'
 * '<S74>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem2/calculate  row and column2'
 * '<S75>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem2/m<n-2'
 * '<S76>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Compute unnormalized associated  legendre polynomials and  derivatives via recursion relations/If Action Subsystem2/m<n-2 '
 * '<S77>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Time adjust the gauss coefficients/If Action Subsystem'
 * '<S78>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Time adjust the gauss coefficients/if (m~=0)'
 * '<S79>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Time adjust the gauss coefficients/if (m~=0)/If Action Subsystem1'
 * '<S80>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Compute magnetic vector in spherical coordinates/For Iterator Subsystem/Time adjust the gauss coefficients/if (m~=0)/If Action Subsystem2'
 * '<S81>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Convert from geodetic to  spherical coordinates/calculate ca'
 * '<S82>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Convert from geodetic to  spherical coordinates/calculate ct'
 * '<S83>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Convert from geodetic to  spherical coordinates/calculate d'
 * '<S84>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Convert from geodetic to  spherical coordinates/calculate q'
 * '<S85>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Convert from geodetic to  spherical coordinates/calculate q2'
 * '<S86>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Convert from geodetic to  spherical coordinates/calculate r2'
 * '<S87>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Convert from geodetic to  spherical coordinates/calculate sa'
 * '<S88>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Convert from geodetic to  spherical coordinates/calculate st'
 * '<S89>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Convert from geodetic to  spherical coordinates /For Iterator Subsystem'
 * '<S90>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Get Cosine and Sine  of Latitude and Longitude/Angle Conversion2'
 * '<S91>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity/Calculate bx'
 * '<S92>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity/Calculate by'
 * '<S93>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity/Calculate bz'
 * '<S94>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity/Compute declination, inclination,  and total intensity'
 * '<S95>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity/Compute declination, inclination,  and total intensity/Angle Conversion'
 * '<S96>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Magnetic_Field_Model/World Magnetic Model/geomag/Rotate magnetic vector components  to geodetic from spherical and  compute declination, inclination  and total intensity/Compute declination, inclination,  and total intensity/Angle Conversion1'
 * '<S97>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Subsystem2/MATLAB Function'
 * '<S98>'  : 'WaskiriSat_COGU_orb_prop_WMM/Orbit_Prop_Magnetic_Field/Subsystem2/MATLAB Function8'
 */
#endif                                 /* Orbit_Prop_Magnetic_Field_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
