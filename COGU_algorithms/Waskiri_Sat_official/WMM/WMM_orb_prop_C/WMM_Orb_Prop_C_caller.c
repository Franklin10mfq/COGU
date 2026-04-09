#include "Orbit_Prop_Magnetic_Field.h"
#include <stdio.h>
#include <time.h>
void WMM_ECI_B(real_T rt2v_pred_init[3], real_T rt2r_pred_init[3], real_T rt2JD, real_T rt2dt, int T, real_T WMM_ECI_B_pred_T[3][600])
{
  Orbit_Prop_Magnetic_Field_initialize();
  rt2U.rt2v_flag=0;
  rt2U.rt2v_pred_init[0]=rt2v_pred_init[0];
  rt2U.rt2v_pred_init[1]=rt2v_pred_init[1];
  rt2U.rt2v_pred_init[2]=rt2v_pred_init[2];
  rt2U.rt2r_flag=0;
  printf("vz = %.2f\n",rt2U.rt2v_pred_init[2]);
  rt2U.rt2r_pred_init[0]=rt2r_pred_init[0];
  rt2U.rt2r_pred_init[1]=rt2r_pred_init[1];
  rt2U.rt2r_pred_init[2]=rt2r_pred_init[2];
  rt2U.rt2accel_thruster[0]=0;
  rt2U.rt2accel_thruster[1]=0;
  rt2U.rt2accel_thruster[2]=0;
  rt2U.rt2JD=rt2JD;
  rt2U.rt2dt=rt2dt;
  
  for(int i=0;i<=T;i++)
  {
    Orbit_Prop_Magnetic_Field_step();
    WMM_ECI_B_pred_T[0][i]=rt2Y.rt2ECI_B_ECI2SC[0];
    WMM_ECI_B_pred_T[1][i]=rt2Y.rt2ECI_B_ECI2SC[1];
    WMM_ECI_B_pred_T[2][i]=rt2Y.rt2ECI_B_ECI2SC[2];
  }
}