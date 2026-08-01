/* cpg_compat.h — Adapta API stacked a per-timestep (generado automaticamente)
 *
 * Dimensiones: NS=7, MI=3, NG=0, T=50
 */
#ifndef CPG_COMPAT_H
#define CPG_COMPAT_H

#include "cpg_solve.h"

/* Parametros horneados como constantes — no-op */
static inline void cpg_update_sqrt_tau(double val)              { (void)val; }
static inline void cpg_update_tau_lamb(double val)              { (void)val; }
static inline void cpg_update_S_x_scaling(int idx, double val)  { (void)idx; (void)val; }
static inline void cpg_update_c_x_scaling(int idx, double val)  { (void)idx; (void)val; }
static inline void cpg_update_S_u_scaling(int idx, double val)  { (void)idx; (void)val; }
static inline void cpg_update_c_u_scaling(int idx, double val)  { (void)idx; (void)val; }


static inline void cpg_update_A_discrete(int idx, double val) {
    int k = idx / 49;
    int inner = idx % 49;
    switch(k) {
        case 0: cpg_update_A_0(inner, val); break;
        case 1: cpg_update_A_1(inner, val); break;
        case 2: cpg_update_A_2(inner, val); break;
        case 3: cpg_update_A_3(inner, val); break;
        case 4: cpg_update_A_4(inner, val); break;
        case 5: cpg_update_A_5(inner, val); break;
        case 6: cpg_update_A_6(inner, val); break;
        case 7: cpg_update_A_7(inner, val); break;
        case 8: cpg_update_A_8(inner, val); break;
        case 9: cpg_update_A_9(inner, val); break;
        case 10: cpg_update_A_10(inner, val); break;
        case 11: cpg_update_A_11(inner, val); break;
        case 12: cpg_update_A_12(inner, val); break;
        case 13: cpg_update_A_13(inner, val); break;
        case 14: cpg_update_A_14(inner, val); break;
        case 15: cpg_update_A_15(inner, val); break;
        case 16: cpg_update_A_16(inner, val); break;
        case 17: cpg_update_A_17(inner, val); break;
        case 18: cpg_update_A_18(inner, val); break;
        case 19: cpg_update_A_19(inner, val); break;
        case 20: cpg_update_A_20(inner, val); break;
        case 21: cpg_update_A_21(inner, val); break;
        case 22: cpg_update_A_22(inner, val); break;
        case 23: cpg_update_A_23(inner, val); break;
        case 24: cpg_update_A_24(inner, val); break;
        case 25: cpg_update_A_25(inner, val); break;
        case 26: cpg_update_A_26(inner, val); break;
        case 27: cpg_update_A_27(inner, val); break;
        case 28: cpg_update_A_28(inner, val); break;
        case 29: cpg_update_A_29(inner, val); break;
        case 30: cpg_update_A_30(inner, val); break;
        case 31: cpg_update_A_31(inner, val); break;
        case 32: cpg_update_A_32(inner, val); break;
        case 33: cpg_update_A_33(inner, val); break;
        case 34: cpg_update_A_34(inner, val); break;
        case 35: cpg_update_A_35(inner, val); break;
        case 36: cpg_update_A_36(inner, val); break;
        case 37: cpg_update_A_37(inner, val); break;
        case 38: cpg_update_A_38(inner, val); break;
        case 39: cpg_update_A_39(inner, val); break;
        case 40: cpg_update_A_40(inner, val); break;
        case 41: cpg_update_A_41(inner, val); break;
        case 42: cpg_update_A_42(inner, val); break;
        case 43: cpg_update_A_43(inner, val); break;
        case 44: cpg_update_A_44(inner, val); break;
        case 45: cpg_update_A_45(inner, val); break;
        case 46: cpg_update_A_46(inner, val); break;
        case 47: cpg_update_A_47(inner, val); break;
        case 48: cpg_update_A_48(inner, val); break;
        case 49: cpg_update_A_49(inner, val); break;
    }
}

static inline void cpg_update_B_discrete(int idx, double val) {
    int k = idx / 21;
    int inner = idx % 21;
    switch(k) {
        case 0: cpg_update_B_0(inner, val); break;
        case 1: cpg_update_B_1(inner, val); break;
        case 2: cpg_update_B_2(inner, val); break;
        case 3: cpg_update_B_3(inner, val); break;
        case 4: cpg_update_B_4(inner, val); break;
        case 5: cpg_update_B_5(inner, val); break;
        case 6: cpg_update_B_6(inner, val); break;
        case 7: cpg_update_B_7(inner, val); break;
        case 8: cpg_update_B_8(inner, val); break;
        case 9: cpg_update_B_9(inner, val); break;
        case 10: cpg_update_B_10(inner, val); break;
        case 11: cpg_update_B_11(inner, val); break;
        case 12: cpg_update_B_12(inner, val); break;
        case 13: cpg_update_B_13(inner, val); break;
        case 14: cpg_update_B_14(inner, val); break;
        case 15: cpg_update_B_15(inner, val); break;
        case 16: cpg_update_B_16(inner, val); break;
        case 17: cpg_update_B_17(inner, val); break;
        case 18: cpg_update_B_18(inner, val); break;
        case 19: cpg_update_B_19(inner, val); break;
        case 20: cpg_update_B_20(inner, val); break;
        case 21: cpg_update_B_21(inner, val); break;
        case 22: cpg_update_B_22(inner, val); break;
        case 23: cpg_update_B_23(inner, val); break;
        case 24: cpg_update_B_24(inner, val); break;
        case 25: cpg_update_B_25(inner, val); break;
        case 26: cpg_update_B_26(inner, val); break;
        case 27: cpg_update_B_27(inner, val); break;
        case 28: cpg_update_B_28(inner, val); break;
        case 29: cpg_update_B_29(inner, val); break;
        case 30: cpg_update_B_30(inner, val); break;
        case 31: cpg_update_B_31(inner, val); break;
        case 32: cpg_update_B_32(inner, val); break;
        case 33: cpg_update_B_33(inner, val); break;
        case 34: cpg_update_B_34(inner, val); break;
        case 35: cpg_update_B_35(inner, val); break;
        case 36: cpg_update_B_36(inner, val); break;
        case 37: cpg_update_B_37(inner, val); break;
        case 38: cpg_update_B_38(inner, val); break;
        case 39: cpg_update_B_39(inner, val); break;
        case 40: cpg_update_B_40(inner, val); break;
        case 41: cpg_update_B_41(inner, val); break;
        case 42: cpg_update_B_42(inner, val); break;
        case 43: cpg_update_B_43(inner, val); break;
        case 44: cpg_update_B_44(inner, val); break;
        case 45: cpg_update_B_45(inner, val); break;
        case 46: cpg_update_B_46(inner, val); break;
        case 47: cpg_update_B_47(inner, val); break;
        case 48: cpg_update_B_48(inner, val); break;
        case 49: cpg_update_B_49(inner, val); break;
    }
}

static inline void cpg_update_y_discrete(int idx, double val) {
    int k = idx / 7;
    int inner = idx % 7;
    switch(k) {
        case 0: cpg_update_y_0(inner, val); break;
        case 1: cpg_update_y_1(inner, val); break;
        case 2: cpg_update_y_2(inner, val); break;
        case 3: cpg_update_y_3(inner, val); break;
        case 4: cpg_update_y_4(inner, val); break;
        case 5: cpg_update_y_5(inner, val); break;
        case 6: cpg_update_y_6(inner, val); break;
        case 7: cpg_update_y_7(inner, val); break;
        case 8: cpg_update_y_8(inner, val); break;
        case 9: cpg_update_y_9(inner, val); break;
        case 10: cpg_update_y_10(inner, val); break;
        case 11: cpg_update_y_11(inner, val); break;
        case 12: cpg_update_y_12(inner, val); break;
        case 13: cpg_update_y_13(inner, val); break;
        case 14: cpg_update_y_14(inner, val); break;
        case 15: cpg_update_y_15(inner, val); break;
        case 16: cpg_update_y_16(inner, val); break;
        case 17: cpg_update_y_17(inner, val); break;
        case 18: cpg_update_y_18(inner, val); break;
        case 19: cpg_update_y_19(inner, val); break;
        case 20: cpg_update_y_20(inner, val); break;
        case 21: cpg_update_y_21(inner, val); break;
        case 22: cpg_update_y_22(inner, val); break;
        case 23: cpg_update_y_23(inner, val); break;
        case 24: cpg_update_y_24(inner, val); break;
        case 25: cpg_update_y_25(inner, val); break;
        case 26: cpg_update_y_26(inner, val); break;
        case 27: cpg_update_y_27(inner, val); break;
        case 28: cpg_update_y_28(inner, val); break;
        case 29: cpg_update_y_29(inner, val); break;
        case 30: cpg_update_y_30(inner, val); break;
        case 31: cpg_update_y_31(inner, val); break;
        case 32: cpg_update_y_32(inner, val); break;
        case 33: cpg_update_y_33(inner, val); break;
        case 34: cpg_update_y_34(inner, val); break;
        case 35: cpg_update_y_35(inner, val); break;
        case 36: cpg_update_y_36(inner, val); break;
        case 37: cpg_update_y_37(inner, val); break;
        case 38: cpg_update_y_38(inner, val); break;
        case 39: cpg_update_y_39(inner, val); break;
        case 40: cpg_update_y_40(inner, val); break;
        case 41: cpg_update_y_41(inner, val); break;
        case 42: cpg_update_y_42(inner, val); break;
        case 43: cpg_update_y_43(inner, val); break;
        case 44: cpg_update_y_44(inner, val); break;
        case 45: cpg_update_y_45(inner, val); break;
        case 46: cpg_update_y_46(inner, val); break;
        case 47: cpg_update_y_47(inner, val); break;
        case 48: cpg_update_y_48(inner, val); break;
        case 49: cpg_update_y_49(inner, val); break;
    }
}

static inline void cpg_update_ox_cvxpy(int idx, double val) {
    int k = idx / 7;
    int inner = idx % 7;
    switch(k) {
        case 0: cpg_update_ox_0(inner, val); break;
        case 1: cpg_update_ox_1(inner, val); break;
        case 2: cpg_update_ox_2(inner, val); break;
        case 3: cpg_update_ox_3(inner, val); break;
        case 4: cpg_update_ox_4(inner, val); break;
        case 5: cpg_update_ox_5(inner, val); break;
        case 6: cpg_update_ox_6(inner, val); break;
        case 7: cpg_update_ox_7(inner, val); break;
        case 8: cpg_update_ox_8(inner, val); break;
        case 9: cpg_update_ox_9(inner, val); break;
        case 10: cpg_update_ox_10(inner, val); break;
        case 11: cpg_update_ox_11(inner, val); break;
        case 12: cpg_update_ox_12(inner, val); break;
        case 13: cpg_update_ox_13(inner, val); break;
        case 14: cpg_update_ox_14(inner, val); break;
        case 15: cpg_update_ox_15(inner, val); break;
        case 16: cpg_update_ox_16(inner, val); break;
        case 17: cpg_update_ox_17(inner, val); break;
        case 18: cpg_update_ox_18(inner, val); break;
        case 19: cpg_update_ox_19(inner, val); break;
        case 20: cpg_update_ox_20(inner, val); break;
        case 21: cpg_update_ox_21(inner, val); break;
        case 22: cpg_update_ox_22(inner, val); break;
        case 23: cpg_update_ox_23(inner, val); break;
        case 24: cpg_update_ox_24(inner, val); break;
        case 25: cpg_update_ox_25(inner, val); break;
        case 26: cpg_update_ox_26(inner, val); break;
        case 27: cpg_update_ox_27(inner, val); break;
        case 28: cpg_update_ox_28(inner, val); break;
        case 29: cpg_update_ox_29(inner, val); break;
        case 30: cpg_update_ox_30(inner, val); break;
        case 31: cpg_update_ox_31(inner, val); break;
        case 32: cpg_update_ox_32(inner, val); break;
        case 33: cpg_update_ox_33(inner, val); break;
        case 34: cpg_update_ox_34(inner, val); break;
        case 35: cpg_update_ox_35(inner, val); break;
        case 36: cpg_update_ox_36(inner, val); break;
        case 37: cpg_update_ox_37(inner, val); break;
        case 38: cpg_update_ox_38(inner, val); break;
        case 39: cpg_update_ox_39(inner, val); break;
        case 40: cpg_update_ox_40(inner, val); break;
        case 41: cpg_update_ox_41(inner, val); break;
        case 42: cpg_update_ox_42(inner, val); break;
        case 43: cpg_update_ox_43(inner, val); break;
        case 44: cpg_update_ox_44(inner, val); break;
        case 45: cpg_update_ox_45(inner, val); break;
        case 46: cpg_update_ox_46(inner, val); break;
        case 47: cpg_update_ox_47(inner, val); break;
        case 48: cpg_update_ox_48(inner, val); break;
        case 49: cpg_update_ox_49(inner, val); break;
        case 50: cpg_update_ox_50(inner, val); break;
    }
}

static inline void cpg_update_ou_cvxpy(int idx, double val) {
    int k = idx / 3;
    int inner = idx % 3;
    switch(k) {
        case 0: cpg_update_ou_0(inner, val); break;
        case 1: cpg_update_ou_1(inner, val); break;
        case 2: cpg_update_ou_2(inner, val); break;
        case 3: cpg_update_ou_3(inner, val); break;
        case 4: cpg_update_ou_4(inner, val); break;
        case 5: cpg_update_ou_5(inner, val); break;
        case 6: cpg_update_ou_6(inner, val); break;
        case 7: cpg_update_ou_7(inner, val); break;
        case 8: cpg_update_ou_8(inner, val); break;
        case 9: cpg_update_ou_9(inner, val); break;
        case 10: cpg_update_ou_10(inner, val); break;
        case 11: cpg_update_ou_11(inner, val); break;
        case 12: cpg_update_ou_12(inner, val); break;
        case 13: cpg_update_ou_13(inner, val); break;
        case 14: cpg_update_ou_14(inner, val); break;
        case 15: cpg_update_ou_15(inner, val); break;
        case 16: cpg_update_ou_16(inner, val); break;
        case 17: cpg_update_ou_17(inner, val); break;
        case 18: cpg_update_ou_18(inner, val); break;
        case 19: cpg_update_ou_19(inner, val); break;
        case 20: cpg_update_ou_20(inner, val); break;
        case 21: cpg_update_ou_21(inner, val); break;
        case 22: cpg_update_ou_22(inner, val); break;
        case 23: cpg_update_ou_23(inner, val); break;
        case 24: cpg_update_ou_24(inner, val); break;
        case 25: cpg_update_ou_25(inner, val); break;
        case 26: cpg_update_ou_26(inner, val); break;
        case 27: cpg_update_ou_27(inner, val); break;
        case 28: cpg_update_ou_28(inner, val); break;
        case 29: cpg_update_ou_29(inner, val); break;
        case 30: cpg_update_ou_30(inner, val); break;
        case 31: cpg_update_ou_31(inner, val); break;
        case 32: cpg_update_ou_32(inner, val); break;
        case 33: cpg_update_ou_33(inner, val); break;
        case 34: cpg_update_ou_34(inner, val); break;
        case 35: cpg_update_ou_35(inner, val); break;
        case 36: cpg_update_ou_36(inner, val); break;
        case 37: cpg_update_ou_37(inner, val); break;
        case 38: cpg_update_ou_38(inner, val); break;
        case 39: cpg_update_ou_39(inner, val); break;
        case 40: cpg_update_ou_40(inner, val); break;
        case 41: cpg_update_ou_41(inner, val); break;
        case 42: cpg_update_ou_42(inner, val); break;
        case 43: cpg_update_ou_43(inner, val); break;
        case 44: cpg_update_ou_44(inner, val); break;
        case 45: cpg_update_ou_45(inner, val); break;
        case 46: cpg_update_ou_46(inner, val); break;
        case 47: cpg_update_ou_47(inner, val); break;
        case 48: cpg_update_ou_48(inner, val); break;
        case 49: cpg_update_ou_49(inner, val); break;
    }
}

#endif /* CPG_COMPAT_H */
