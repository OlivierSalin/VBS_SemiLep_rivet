#include "EventWeights.h"

double eventWeight;
double eventWeight_fm0_quad, eventWeight_fm1_quad, eventWeight_fm2_quad, eventWeight_fm3_quad, eventWeight_fm4_quad;
double eventWeight_fm5_quad, eventWeight_fm7_quad, eventWeight_fm8_quad, eventWeight_fm9_quad;
double eventWeight_fs0_quad, eventWeight_fs1_quad, eventWeight_fs2_quad;
double eventWeight_ft0_quad, eventWeight_ft1_quad, eventWeight_ft2_quad, eventWeight_ft3_quad, eventWeight_ft4_quad, eventWeight_ft5_quad, eventWeight_ft6_quad;

double eventWeight_fs0_fs1_cross, eventWeight_fs0_fs2_cross;
double eventWeight_fs1_fs2_cross;

double eventWeight_ft0_ft1_cross, eventWeight_ft0_ft2_cross, eventWeight_ft0_ft3_cross, eventWeight_ft0_ft4_cross, eventWeight_ft0_ft5_cross, eventWeight_ft0_ft6_cross;
double eventWeight_ft1_ft2_cross, eventWeight_ft1_ft3_cross, eventWeight_ft1_ft4_cross, eventWeight_ft1_ft5_cross, eventWeight_ft1_ft6_cross;
double eventWeight_ft2_ft3_cross, eventWeight_ft2_ft4_cross, eventWeight_ft2_ft5_cross, eventWeight_ft2_ft6_cross;
double eventWeight_ft3_ft4_cross, eventWeight_ft3_ft5_cross, eventWeight_ft3_ft6_cross;
double eventWeight_ft4_ft5_cross, eventWeight_ft4_ft6_cross;
double eventWeight_ft5_ft6_cross;

double eventWeight_fm0_fm1_cross, eventWeight_fm0_fm2_cross, eventWeight_fm0_fm3_cross, eventWeight_fm0_fm4_cross, eventWeight_fm0_fm5_cross, eventWeight_fm0_fm7_cross, eventWeight_fm0_fm8_cross, eventWeight_fm0_fm9_cross;
double eventWeight_fm1_fm2_cross, eventWeight_fm1_fm3_cross, eventWeight_fm1_fm4_cross, eventWeight_fm1_fm5_cross, eventWeight_fm1_fm7_cross, eventWeight_fm1_fm8_cross, eventWeight_fm1_fm9_cross;
double eventWeight_fm2_fm3_cross, eventWeight_fm2_fm4_cross, eventWeight_fm2_fm5_cross, eventWeight_fm2_fm7_cross, eventWeight_fm2_fm8_cross, eventWeight_fm2_fm9_cross;
double eventWeight_fm3_fm4_cross, eventWeight_fm3_fm5_cross, eventWeight_fm3_fm7_cross, eventWeight_fm3_fm8_cross, eventWeight_fm3_fm9_cross;
double eventWeight_fm4_fm5_cross, eventWeight_fm4_fm7_cross, eventWeight_fm4_fm8_cross, eventWeight_fm4_fm9_cross;
double eventWeight_fm5_fm7_cross, eventWeight_fm5_fm8_cross, eventWeight_fm5_fm9_cross;
double eventWeight_fm7_fm8_cross, eventWeight_fm7_fm9_cross;
double eventWeight_fm8_fm9_cross;

double eventWeight_fm0_quad_ll, eventWeight_fm1_quad_ll, eventWeight_fm2_quad_ll, eventWeight_fm3_quad_ll, eventWeight_fm4_quad_ll;
double eventWeight_fm5_quad_ll, eventWeight_fm7_quad_ll, eventWeight_fm8_quad_ll, eventWeight_fm9_quad_ll;
double eventWeight_fs0_quad_ll, eventWeight_fs1_quad_ll, eventWeight_fs2_quad_ll;
double eventWeight_ft0_quad_ll, eventWeight_ft1_quad_ll, eventWeight_ft2_quad_ll, eventWeight_ft3_quad_ll, eventWeight_ft4_quad_ll, eventWeight_ft5_quad_ll, eventWeight_ft6_quad_ll;

double eventWeight_fm0_quad_lt, eventWeight_fm1_quad_lt, eventWeight_fm2_quad_lt, eventWeight_fm3_quad_lt, eventWeight_fm4_quad_lt;
double eventWeight_fm5_quad_lt, eventWeight_fm7_quad_lt, eventWeight_fm8_quad_lt, eventWeight_fm9_quad_lt;
double eventWeight_fs0_quad_lt, eventWeight_fs1_quad_lt, eventWeight_fs2_quad_lt;
double eventWeight_ft0_quad_lt, eventWeight_ft1_quad_lt, eventWeight_ft2_quad_lt, eventWeight_ft3_quad_lt, eventWeight_ft4_quad_lt, eventWeight_ft5_quad_lt, eventWeight_ft6_quad_lt;

double eventWeight_fm0_quad_tl, eventWeight_fm1_quad_tl, eventWeight_fm2_quad_tl, eventWeight_fm3_quad_tl, eventWeight_fm4_quad_tl;
double eventWeight_fm5_quad_tl, eventWeight_fm7_quad_tl, eventWeight_fm8_quad_tl, eventWeight_fm9_quad_tl;
double eventWeight_fs0_quad_tl, eventWeight_fs1_quad_tl, eventWeight_fs2_quad_tl;
double eventWeight_ft0_quad_tl, eventWeight_ft1_quad_tl, eventWeight_ft2_quad_tl, eventWeight_ft3_quad_tl, eventWeight_ft4_quad_tl, eventWeight_ft5_quad_tl, eventWeight_ft6_quad_tl;

double eventWeight_fm0_quad_tt, eventWeight_fm1_quad_tt, eventWeight_fm2_quad_tt, eventWeight_fm3_quad_tt, eventWeight_fm4_quad_tt;
double eventWeight_fm5_quad_tt, eventWeight_fm7_quad_tt, eventWeight_fm8_quad_tt, eventWeight_fm9_quad_tt;
double eventWeight_fs0_quad_tt, eventWeight_fs1_quad_tt, eventWeight_fs2_quad_tt;
double eventWeight_ft0_quad_tt, eventWeight_ft1_quad_tt, eventWeight_ft2_quad_tt, eventWeight_ft3_quad_tt, eventWeight_ft4_quad_tt, eventWeight_ft5_quad_tt, eventWeight_ft6_quad_tt;

std::map<std::string, double> weightMap = {
    {"EventWeight_0", eventWeight},
    {"EventWeight_fm0_quad", eventWeight_fm0_quad},
    {"EventWeight_fm1_quad", eventWeight_fm1_quad},
    {"EventWeight_fm2_quad", eventWeight_fm2_quad},
    {"EventWeight_fm3_quad", eventWeight_fm3_quad},
    {"EventWeight_fm4_quad", eventWeight_fm4_quad},
    {"EventWeight_fm5_quad", eventWeight_fm5_quad},
    {"EventWeight_fm7_quad", eventWeight_fm7_quad},
    {"EventWeight_fm8_quad", eventWeight_fm8_quad},
    {"EventWeight_fm9_quad", eventWeight_fm9_quad},
    {"EventWeight_fs0_quad", eventWeight_fs0_quad},
    {"EventWeight_fs1_quad", eventWeight_fs1_quad},
    {"EventWeight_fs2_quad", eventWeight_fs2_quad},
    {"EventWeight_ft0_quad", eventWeight_ft0_quad},
    {"EventWeight_ft1_quad", eventWeight_ft1_quad},
    {"EventWeight_ft2_quad", eventWeight_ft2_quad},
    {"EventWeight_ft3_quad", eventWeight_ft3_quad},
    {"EventWeight_ft4_quad", eventWeight_ft4_quad},
    {"EventWeight_ft5_quad", eventWeight_ft5_quad},
    {"EventWeight_ft6_quad", eventWeight_ft6_quad}
};

std::map<std::string, double> weightMap_cross = {
    {"EventWeight_fm0_fm1_cross", eventWeight_fm0_fm1_cross},
    {"EventWeight_fm0_fm2_cross", eventWeight_fm0_fm2_cross},
    {"EventWeight_fm0_fm3_cross", eventWeight_fm0_fm3_cross},
    {"EventWeight_fm0_fm4_cross", eventWeight_fm0_fm4_cross},
    {"EventWeight_fm0_fm5_cross", eventWeight_fm0_fm5_cross},
    {"EventWeight_fm0_fm7_cross", eventWeight_fm0_fm7_cross},
    {"EventWeight_fm0_fm8_cross", eventWeight_fm0_fm8_cross},
    {"EventWeight_fm0_fm9_cross", eventWeight_fm0_fm9_cross},
    {"EventWeight_fm1_fm2_cross", eventWeight_fm1_fm2_cross},
    {"EventWeight_fm1_fm3_cross", eventWeight_fm1_fm3_cross},
    {"EventWeight_fm1_fm4_cross", eventWeight_fm1_fm4_cross},
    {"EventWeight_fm1_fm5_cross", eventWeight_fm1_fm5_cross},
    {"EventWeight_fm1_fm7_cross", eventWeight_fm1_fm7_cross},
    {"EventWeight_fm1_fm8_cross", eventWeight_fm1_fm8_cross},
    {"EventWeight_fm1_fm9_cross", eventWeight_fm1_fm9_cross},
    {"EventWeight_fm2_fm3_cross", eventWeight_fm2_fm3_cross},
    {"EventWeight_fm2_fm4_cross", eventWeight_fm2_fm4_cross},
    {"EventWeight_fm2_fm5_cross", eventWeight_fm2_fm5_cross},
    {"EventWeight_fm2_fm7_cross", eventWeight_fm2_fm7_cross},
    {"EventWeight_fm2_fm8_cross", eventWeight_fm2_fm8_cross},
    {"EventWeight_fm2_fm9_cross", eventWeight_fm2_fm9_cross},
    {"EventWeight_fm3_fm4_cross", eventWeight_fm3_fm4_cross},
    {"EventWeight_fm3_fm5_cross", eventWeight_fm3_fm5_cross},
    {"EventWeight_fm3_fm7_cross", eventWeight_fm3_fm7_cross},
    {"EventWeight_fm3_fm8_cross", eventWeight_fm3_fm8_cross},
    {"EventWeight_fm3_fm9_cross", eventWeight_fm3_fm9_cross},
    {"EventWeight_fm4_fm5_cross", eventWeight_fm4_fm5_cross},
    {"EventWeight_fm4_fm7_cross", eventWeight_fm4_fm7_cross},
    {"EventWeight_fm4_fm8_cross", eventWeight_fm4_fm8_cross},
    {"EventWeight_fm4_fm9_cross", eventWeight_fm4_fm9_cross},
    {"EventWeight_fm5_fm7_cross", eventWeight_fm5_fm7_cross},
    {"EventWeight_fm5_fm8_cross", eventWeight_fm5_fm8_cross},
    {"EventWeight_fm5_fm9_cross", eventWeight_fm5_fm9_cross},
    {"EventWeight_fm7_fm8_cross", eventWeight_fm7_fm8_cross},
    {"EventWeight_fm7_fm9_cross", eventWeight_fm7_fm9_cross},
    {"EventWeight_fm8_fm9_cross", eventWeight_fm8_fm9_cross},
    {"EventWeight_fs0_fs1_cross", eventWeight_fs0_fs1_cross},
    {"EventWeight_fs0_fs2_cross", eventWeight_fs0_fs2_cross},
    {"EventWeight_fs1_fs2_cross", eventWeight_fs1_fs2_cross},
    {"EventWeight_ft0_ft1_cross", eventWeight_ft0_ft1_cross},
    {"EventWeight_ft0_ft2_cross", eventWeight_ft0_ft2_cross},
    {"EventWeight_ft0_ft3_cross", eventWeight_ft0_ft3_cross},
    {"EventWeight_ft0_ft4_cross", eventWeight_ft0_ft4_cross},
    {"EventWeight_ft0_ft5_cross", eventWeight_ft0_ft5_cross},
    {"EventWeight_ft0_ft6_cross", eventWeight_ft0_ft6_cross},
    {"EventWeight_ft1_ft2_cross", eventWeight_ft1_ft2_cross},
    {"EventWeight_ft1_ft3_cross", eventWeight_ft1_ft3_cross},
    {"EventWeight_ft1_ft4_cross", eventWeight_ft1_ft4_cross},
    {"EventWeight_ft1_ft5_cross", eventWeight_ft1_ft5_cross},
    {"EventWeight_ft1_ft6_cross", eventWeight_ft1_ft6_cross},
    {"EventWeight_ft2_ft3_cross", eventWeight_ft2_ft3_cross},
    {"EventWeight_ft2_ft4_cross", eventWeight_ft2_ft4_cross},
    {"EventWeight_ft2_ft5_cross", eventWeight_ft2_ft5_cross},
    {"EventWeight_ft2_ft6_cross", eventWeight_ft2_ft6_cross},
    {"EventWeight_ft3_ft4_cross", eventWeight_ft3_ft4_cross},
    {"EventWeight_ft3_ft5_cross", eventWeight_ft3_ft5_cross},
    {"EventWeight_ft3_ft6_cross", eventWeight_ft3_ft6_cross},
    {"EventWeight_ft4_ft5_cross", eventWeight_ft4_ft5_cross},
    {"EventWeight_ft4_ft6_cross", eventWeight_ft4_ft6_cross},
    {"EventWeight_ft5_ft6_cross", eventWeight_ft5_ft6_cross}
};

std::map<std::string, double> weightMap_Polarisation = {
    // ll
    {"EventWeight_fm0_quad_ll", eventWeight_fm0_quad_ll},
    {"EventWeight_fm1_quad_ll", eventWeight_fm1_quad_ll},
    {"EventWeight_fm2_quad_ll", eventWeight_fm2_quad_ll},
    {"EventWeight_fm3_quad_ll", eventWeight_fm3_quad_ll},
    {"EventWeight_fm4_quad_ll", eventWeight_fm4_quad_ll},
    {"EventWeight_fm5_quad_ll", eventWeight_fm5_quad_ll},
    {"EventWeight_fm7_quad_ll", eventWeight_fm7_quad_ll},
    {"EventWeight_fm8_quad_ll", eventWeight_fm8_quad_ll},
    {"EventWeight_fm9_quad_ll", eventWeight_fm9_quad_ll},
    {"EventWeight_fs0_quad_ll", eventWeight_fs0_quad_ll},
    {"EventWeight_fs1_quad_ll", eventWeight_fs1_quad_ll},
    {"EventWeight_fs2_quad_ll", eventWeight_fs2_quad_ll},
    {"EventWeight_ft0_quad_ll", eventWeight_ft0_quad_ll},
    {"EventWeight_ft1_quad_ll", eventWeight_ft1_quad_ll},
    {"EventWeight_ft2_quad_ll", eventWeight_ft2_quad_ll},
    {"EventWeight_ft3_quad_ll", eventWeight_ft3_quad_ll},
    {"EventWeight_ft4_quad_ll", eventWeight_ft4_quad_ll},
    {"EventWeight_ft5_quad_ll", eventWeight_ft5_quad_ll},
    {"EventWeight_ft6_quad_ll", eventWeight_ft6_quad_ll},
    // lt
    {"EventWeight_fm0_quad_lt", eventWeight_fm0_quad_lt},
    {"EventWeight_fm1_quad_lt", eventWeight_fm1_quad_lt},
    {"EventWeight_fm2_quad_lt", eventWeight_fm2_quad_lt},
    {"EventWeight_fm3_quad_lt", eventWeight_fm3_quad_lt},
    {"EventWeight_fm4_quad_lt", eventWeight_fm4_quad_lt},
    {"EventWeight_fm5_quad_lt", eventWeight_fm5_quad_lt},
    {"EventWeight_fm7_quad_lt", eventWeight_fm7_quad_lt},
    {"EventWeight_fm8_quad_lt", eventWeight_fm8_quad_lt},
    {"EventWeight_fm9_quad_lt", eventWeight_fm9_quad_lt},
    {"EventWeight_fs0_quad_lt", eventWeight_fs0_quad_lt},
    {"EventWeight_fs1_quad_lt", eventWeight_fs1_quad_lt},
    {"EventWeight_fs2_quad_lt", eventWeight_fs2_quad_lt},
    {"EventWeight_ft0_quad_lt", eventWeight_ft0_quad_lt},
    {"EventWeight_ft1_quad_lt", eventWeight_ft1_quad_lt},
    {"EventWeight_ft2_quad_lt", eventWeight_ft2_quad_lt},
    {"EventWeight_ft3_quad_lt", eventWeight_ft3_quad_lt},
    {"EventWeight_ft4_quad_lt", eventWeight_ft4_quad_lt},
    {"EventWeight_ft5_quad_lt", eventWeight_ft5_quad_lt},
    {"EventWeight_ft6_quad_lt", eventWeight_ft6_quad_lt},
    // tl
    {"EventWeight_fm0_quad_tl", eventWeight_fm0_quad_tl},
    {"EventWeight_fm1_quad_tl", eventWeight_fm1_quad_tl},
    {"EventWeight_fm2_quad_tl", eventWeight_fm2_quad_tl},
    {"EventWeight_fm3_quad_tl", eventWeight_fm3_quad_tl},
    {"EventWeight_fm4_quad_tl", eventWeight_fm4_quad_tl},
    {"EventWeight_fm5_quad_tl", eventWeight_fm5_quad_tl},
    {"EventWeight_fm7_quad_tl", eventWeight_fm7_quad_tl},
    {"EventWeight_fm8_quad_tl", eventWeight_fm8_quad_tl},
    {"EventWeight_fm9_quad_tl", eventWeight_fm9_quad_tl},
    {"EventWeight_fs0_quad_tl", eventWeight_fs0_quad_tl},
    {"EventWeight_fs1_quad_tl", eventWeight_fs1_quad_tl},
    {"EventWeight_fs2_quad_tl", eventWeight_fs2_quad_tl},
    {"EventWeight_ft0_quad_tl", eventWeight_ft0_quad_tl},
    {"EventWeight_ft1_quad_tl", eventWeight_ft1_quad_tl},
    {"EventWeight_ft2_quad_tl", eventWeight_ft2_quad_tl},
    {"EventWeight_ft3_quad_tl", eventWeight_ft3_quad_tl},
    {"EventWeight_ft4_quad_tl", eventWeight_ft4_quad_tl},
    {"EventWeight_ft5_quad_tl", eventWeight_ft5_quad_tl},
    {"EventWeight_ft6_quad_tl", eventWeight_ft6_quad_tl},
    // tt
    {"EventWeight_fm0_quad_tt", eventWeight_fm0_quad_tt},
    {"EventWeight_fm1_quad_tt", eventWeight_fm1_quad_tt},
    {"EventWeight_fm2_quad_tt", eventWeight_fm2_quad_tt},
    {"EventWeight_fm3_quad_tt", eventWeight_fm3_quad_tt},
    {"EventWeight_fm4_quad_tt", eventWeight_fm4_quad_tt},
    {"EventWeight_fm5_quad_tt", eventWeight_fm5_quad_tt},
    {"EventWeight_fm7_quad_tt", eventWeight_fm7_quad_tt},
    {"EventWeight_fm8_quad_tt", eventWeight_fm8_quad_tt},
    {"EventWeight_fm9_quad_tt", eventWeight_fm9_quad_tt},
    {"EventWeight_fs0_quad_tt", eventWeight_fs0_quad_tt},
    {"EventWeight_fs1_quad_tt", eventWeight_fs1_quad_tt},
    {"EventWeight_fs2_quad_tt", eventWeight_fs2_quad_tt},
    {"EventWeight_ft0_quad_tt", eventWeight_ft0_quad_tt},
    {"EventWeight_ft1_quad_tt", eventWeight_ft1_quad_tt},
    {"EventWeight_ft2_quad_tt", eventWeight_ft2_quad_tt},
    {"EventWeight_ft3_quad_tt", eventWeight_ft3_quad_tt},
    {"EventWeight_ft4_quad_tt", eventWeight_ft4_quad_tt},
    {"EventWeight_ft5_quad_tt", eventWeight_ft5_quad_tt},
    {"EventWeight_ft6_quad_tt", eventWeight_ft6_quad_tt}
};
