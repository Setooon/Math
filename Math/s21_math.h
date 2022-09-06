#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_

#define S21_LLONG_MAX __LONG_LONG_MAX__
#define S21_LLONG_MIN (-__LONG_LONG_MAX__ - 1)
#define S21_E 2.71828182845904523536028747135266250
#define S21_M_PI \
  3.1415926535897932384626433832795028841971693993751058209749445923078164
#define S21_NAN 0.0 / 0.0
#define S21_INF 1.0 / 0.0
#define EPS_6 1e-06

long double s21_pow(double base, double expi);
long double s21_sin(double x);
long double s21_sqrt(double x);
long double s21_tan(double x);
long double s21_floor(double num);
long double s21_log(double num);
long double s21_fmod(double num_1, double num_2);
long double s21_trunc(double x);
long double s21_cos(double x);
long double s21_fabs(double x);
long double s21_ceil(double x);
long double s21_exp(double x);
int s21_abs(int a);
long double s21_acos(double x);
long double s21_asin(double a);
long double s21_atan(double x);

#endif  // SRC_S21_MATH_H_
