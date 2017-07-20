/*
 * Copyright (c) 1999
 * Silicon Graphics Computer Systems, Inc.
 *
 * Copyright (c) 1999
 * Boris Fomitchev
 *
 * Copyright (c) 2012
 * Andreas Kloeckner
 *
 * This material is provided "as is", with absolutely no warranty expressed
 * or implied. Any use is at your own risk.
 *
 * Permission to use or copy this software for any purpose is hereby granted
 * without fee, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is granted,
 * provided the above notices are retained, and a notice that the code was
 * modified is included with the above copyright notice.
 *
 */

// This file is available for inclusion in pyopencl kernels and provides
// complex types 'cfloat_t' and 'cdouble_t', along with a number of special
// functions as visible below, e.g. cdouble_log(z).
//
// Under the hood, the complex types are simply float2 and double2.
// Note that native (operator-based) addition (float + float2) and
// multiplication (float2*float1) is defined for these types,
// but do not match the rules of complex arithmetic.

#pragma once

#ifdef PYOPENCL_COMPLEX_USE_VECTOR
#define PYOPENCL_COMPLEX_REAL x
#define PYOPENCL_COMPLEX_IMAG y
#else
#define PYOPENCL_COMPLEX_REAL real
#define PYOPENCL_COMPLEX_IMAG imag
#endif

#define PYOPENCL_DECLARE_COMPLEX_TYPE_INT(REAL_TP, REAL_3LTR, TPROOT, TP) \
  \
  inline REAL_TP TPROOT##_real(TP a) { return a.PYOPENCL_COMPLEX_REAL; } \
  inline REAL_TP TPROOT##_imag(TP a) { return a.PYOPENCL_COMPLEX_IMAG; }        \
  inline REAL_TP TPROOT##_abs(TP a) { return hypot(TPROOT##_real(a), TPROOT##_imag(a)); }   \
  inline REAL_TP TPROOT##_abs_squared(TP a) { return TPROOT##_real(a) * TPROOT##_real(a) + TPROOT##_imag(a) * TPROOT##_imag(a); } \
  \
  inline TP TPROOT##_new(REAL_TP real, REAL_TP imag)  \
  { \
    TP result; \
    result.PYOPENCL_COMPLEX_REAL = real;		\
    result.PYOPENCL_COMPLEX_IMAG = imag;		\
    return result; \
  } \
  \
  inline TP TPROOT##_fromreal(REAL_TP real)     \
  { \
    TP result; \
    result.PYOPENCL_COMPLEX_REAL = real;		\
    result.PYOPENCL_COMPLEX_IMAG = 0;			\
    return result; \
  } \
  \
  \
  inline TP TPROOT##_neg(TP a) { return TPROOT##_new(-TPROOT##_real(a), -TPROOT##_imag(a)); } \
  inline TP TPROOT##_conj(TP a) { return TPROOT##_new(TPROOT##_real(a), -TPROOT##_imag(a)); } \
  \
  inline TP TPROOT##_add(TP a, TP b)            \
  { \
      return TPROOT##_new(TPROOT##_real(a) + TPROOT##_real(b), TPROOT##_imag(a) + TPROOT##_imag(b)); \
    ; \
  } \
  inline TP TPROOT##_addr(TP a, REAL_TP b)      \
  { \
    return TPROOT##_new(b+TPROOT##_real(a), TPROOT##_imag(a)); \
  } \
  inline TP TPROOT##_radd(REAL_TP a, TP b)      \
  { \
      return TPROOT##_new(a+TPROOT##_real(b), TPROOT##_imag(b));	\
  } \
  \
  inline TP TPROOT##_sub(TP a, TP b)            \
  { \
    return TPROOT##_new(TPROOT##_real(a) - TPROOT##_real(b), TPROOT##_imag(a) - TPROOT##_imag(b)); \
    ; \
  } \
  \
  inline TP TPROOT##_mul(TP a, TP b)            \
  { \
    return TPROOT##_new( \
        TPROOT##_real(a)*TPROOT##_real(b) - TPROOT##_imag(a)*TPROOT##_imag(b), \
        TPROOT##_real(a)*TPROOT##_imag(b) + TPROOT##_imag(a)*TPROOT##_real(b)); \
  } \
  \
  inline TP TPROOT##_mulr(TP a, REAL_TP b)      \
  { \
    return TPROOT##_new(TPROOT##_real(a)*b, TPROOT##_imag(a)*b); \
  } \
  \
  inline TP TPROOT##_rmul(REAL_TP a, TP b)      \
  { \
    return TPROOT##_new(a*TPROOT##_real(b), a*TPROOT##_imag(b)); \
  } \
  \
  inline TP TPROOT##_rdivide(REAL_TP z1, TP z2) \
  { \
    if (fabs(TPROOT##_real(z2)) <= fabs(TPROOT##_imag(z2))) { \
      REAL_TP ratio = TPROOT##_real(z2) / TPROOT##_imag(z2); \
      REAL_TP denom = TPROOT##_imag(z2) * (1 + ratio * ratio); \
      return TPROOT##_new((z1 * ratio) / denom, - z1 / denom); \
    } \
    else { \
      REAL_TP ratio = TPROOT##_imag(z2) / TPROOT##_real(z2); \
      REAL_TP denom = TPROOT##_real(z2) * (1 + ratio * ratio); \
      return TPROOT##_new(z1 / denom, - (z1 * ratio) / denom); \
    } \
  } \
  \
  inline TP TPROOT##_divide(TP z1, TP z2)       \
  { \
    REAL_TP ratio, denom, a, b, c, d; \
    \
    if (fabs(TPROOT##_real(z2)) <= fabs(TPROOT##_imag(z2))) { \
      ratio = TPROOT##_real(z2) / TPROOT##_imag(z2); \
      denom = TPROOT##_imag(z2); \
      a = TPROOT##_imag(z1); \
      b = TPROOT##_real(z1); \
      c = -TPROOT##_real(z1); \
      d = TPROOT##_imag(z1); \
    } \
    else { \
      ratio = TPROOT##_imag(z2) / TPROOT##_real(z2); \
      denom = TPROOT##_real(z2); \
      a = TPROOT##_real(z1); \
      b = TPROOT##_imag(z1); \
      c = TPROOT##_imag(z1); \
      d = -TPROOT##_real(z1); \
    } \
    denom *= (1 + ratio * ratio); \
    return TPROOT##_new( \
       (a + b * ratio) / denom, \
       (c + d * ratio) / denom); \
  } \
  \
  inline TP TPROOT##_divider(TP a, REAL_TP b)   \
  { \
    return TPROOT##_new(TPROOT##_real(a)/b, TPROOT##_imag(a)/b); \
  } \
  \
  inline TP TPROOT##_pow(TP a, TP b)            \
  { \
    REAL_TP logr = log(hypot(TPROOT##_real(a), TPROOT##_imag(a))); \
    REAL_TP logi = atan2(TPROOT##_imag(a), TPROOT##_real(a)); \
    REAL_TP x = exp(logr * TPROOT##_real(b) - logi * TPROOT##_imag(b)); \
    REAL_TP y = logr * TPROOT##_imag(b) + logi * TPROOT##_real(b); \
    \
    REAL_TP cosy; \
    REAL_TP siny = sincos(y, &cosy); \
    return TPROOT##_new(x*cosy, x*siny); \
  } \
  \
  inline TP TPROOT##_powr(TP a, REAL_TP b)      \
  { \
    REAL_TP logr = log(hypot(TPROOT##_real(a), TPROOT##_imag(a))); \
    REAL_TP logi = atan2(TPROOT##_imag(a), TPROOT##_real(a)); \
    REAL_TP x = exp(logr * b); \
    REAL_TP y = logi * b; \
    \
    REAL_TP cosy; \
    REAL_TP siny = sincos(y, &cosy); \
    \
    return TPROOT##_new(x * cosy, x*siny); \
  } \
  \
  inline TP TPROOT##_rpow(REAL_TP a, TP b)      \
  { \
    REAL_TP logr = log(a); \
    REAL_TP x = exp(logr * TPROOT##_real(b)); \
    REAL_TP y = logr * TPROOT##_imag(b); \
    \
    REAL_TP cosy; \
    REAL_TP siny = sincos(y, &cosy); \
    return TPROOT##_new(x * cosy, x * siny); \
  } \
  \
  inline TP TPROOT##_sqrt(TP a)                 \
  { \
    REAL_TP re = TPROOT##_real(a); \
    REAL_TP im = TPROOT##_imag(a); \
    REAL_TP mag = hypot(re, im); \
    TP result; \
    \
    if (mag == 0.f) { \
      result.PYOPENCL_COMPLEX_REAL = result.PYOPENCL_COMPLEX_IMAG = 0.f; \
    } else if (re > 0.f) { \
      result.PYOPENCL_COMPLEX_REAL = sqrt(0.5f * (mag + re)); \
      result.PYOPENCL_COMPLEX_IMAG = im/TPROOT##_real(result)/2.f; \
    } else { \
      result.PYOPENCL_COMPLEX_IMAG = sqrt(0.5f * (mag - re)); \
      if (im < 0.f) \
        result.PYOPENCL_COMPLEX_IMAG = - TPROOT##_imag(result); \
      result.PYOPENCL_COMPLEX_REAL = im/TPROOT##_imag(result)/2.f; \
    } \
    return result; \
  } \
  \
  inline TP TPROOT##_exp(TP a) \
  { \
    REAL_TP expr = exp(TPROOT##_real(a)); \
    REAL_TP cosi; \
    REAL_TP sini = sincos(TPROOT##_imag(a), &cosi); \
    return TPROOT##_new(expr * cosi, expr * sini); \
  } \
  \
  inline TP TPROOT##_log(TP a)                                                 \
  { return TPROOT##_new(log(hypot(TPROOT##_real(a), TPROOT##_imag(a))), atan2(TPROOT##_imag(a), TPROOT##_real(a))); } \
  \
  inline TP TPROOT##_sin(TP a) \
  { \
    REAL_TP cosr; \
    REAL_TP sinr = sincos(TPROOT##_real(a), &cosr); \
    return TPROOT##_new(sinr*cosh(TPROOT##_imag(a)), cosr*sinh(TPROOT##_imag(a))); \
  } \
  \
  inline TP TPROOT##_cos(TP a) \
  { \
    REAL_TP cosr; \
    REAL_TP sinr = sincos(TPROOT##_real(a), &cosr); \
    return TPROOT##_new(cosr*cosh(TPROOT##_imag(a)), -sinr*sinh(TPROOT##_imag(a))); \
  } \
  \
  inline TP TPROOT##_tan(TP a) \
  { \
    REAL_TP re2 = 2.f * TPROOT##_real(a); \
    REAL_TP im2 = 2.f * TPROOT##_imag(a); \
    \
    const REAL_TP limit = log(REAL_3LTR##_MAX); \
    \
    if (fabs(im2) > limit) \
      return TPROOT##_new(0.f, (im2 > 0 ? 1.f : -1.f)); \
    else \
    { \
      REAL_TP den = cos(re2) + cosh(im2); \
      return TPROOT##_new(sin(re2) / den, sinh(im2) / den); \
    } \
  } \
  \
  inline TP TPROOT##_sinh(TP a) \
  { \
    REAL_TP cosi; \
    REAL_TP sini = sincos(TPROOT##_imag(a), &cosi); \
    return TPROOT##_new(sinh(TPROOT##_real(a))*cosi, cosh(TPROOT##_real(a))*sini); \
  } \
  \
  inline TP TPROOT##_cosh(TP a) \
  { \
    REAL_TP cosi; \
    REAL_TP sini = sincos(TPROOT##_imag(a), &cosi); \
    return TPROOT##_new(cosh(TPROOT##_real(a))*cosi, sinh(TPROOT##_real(a))*sini); \
  } \
  \
  inline TP TPROOT##_tanh(TP a) \
  { \
    REAL_TP re2 = 2.f * TPROOT##_real(a); \
    REAL_TP im2 = 2.f * TPROOT##_imag(a); \
    \
    const REAL_TP limit = log(REAL_3LTR##_MAX); \
    \
    if (fabs(re2) > limit) \
      return TPROOT##_new((re2 > 0 ? 1.f : -1.f), 0.f); \
    else \
    { \
      REAL_TP den = cosh(re2) + cos(im2); \
      return TPROOT##_new(sinh(re2) / den, sin(im2) / den); \
    } \
  } \


#ifdef PYOPENCL_COMPLEX_USE_VECTOR
#define PYOPENCL_DECLARE_COMPLEX_TYPE_INT0(BASE) \
  typedef BASE##2 c##BASE##_t;
#else
#define PYOPENCL_DECLARE_COMPLEX_TYPE_INT0(BASE) \
  typedef union \
  { \
    struct { BASE x, y; }; \
    struct { BASE real, imag; }; \
  } c##BASE##_t;
#endif

#define PYOPENCL_DECLARE_COMPLEX_TYPE(BASE, BASE_3LTR) \
  PYOPENCL_DECLARE_COMPLEX_TYPE_INT0(BASE) \
  PYOPENCL_DECLARE_COMPLEX_TYPE_INT(BASE, BASE_3LTR, c##BASE, c##BASE##_t)

PYOPENCL_DECLARE_COMPLEX_TYPE(float, FLT);
#define cfloat_cast(a) cfloat_new((a).PYOPENCL_COMPLEX_REAL, (a).PYOPENCL_COMPLEX_IMAG)

#ifdef PYOPENCL_DEFINE_CDOUBLE
PYOPENCL_DECLARE_COMPLEX_TYPE(double, DBL);
#define cdouble_cast(a) cdouble_new((a).PYOPENCL_COMPLEX_REAL, (a).PYOPENCL_COMPLEX_IMAG)
#endif
