/*
 * Copyright (c) 2007 - 2019 Joseph Gaeddert
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

//
// Numerically-controlled oscillator (nco) with internal phase-locked
// loop (pll) implementation
//

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define NCO_STATIC_LUT_WORDBITS     32
#define NCO_STATIC_LUT_NBITS        10
#define NCO_STATIC_LUT_SIZE         (1LLU << NCO_STATIC_LUT_NBITS)
#define NCO_STATIC_LUT_INDEX_SHIFTED_PI2(index) \
    (((index) + (NCO_STATIC_LUT_SIZE >> 2)) & (NCO_STATIC_LUT_SIZE - 1))
#define NCO_STATIC_LUT_THETA_SHIFTED_PI2(theta) \
    ((theta) + (1LLU << (NCO_STATIC_LUT_WORDBITS - 2)))

#define NCO_PLL_BANDWIDTH_DEFAULT   (0.1)
#define NCO_PLL_GAIN_DEFAULT        (1000)

#define LIQUID_DEBUG_NCO            (0)

typedef T (*fp_trigonometric_func)(T);

struct vco_tab_e_s {
    T m2,c;
};
typedef struct vco_tab_e_s vco_tab_e;

struct NCO(_s) {
    liquid_ncotype  type;           // NCO type (e.g. LIQUID_VCO)

    // type == LIQUID_NCO, LIQUID_VCO_INTERP
    uint32_t        theta;          // 32-bit phase     [radians]
    uint32_t        d_theta;        // 32-bit frequency [radians/sample]

    // type == LIQUID_NCO
    T*              nco_sintab;     // sine direct look-up table

    // type == LIQUID_VCO_INTERP
    vco_tab_e*      vcoi_sintab;    // sine interpolated look-up table

    // type == LIQUID_VCO_DIRECT
    int             vcod_n;         // normalized multiplier coefficient
    unsigned int    vcod_m;         // normalized divider coefficient
    T*              vcod_sintab;    // sine direct look-up table
    T*              vcod_costab;    // cosine direct look-up table
    unsigned int    vcod_index;     // direct look-up index [0, m)

    // phase-locked loop
    T               alpha;          // frequency proportion
    T               beta;           // phase proportion
};

static const T TWO_TO_THE_31         = TIL(2147483648);
static const T TWO_TO_THE_32_MINUS_1 = TIL(4294967295);

/* constrain phase (or frequency) and convert to fixed-point */
uint32_t NCO(_constrain)(T _theta);

/* utilities used with type LIQUID_VCO_DIRECT */
void NCO(_constrain_vcod)(int *_n, unsigned int *_m);
inline int NCO(_has_vcod_lookup_table)(NCO() _q);
void NCO(_calc_vcod_lookup_table)(T*                    tab,
                                  unsigned int          size,
                                  int32_t               d_theta,
                                  T*                    scale_factor,
                                  fp_trigonometric_func f);

/* utilities used with types LIQUID_VCO_* */
void NCO(_calc_vco_tab_e)(vco_tab_e*            e,
                          uint32_t              theta,
                          int32_t               d_theta,
                          fp_trigonometric_func f);
T NCO(_fp_sin)(T x);
T NCO(_fp_cos)(T x);

// compute index for sine/cosine look-up table
unsigned int NCO(_index)(NCO() _q);

// create nco/vco object
NCO() NCO(_create)(liquid_ncotype _type)
{
    NCO() q = (NCO()) malloc(sizeof(struct NCO(_s)));
    q->type = _type;

    // initialize sine/cosine tables
    unsigned int i;
    switch (q->type) {
    case LIQUID_NCO: {
        q->vcoi_sintab = NULL;
        q->vcod_sintab = NULL;
        q->vcod_costab = NULL;
        q->nco_sintab  = (T*)malloc(NCO_STATIC_LUT_SIZE*sizeof(T));
        for (i=0; i<NCO_STATIC_LUT_SIZE; i++)
            q->nco_sintab[i] = SIN(TIL(2)*TFL(M_PI)*(T)(i)/(T)(NCO_STATIC_LUT_SIZE));
        break;
    }
    case LIQUID_VCO_INTERP: {
        q->nco_sintab  = NULL;
        q->vcod_sintab = NULL;
        q->vcod_costab = NULL;
        q->vcoi_sintab = (vco_tab_e*)malloc(NCO_STATIC_LUT_SIZE*sizeof(vco_tab_e));
        uint32_t theta = 0;
        const int32_t d_theta = (uint32_t)INT32_MAX*2/NCO_STATIC_LUT_SIZE;
        unsigned int i;
        for (i=0; i<NCO_STATIC_LUT_SIZE; i++) {
            NCO(_calc_vco_tab_e)(&(q->vcoi_sintab[i]), theta, d_theta, NCO(_fp_sin));
            theta += d_theta;
        }
        /* TODO: compensate range overflow caused by small undefined error
         *       ("scale_factor" solution from LIQUID_VCO_DIRECT doesn't fit here)
         */
        break;
    }
    case LIQUID_VCO_DIRECT: {
        q->nco_sintab  = NULL;
        q->vcoi_sintab = NULL;
        q->vcod_sintab = NULL;
        q->vcod_costab = NULL;
        break;
    }
    default:
        fprintf(stderr,"error: NCO(_create)(), unknown type : %u\n", q->type);
        exit(1);
    }

    q->vcod_n = 0;
    q->vcod_m = 0;

    // set default pll bandwidth
    NCO(_pll_set_bandwidth)(q, NCO_PLL_BANDWIDTH_DEFAULT);

    // reset object and return
    NCO(_reset)(q);
    return q;
}

// destroy nco object
void NCO(_destroy)(NCO() _q)
{
    if (!_q) {
        return;
    }

    switch (_q->type) {
    case LIQUID_NCO:
        free(_q->nco_sintab);
        break;
    case LIQUID_VCO_INTERP:
        free(_q->vcoi_sintab);
        break;
    case LIQUID_VCO_DIRECT:
        if (NCO(_has_vcod_lookup_table)(_q)) {
            free(_q->vcod_sintab);
            free(_q->vcod_costab);
        }
        break;
    default:
        break;
    }

    free(_q);
}

// Print nco object internals to stdout
void NCO(_print)(NCO() _q)
{
    if ((_q->type == LIQUID_NCO) || (_q->type == LIQUID_VCO_INTERP)) {
        const char * typestr = (_q->type == LIQUID_NCO) ? "NCO" : "VCO-interp";
        printf("nco [type: %s, phase: 0x%.8x rad, freq: 0x%.8x rad/sample]\n",
                typestr, _q->theta, _q->d_theta);
    } else if (_q->type == LIQUID_VCO_DIRECT) {
        printf("nco [type: VCO-direct, n: %11i, m: %10u, phase index: %10u]\n",
                _q->vcod_n, _q->vcod_m, _q->vcod_index);
    } else {
        printf("nco [type: <INVALID>]\n");
    }
#if LIQUID_DEBUG_NCO
    // print entire table(s)
    unsigned int i;
    switch (_q->type) {
    case LIQUID_NCO:
        for (i=0; i<NCO_STATIC_LUT_SIZE; i++)
            printf("  sintab[%4u] = %16.12f\n", i, _q->nco_sintab[i]);
        break;
    case LIQUID_VCO_INTERP:
        for (i=0; i<NCO_STATIC_LUT_SIZE; i++)
            printf("  sintab[%4u]  .m2 = %19.15e  .c = %19.15e\n",
                   i, _q->vcoi_sintab[i].m2, _q->vcoi_sintab[i].c);
        break;
    case LIQUID_VCO_DIRECT:
        if (!NCO(_has_vcod_lookup_table)(_q))
            break;
        for (i=0; i<_q->vcod_m; i++)
            printf("  [%10u]  sintab[] = %16.12f  costab[] = %16.12f\n",
                   i, _q->vcod_sintab[i], _q->vcod_costab[i]);
        break;
    default:
        break;
    }
#endif
}

// reset internal state of nco object
void NCO(_reset)(NCO() _q)
{
    switch (_q->type) {
    case LIQUID_NCO:
    case LIQUID_VCO_INTERP:
        // reset phase and frequency states
        _q->theta   = 0;
        _q->d_theta = 0;
        break;
    case LIQUID_VCO_DIRECT:
        // destroy table (if any) and reset all parameters
        if (NCO(_has_vcod_lookup_table)(_q)) {
            free(_q->vcod_sintab);
            _q->vcod_sintab = NULL;
            free(_q->vcod_costab);
            _q->vcod_costab = NULL;
        }
        _q->vcod_n     = 0;
        _q->vcod_m     = 0;
        _q->vcod_index = 0;
        break;
    default:
        break;
    }

    // reset pll filter state
    NCO(_pll_reset)(_q);
}

// set frequency of nco object
void NCO(_set_frequency)(NCO() _q,
                         T     _dtheta)
{
    if (_q->type == LIQUID_VCO_DIRECT) {
        fprintf(stderr,"error: nco_set_frequency(), "
                       "cannot be used with object type == LIQUID_VCO_DIRECT\n");
        exit(1);
    }
    _q->d_theta = NCO(_constrain)(_dtheta);
}

// adjust frequency of nco object
void NCO(_adjust_frequency)(NCO() _q,
                            T     _df)
{
    if (_q->type == LIQUID_VCO_DIRECT) {
        fprintf(stderr,"error: nco_adjust_frequency(), "
                       "cannot be used with object type == LIQUID_VCO_DIRECT\n");
        exit(1);
    }
    _q->d_theta += NCO(_constrain)(_df);
}

// set phase of nco object, constraining phase
void NCO(_set_phase)(NCO() _q,
                     T     _phi)
{
    if (_q->type == LIQUID_VCO_DIRECT) {
        fprintf(stderr,"error: nco_set_phase(), "
                       "cannot be used with object type == LIQUID_VCO_DIRECT\n");
        exit(1);
    }
    _q->theta = NCO(_constrain)(_phi);
}

// adjust phase of nco object, constraining phase
void NCO(_adjust_phase)(NCO() _q,
                        T     _dphi)
{
    if (_q->type == LIQUID_VCO_DIRECT) {
        fprintf(stderr,"error: nco_adjust_phase(), "
                       "cannot be used with object type == LIQUID_VCO_DIRECT\n");
        exit(1);
    }
    _q->theta += NCO(_constrain)(_dphi);
}

// increment internal phase of nco object
void NCO(_step)(NCO() _q)
{
    if ((_q->type == LIQUID_NCO) || (_q->type == LIQUID_VCO_INTERP)) {
        _q->theta += _q->d_theta;
    } else if ((_q->type == LIQUID_VCO_DIRECT) && NCO(_has_vcod_lookup_table)(_q)) {
        (_q->vcod_index)++;
        if (_q->vcod_index == _q->vcod_m)
            _q->vcod_index = 0;
    }
}

// get phase [radians]
T NCO(_get_phase)(NCO() _q)
{
    if (_q->type == LIQUID_VCO_DIRECT) {
        fprintf(stderr,"error: nco_get_phase(), "
                       "cannot be used with object type == LIQUID_VCO_DIRECT\n");
        exit(1);
    }
    return TIL(2)*TFL(M_PI)*(T)_q->theta / (T)(1LLU<<32);
}

// get frequency [radians/sample]
T NCO(_get_frequency)(NCO() _q)
{
    if (_q->type == LIQUID_VCO_DIRECT) {
        fprintf(stderr,"error: nco_get_frequency(), "
                       "cannot be used with object type == LIQUID_VCO_DIRECT\n");
        exit(1);
    }
    T d_theta = TIL(2)*TFL(M_PI)*(T)_q->d_theta / (T)(1LLU<<32);
    return d_theta > TFL(M_PI) ? d_theta - TIL(2)*TFL(M_PI) : d_theta;
}

// get frequency of LIQUID_VCO_DIRECT object
// [fraction defined by normalized multiplier and divider coefficients]
void NCO(_get_vcodirect_frequency)(NCO()         _q,
                                   int*          _n,
                                   unsigned int* _m)
{
    if (_q->type != LIQUID_VCO_DIRECT) {
        fprintf(stderr,"error: nco_get_vcodirect_frequency(), "
                       "cannot be used with object type != LIQUID_VCO_DIRECT\n");
        exit(1);
    }
    *_n = _q->vcod_n;
    *_m = _q->vcod_m;
}

// set frequency of LIQUID_VCO_DIRECT object
// [fraction defined by multiplier and divider coefficients]
void NCO(_set_vcodirect_frequency)(NCO()        _q,
                                   int          _n,
                                   unsigned int _m)
{
    if (_q->type != LIQUID_VCO_DIRECT) {
        fprintf(stderr,"error: nco_set_vcodirect_frequency(), "
                       "cannot be used with object type != LIQUID_VCO_DIRECT\n");
        exit(1);
    }

    if (NCO(_has_vcod_lookup_table)(_q)) {
        free(_q->vcod_sintab);
        free(_q->vcod_costab);
    }

    _q->vcod_n     = _n;
    _q->vcod_m     = _m;
    _q->vcod_index = 0;

    NCO(_constrain_vcod)(&(_q->vcod_n), &(_q->vcod_m));

    if (NCO(_has_vcod_lookup_table)(_q)) {
        _q->vcod_sintab = (T*)malloc(_q->vcod_m*sizeof(T));
        _q->vcod_costab = (T*)malloc(_q->vcod_m*sizeof(T));
        const int32_t d_theta = (int32_t)((double)TWO_TO_THE_31 * 2.0
                                          * (double)(_q->vcod_n)/(double)(_q->vcod_m));
        T scale_factor = TIL(1);
        NCO(_calc_vcod_lookup_table)(_q->vcod_sintab, _q->vcod_m,
                                     d_theta, &scale_factor, NCO(_fp_sin));
        NCO(_calc_vcod_lookup_table)(_q->vcod_costab, _q->vcod_m,
                                     d_theta, &scale_factor, NCO(_fp_cos));
        unsigned int i;
        for (i=0; i<_q->vcod_m; i++) {
            _q->vcod_sintab[i] /= scale_factor;
            _q->vcod_costab[i] /= scale_factor;
        }
    } else {
        _q->vcod_sintab = NULL;
        _q->vcod_costab = NULL;
    }
}

// compute sine, cosine internally
T NCO(_sin)(NCO() _q)
{
    unsigned int index = NCO(_index)(_q);
    T value = TIL(0);

    switch (_q->type) {
    case LIQUID_NCO: {
        value = _q->nco_sintab[index];
        break;
    }
    case LIQUID_VCO_INTERP: {
        value = _q->vcoi_sintab[index].m2 * (_q->theta >> 1)
                + _q->vcoi_sintab[index].c;
        break;
    }
    case LIQUID_VCO_DIRECT: {
        if (!NCO(_has_vcod_lookup_table)(_q))
            break;
        value = _q->vcod_sintab[index];
        break;
    }
    default:
        break;
    }

    return value;
}

T NCO(_cos)(NCO() _q)
{
    unsigned int index = NCO(_index)(_q);
    T value = TIL(1);

    switch (_q->type) {
    case LIQUID_NCO:
    case LIQUID_VCO_INTERP: {
        /* add pi/2 phase shift */
        index = NCO_STATIC_LUT_INDEX_SHIFTED_PI2(index);
        if (_q->type == LIQUID_NCO) {
            value = _q->nco_sintab[index];
        } else {
            uint32_t theta = NCO_STATIC_LUT_THETA_SHIFTED_PI2(_q->theta);
            value = _q->vcoi_sintab[index].m2 * (theta >> 1)
                    + _q->vcoi_sintab[index].c;
        }
        break;
    }
    case LIQUID_VCO_DIRECT: {
        if (!NCO(_has_vcod_lookup_table)(_q))
            break;
        value = _q->vcod_costab[index];
        break;
    }
    default:
        break;
    }

    return value;
}

// compute sin, cos of internal phase
void NCO(_sincos)(NCO() _q,
                  T *   _s,
                  T *   _c)
{
    unsigned int index = NCO(_index)(_q);

    // return result
    switch (_q->type) {
    case LIQUID_NCO:
    case LIQUID_VCO_INTERP: {
        unsigned int shifted_index = NCO_STATIC_LUT_INDEX_SHIFTED_PI2(index);
        if (_q->type == LIQUID_NCO) {
            *_s = _q->nco_sintab[index];
            *_c = _q->nco_sintab[shifted_index];
        } else {
            uint32_t shifted_theta = NCO_STATIC_LUT_THETA_SHIFTED_PI2(_q->theta);
            *_s = _q->vcoi_sintab[index].m2 * (_q->theta >> 1)
                  + _q->vcoi_sintab[index].c;
            *_c = _q->vcoi_sintab[shifted_index].m2 * (shifted_theta >> 1)
                  + _q->vcoi_sintab[shifted_index].c;
        }
        break;
    }
    case LIQUID_VCO_DIRECT: {
        if (NCO(_has_vcod_lookup_table)(_q)) {
            *_s = _q->vcod_sintab[index];
            *_c = _q->vcod_costab[index];
            break;
        }
        /* no break */
    }
    default:
        *_s = TIL(0);
        *_c = TIL(1);
        break;
    }
}

// compute complex exponential of internal phase
void NCO(_cexpf)(NCO() _q,
                 TC *  _y)
{
    T vsin;
    T vcos;
    NCO(_sincos)(_q, &vsin, &vcos);
    *_y = vcos + _Complex_I*vsin;
}

// pll methods

// reset pll state, retaining base frequency
void NCO(_pll_reset)(NCO() _q)
{
}

// set pll bandwidth
void NCO(_pll_set_bandwidth)(NCO() _q,
                             T     _bw)
{
    // validate input
    if (_bw < TIL(0)) {
        fprintf(stderr,"error: nco_pll_set_bandwidth(), bandwidth must be positive\n");
        exit(1);
    }

    _q->alpha = _bw;                // frequency proportion
    _q->beta  = sqrtf(_q->alpha);   // phase proportion
}

// advance pll phase
//  _q      :   nco object
//  _dphi   :   phase error
void NCO(_pll_step)(NCO() _q,
                    T     _dphi)
{
    if (_q->type == LIQUID_VCO_DIRECT) {
        fprintf(stderr,"error: nco_pll_step(), "
                       "cannot be used with object type == LIQUID_VCO_DIRECT\n");
        exit(1);
    }

    // increase frequency proportional to error
    NCO(_adjust_frequency)(_q, _dphi*_q->alpha);

    // increase phase proportional to error
    NCO(_adjust_phase)(_q, _dphi*_q->beta);

    // constrain frequency
    //NCO(_constrain_frequency)(_q);
}

// mixing functions

// Rotate input vector up by NCO angle, y = x exp{+j theta}
//  _q      :   nco object
//  _x      :   input sample
//  _y      :   output sample
void NCO(_mix_up)(NCO() _q,
                  TC    _x,
                  TC *  _y)
{
    // compute complex phasor
    TC v;
    NCO(_cexpf)(_q, &v);

    // rotate input
    *_y = _x * v;
}

// Rotate input vector down by NCO angle, y = x exp{-j theta}
//  _q      :   nco object
//  _x      :   input sample
//  _y      :   output sample
void NCO(_mix_down)(NCO() _q,
                    TC _x,
                    TC *_y)
{
    // compute complex phasor
    TC v;
    NCO(_cexpf)(_q, &v);

    // rotate input (negative direction)
    *_y = _x * conj(v);
}


// Rotate input vector array up by NCO angle:
//      y(t) = x(t) exp{+j (f*t + theta)}
// TODO : implement NCO/VCO-specific versions
//  _q      :   nco object
//  _x      :   input array [size: _n x 1]
//  _y      :   output sample [size: _n x 1]
//  _n      :   number of input, output samples
void NCO(_mix_block_up)(NCO() _q,
                        TC *_x,
                        TC *_y,
                        unsigned int _n)
{
    unsigned int i;
    // FIXME: this method should be more efficient but is causing occasional
    //        errors so instead favor slower but more reliable algorithm
    //        (anyway it must be rewritten to work with LIQUID_VCO_DIRECT type)
#if 0
    T theta =   _q->theta;
    T d_theta = _q->d_theta;
    for (i=0; i<_n; i++) {
        // multiply _x[i] by [cos(theta) + _Complex_I*sin(theta)]
        _y[i] = _x[i] * liquid_cexpjf(theta);
        
        theta += d_theta;
    }

    NCO(_set_phase)(_q, theta);
#else
    for (i=0; i<_n; i++) {
        // mix single sample up
        NCO(_mix_up)(_q, _x[i], &_y[i]);

        // step NCO phase
        NCO(_step)(_q);
    }
#endif
}

// Rotate input vector array down by NCO angle:
//      y(t) = x(t) exp{-j (f*t + theta)}
// TODO : implement NCO/VCO-specific versions
//  _q      :   nco object
//  _x      :   input array [size: _n x 1]
//  _y      :   output sample [size: _n x 1]
//  _n      :   number of input, output samples
void NCO(_mix_block_down)(NCO() _q,
                          TC *_x,
                          TC *_y,
                          unsigned int _n)
{
    unsigned int i;
    // FIXME: this method should be more efficient but is causing occasional
    //        errors so instead favor slower but more reliable algorithm
    //        (anyway it must be rewritten to work with LIQUID_VCO_DIRECT type)
#if 0
    T theta =   _q->theta;
    T d_theta = _q->d_theta;
    for (i=0; i<_n; i++) {
        // multiply _x[i] by [cos(-theta) + _Complex_I*sin(-theta)]
        _y[i] = _x[i] * liquid_cexpjf(-theta);
        
        theta += d_theta;
    }

    NCO(_set_phase)(_q, theta);
#else
    for (i=0; i<_n; i++) {
        // mix single sample down
        NCO(_mix_down)(_q, _x[i], &_y[i]);

        // step NCO phase
        NCO(_step)(_q);
    }
#endif
}

//
// internal methods
//

uint32_t NCO(_constrain)(T _theta)
{
    // divide value by 2*pi and compute modulo
    T p = _theta * 0.159154943091895;   // 1/(2 pi) ~ 0.159154943091895

    // extract fractional part of p
    T fpart = p - ((long)p);    // fpart is in (-1,1)

    // ensure fpart is in [0,1)
    if (fpart < 0.) fpart += 1.;

    // map to range of precision needed
    return (uint32_t)(fpart * 0xffffffff);
}

void NCO(_constrain_vcod)(int *_n, unsigned int *_m)
{
    if ((*_m) == 0)
        return;

    /* fold 'n' into [-'m'/2, 'm'/2) range */
    *_n %= (int)(*_m);
    if ((T)(abs(*_n)) >= (T)(*_m)/TIL(2)) {
        int sign = ((*_n) > 0) ? -1 : 1;
        *_n = sign*((*_m) - (unsigned int)abs(*_n));
    }

    /* try optimize values via reducing them by common denominators */
    // with base 10
    while ((((*_n) % 10) == 0) && (((*_m) % 10) == 0)) {
        *_n /= 10;
        *_m /= 10;
    }
    // with base 2
    while ((((*_n) & 1) == 0) && (((*_m) & 1) == 0)) {
        *_n >>= 1;
        *_m >>= 1;
    }
}

inline int NCO(_has_vcod_lookup_table)(NCO() _q)
{
    return ((_q->vcod_n != 0) && (_q->vcod_m > 0));
}

void NCO(_calc_vcod_lookup_table)(T*                    tab,
                                  unsigned int          size,
                                  int32_t               d_theta,
                                  T*                    scale_factor,
                                  fp_trigonometric_func f)
{
    uint32_t theta = 0;
    unsigned int i;
    vco_tab_e e;
    for (i=0; i<size; i++) {
        NCO(_calc_vco_tab_e)(&e, theta, d_theta, f);

        tab[i] = e.m2 * (theta >> 1) + e.c;
        if (FABS(tab[i]) > (*scale_factor))
            *scale_factor = FABS(tab[i]);

        theta += d_theta;
    }
}

void NCO(_calc_vco_tab_e)(vco_tab_e*            e,
                          uint32_t              theta,
                          int32_t               d_theta,
                          fp_trigonometric_func f)
{
    const T a = (T)(theta);
    const T b = (T)((int64_t)(theta) + d_theta);
    e->m2 = TIL(2) * ((f(b)-f(a)) / (b-a));
    e->c  = (TIL(3)*a+b)*(f(a)-f(b))/(TIL(4)*(b-a)) + (f((a+b)/TIL(2)) + f(a))/TIL(2);
}

T NCO(_fp_sin)(T x)
{
    return SIN(x * TFL(M_PI) / TWO_TO_THE_31);
}

T NCO(_fp_cos)(T x)
{
    return COS(x * TFL(M_PI) / TWO_TO_THE_31);
}

unsigned int NCO(_index)(NCO() _q)
{
    /* TODO: LIQUID_NCO and LIQUID_VCO_INTERP are expected to share same code.
     *       But "appropriate" rounding for LIQUID_VCO_INTERP type causes
     *       phase breaks at some wrap points.
     *       Not sure, so just keep it for type it originated from...
     */
    switch (_q->type) {
    case LIQUID_NCO:
        //round down
        //return (_q->theta >> (NCO_STATIC_LUT_WORDBITS-NCO_STATIC_LUT_NBITS))
        //        & (NCO_STATIC_LUT_SIZE - 1);
        // round appropriately
        return ((_q->theta + (1<<(NCO_STATIC_LUT_WORDBITS-NCO_STATIC_LUT_NBITS-1)))
                 >> (NCO_STATIC_LUT_WORDBITS-NCO_STATIC_LUT_NBITS))
                & (NCO_STATIC_LUT_SIZE-1);
    case LIQUID_VCO_INTERP:
        return (_q->theta >> (NCO_STATIC_LUT_WORDBITS-NCO_STATIC_LUT_NBITS))
                & (NCO_STATIC_LUT_SIZE-1);
    case LIQUID_VCO_DIRECT:
        return _q->vcod_index;
    default:
        return 0;
    }
}

