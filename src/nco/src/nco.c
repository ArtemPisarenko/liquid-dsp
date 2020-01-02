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

#define NCO_PLL_BANDWIDTH_DEFAULT   (0.1)
#define NCO_PLL_GAIN_DEFAULT        (1000)

#define LIQUID_DEBUG_NCO            (0)

typedef T (*fp_trigonometric_func)(T);

struct NCO(_s) {
    liquid_ncotype  type;           // NCO type (e.g. LIQUID_VCO)

    // type == LIQUID_NCO
    T*              nco_sintab;     // sine look-up table
    uint32_t        theta;          // 32-bit phase     [radians]
    uint32_t        d_theta;        // 32-bit frequency [radians/sample]

    // type == LIQUID_VCO
    uint32_t        vco_size;      // size of look-up tables
    T*              vco_sintab;    // sine direct look-up table
    T*              vco_costab;    // cosine direct look-up table
    uint32_t        vco_index;     // direct look-up index [0, size]

    // phase-locked loop
    T               alpha;          // frequency proportion
    T               beta;           // phase proportion
};

static const T TWO_TO_THE_31 = 2147483648.0;
static const T TWO_TO_THE_32_MINUS_1 = 4294967295.0;

// constrain phase (or frequency) and convert to fixed-point
uint32_t NCO(_constrain)(float _theta);

void NCO(_constrain_precise)(int *_n, unsigned int *_m);

void NCO(_calc_precise_lookup_table)(T*                    tab,
                                     unsigned int          size,
                                     int32_t               d_theta,
                                     T*                    scale_factor,
                                     fp_trigonometric_func f);

T NCO(_fp_sin)(T x);

T NCO(_fp_cos)(T x);

// compute index for sine look-up table
unsigned int NCO(_index)(NCO() _q);

// create nco/vco object
NCO() NCO(_create)(liquid_ncotype _type)
{
    NCO() q = (NCO()) malloc(sizeof(struct NCO(_s)));
    q->type = _type;

    if (_type == LIQUID_NCO) {
        // initialize sine table
        q->nco_sintab = (T*)malloc(1024*sizeof(T));
        unsigned int i;
        for (i=0; i<1024; i++)
            q->nco_sintab[i] = SIN(2.0f*M_PI*(float)(i)/1024.0f);
    } else {
        q->nco_sintab = NULL;
    }

    q->vco_size   = 0;
    q->vco_sintab = NULL;
    q->vco_costab = NULL;

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

    if (_q->type == LIQUID_NCO)
        free(_q->nco_sintab);
    if (_q->vco_size != 0) {
        free(_q->vco_sintab);
        free(_q->vco_costab);
    }

    free(_q);
}

// Print nco object internals to stdout
void NCO(_print)(NCO() _q)
{
    if (_q->type == LIQUID_NCO) {
        printf("nco [type: NCO, phase: 0x%.8x rad, freq: 0x%.8x rad/sample]\n",
                _q->theta, _q->d_theta);
    } else if (_q->type == LIQUID_VCO) {
        printf("nco [type: VCO, m/size: 0x%.8x, phase index: 0x%.8x]\n",
                _q->vco_size, _q->vco_index);
    } else {
        printf("nco [type: <INVALID>]\n");
    }
#if LIQUID_DEBUG_NCO
    // print entire table
    unsigned int i;
    if (_q->type == LIQUID_NCO) {
        for (i=0; i<1024; i++)
            printf("  sintab[%4u] = %16.12f\n", i, _q->nco_sintab[i]);
    } else if (_q->type == LIQUID_VCO) {
        for (i=0; i<_q->vco_size; i++)
            printf("  %10u  sintab[i] = %16.12f  costab[i] = %16.12f\n",
                   i, _q->vco_sintab[i], _q->vco_costab[i]);
    }
#endif
}

// reset internal state of nco object
void NCO(_reset)(NCO() _q)
{
    // reset phase and frequency states
    _q->theta   = 0;
    _q->d_theta = 0;

    if (_q->vco_size != 0) {
        _q->vco_size = 0;
        free(_q->vco_sintab);
        _q->vco_sintab = NULL;
        free(_q->vco_costab);
        _q->vco_costab = NULL;
    }
    _q->vco_index = 0;

    // reset pll filter state
    NCO(_pll_reset)(_q);
}

// set frequency of nco object
void NCO(_set_frequency)(NCO() _q,
                         T     _dtheta)
{
    if (_q->type != LIQUID_NCO) {
        fprintf(stderr,"error: nco_set_frequency(), cannot be used with object type != LIQUID_NCO\n");
        exit(1);
    }
    _q->d_theta = NCO(_constrain)(_dtheta);
}

// adjust frequency of nco object
void NCO(_adjust_frequency)(NCO() _q,
                            T     _df)
{
    if (_q->type != LIQUID_NCO) {
        fprintf(stderr,"error: nco_adjust_frequency(), cannot be used with object type != LIQUID_NCO\n");
        exit(1);
    }
    _q->d_theta += NCO(_constrain)(_df);
}

// set phase of nco object, constraining phase
void NCO(_set_phase)(NCO() _q,
                     T     _phi)
{
    if (_q->type != LIQUID_NCO) {
        fprintf(stderr,"error: nco_set_phase(), cannot be used with object type != LIQUID_NCO\n");
        exit(1);
    }
    _q->theta = NCO(_constrain)(_phi);
}

// adjust phase of nco object, constraining phase
void NCO(_adjust_phase)(NCO() _q,
                        T     _dphi)
{
    if (_q->type != LIQUID_NCO) {
        fprintf(stderr,"error: nco_adjust_phase(), cannot be used with object type != LIQUID_NCO\n");
        exit(1);
    }
    _q->theta += NCO(_constrain)(_dphi);
}

// increment internal phase of nco object
void NCO(_step)(NCO() _q)
{
    if (_q->type == LIQUID_NCO) {
        _q->theta += _q->d_theta;
    } else if ((_q->type == LIQUID_VCO) && (_q->vco_size != 0)) {
        (_q->vco_index)++;
        if (_q->vco_index == _q->vco_size)
            _q->vco_index = 0;
    }
}

// get phase [radians]
T NCO(_get_phase)(NCO() _q)
{
    if (_q->type != LIQUID_NCO) {
        fprintf(stderr,"error: nco_get_phase(), cannot be used with object type != LIQUID_NCO\n");
        exit(1);
    }
    return 2.0f*M_PI*(float)_q->theta / (float)(1LLU<<32);
}

// get frequency [radians/sample]
T NCO(_get_frequency)(NCO() _q)
{
    if (_q->type != LIQUID_NCO) {
        fprintf(stderr,"error: nco_get_frequency(), cannot be used with object type != LIQUID_NCO\n");
        exit(1);
    }
    float d_theta = 2.0f*M_PI*(float)_q->d_theta / (float)(1LLU<<32);
    return d_theta > M_PI ? d_theta - 2*M_PI : d_theta;
}

void NCO(_set_precise_frequency)(NCO()        _q,
                                 int          _n,
                                 unsigned int _m)
{
    if (_q->type != LIQUID_VCO) {
        fprintf(stderr,"error: nco_set_precise_frequency(), cannot be used with object type != LIQUID_VCO\n");
        exit(1);
    }

    if (_q->vco_size != 0) {
        free(_q->vco_sintab);
        free(_q->vco_costab);
    }

    if (_m > 0)
        NCO(_constrain_precise)(&_n, &_m);

    _q->vco_index = 0;

    if ((_n != 0) && (_m > 0)) {
        _q->vco_size   = _m;
        _q->vco_sintab = (T*)malloc(_m*sizeof(T));
        _q->vco_costab = (T*)malloc(_m*sizeof(T));
        const int32_t d_theta = (int32_t)((T)TWO_TO_THE_31 * ( 2.0 * (T)(_n) / (T)(_m)));
        T scale_factor = 1.0;
        NCO(_calc_precise_lookup_table)(_q->vco_sintab, _m, d_theta, &scale_factor, NCO(_fp_sin));
        NCO(_calc_precise_lookup_table)(_q->vco_costab, _m, d_theta, &scale_factor, NCO(_fp_cos));
        unsigned int i;
        for (i=0; i<_m; i++) {
            _q->vco_sintab[i] /= scale_factor;
            _q->vco_costab[i] /= scale_factor;
        }
    } else {
        _q->vco_size  = 0;
        _q->vco_sintab = NULL;
        _q->vco_costab = NULL;
    }
}

// compute sine, cosine internally
T NCO(_sin)(NCO() _q)
{
    T value = 0.0;
    if (_q->type == LIQUID_NCO) {
        unsigned int index = NCO(_index)(_q);
        value = _q->nco_sintab[index];
    } else if ((_q->type == LIQUID_VCO) && (_q->vco_size != 0)) {
        value = _q->vco_sintab[_q->vco_index];
    }
    return value;
}

T NCO(_cos)(NCO() _q)
{
    T value = 1.0;
    if (_q->type == LIQUID_NCO) {
        // add pi/2 phase shift
        unsigned int index = (NCO(_index)(_q) + 256) & 0x3ff;
        value = _q->nco_sintab[index];
    } else if ((_q->type == LIQUID_VCO) && (_q->vco_size != 0)) {
        value = _q->vco_costab[_q->vco_index];
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
    if (_q->type == LIQUID_NCO) {
        *_s = _q->nco_sintab[(index    )        ];
        *_c = _q->nco_sintab[(index+256) & 0x3ff];
    } else if ((_q->type == LIQUID_VCO) && (_q->vco_size != 0)) {
        *_s = _q->vco_sintab[index];
        *_c = _q->vco_costab[index];
    } else {
        *_s = 0.0;
        *_c = 1.0;
    }
}

// compute complex exponential of internal phase
void NCO(_cexpf)(NCO() _q,
                 TC *  _y)
{
    float vsin;
    float vcos;
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
    if (_bw < 0.0f) {
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
    if (_q->type != LIQUID_NCO) {
        fprintf(stderr,"error: nco_pll_step(), cannot be used with object type != LIQUID_NCO\n");
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

// constrain phase (or frequency) and convert to fixed-point
uint32_t NCO(_constrain)(float _theta)
{
    // divide value by 2*pi and compute modulo
    float p = _theta * 0.159154943091895;   // 1/(2 pi) ~ 0.159154943091895

    // extract fractional part of p
    float fpart = p - ((long)p);    // fpart is in (-1,1)

    // ensure fpart is in [0,1)
    if (fpart < 0.) fpart += 1.;

    // map to range of precision needed
    return (uint32_t)(fpart * 0xffffffff);
}

void NCO(_constrain_precise)(int *_n, unsigned int *_m)
{
    /* fold 'n' into [-'m'/2, 'm'/2) range */
    *_n %= (int)(*_m);
    if ((T)(abs(*_n)) >= (T)(*_m)/2.0) {
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

void NCO(_calc_precise_lookup_table)(T*                    tab,
                                     unsigned int          size,
                                     int32_t               d_theta,
                                     T*                    scale_factor,
                                     fp_trigonometric_func f)
{
    uint32_t theta = 0;
    unsigned int i;
    for (i=0; i<size; i++) {
        const T a = (T)(theta);
        const T b = (T)((int64_t)(theta) + d_theta);

        const T m = ((f(b)-f(a)) / (b-a));
        const T c = (3.0*a+b)*(f(a)-f(b))/(4.0*(b-a)) + (f((a+b)/2.0) + f(a))/2.0;

        tab[i] = m * theta + c;

        if (FABS(tab[i]) > (*scale_factor))
            *scale_factor = FABS(tab[i]);

        theta += d_theta;
    }
}

T NCO(_fp_sin)(T x)
{
    return SIN(x * M_PI / TWO_TO_THE_31);
}

T NCO(_fp_cos)(T x)
{
    return COS(x * M_PI / TWO_TO_THE_31);
}

// compute index for sine look-up table
unsigned int NCO(_index)(NCO() _q)
{
    unsigned int value = 0;
    if (_q->type == LIQUID_NCO) {
        //return (_q->theta >> 22) & 0x3ff; // round down
        value = ((_q->theta + (1<<21)) >> 22) & 0x3ff; // round appropriately
    } else if ((_q->type == LIQUID_VCO) && (_q->vco_size != 0)) {
        value = _q->vco_index;
    }
    return value;
}

