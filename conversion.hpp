#ifndef conversion_hpp_INCLUDED
#define conversion_hpp_INCLUDED

#include "common.hpp"

#include <algorithm>
#include <iostream>

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define DBLFLT_EPS ((double)FLT_EPSILON)
#define MAX_IMBA_DEP ((uint32_t)INT32_MAX)

// http://snap.stanford.edu/class/cs224w-2015/slides/04-powerlaws.pdf
// https://cs.brynmawr.edu/Courses/cs380/spring2013/section02/slides/10_ScaleFreeNetworks.pdf
// #define PLEXP (3.0) // Power-law exponent of the degree distribution of the BA network
// #define SYS_QMAX (25.0) // Is it the PHRED-scaled probability that a germline event is an outlier in terms of allele fraction?

template <class T>
T
div_by_20(T n) {
    return n / 20;    
}

auto
CENTER(auto a, auto b, int center = 0) {
    return ((abs(a - center) < abs(b - center)) ? a : b);
}

auto 
MIN(auto a, auto b) {
    return (a < b ? a : b);
}

auto
MIN3(auto a, auto b, auto c) {
    return MIN(a, MIN(b, c));
}

auto
MIN4(auto a, auto b, auto c, auto d) {
    return MIN(a, MIN3(b, c, d));
}

auto 
MAX(auto a, auto b) {
    return (a > b ? a :b);
}

auto 
MAX3(auto a, auto b, auto c) {
    return MAX(a, MAX(b, c));
}

auto 
MAX4(auto a, auto b, auto c, auto d) {
    return MAX(a, MAX3(b, c, d));
}

auto 
MAX5(auto a, auto b, auto c, auto d, auto e) {
    return MAX(a, MAX4(b, c, d, e));
}

auto 
MAX6(auto a, auto b, auto c, auto d, auto e, auto f) {
    return MAX(a, MAX5(b, c, d, e, f));
}

auto
BETWEEN(auto v, auto a, auto b) {
    return MIN(MAX(a, v), b);
}

auto add01_between_min01_max01(auto a, auto b) {
    return BETWEEN(MIN(a, b), a + b, MAX(a ,b));
}

void 
UPDATE_MIN(auto & a, const auto & b) {
    a = MIN(a, b);
}

void 
UPDATE_MAX(auto & a, const auto & b) {
    a = MAX(a, b);
}

auto 
SUM2(const auto & vec) {
    static_assert (vec.size() == 2);
    return vec[0] + vec[1];
}

auto
SUMVEC(const auto & vec) {
    auto r = 0;
    for (unsigned int i = 0; i < vec.size(); i++) {
        r += vec[i];
    }
    return r;
}

void 
UPDATE_MIN2(auto & a, const auto & b) {
    for (int i = 0; i < 2; i++) { UPDATE_MIN(a[i], b[i]); }
}

void 
UPDATE_MAX2(auto & a, const auto & b) {
    for (int i = 0; i < 2; i++) { UPDATE_MAX(a[i], b[i]); }
}

template <bool TIsSorted = false>
auto 
MEDIAN(auto v) {
    if (!TIsSorted) {
        std::sort(v.begin(), v.end());
    }
    if (v.size() % 2 == 1) {
        return v[v.size() / 2];
    } else {
        return (v[v.size() / 2] + v[v.size() / 2 - 1]) / 2;  
    }
}

template <class V, class W>
V
hmean(V v1, W w1, V v2, W w2) {
    const auto f = (1024UL * 64UL);
    return  f * (w1 + w2) / (f * w1 / (v1+1UL) + f * w2 / (v2+1UL) + 1UL);
}

template <class T>
T
calc_non_negative(const T v, T base = pow(10.0, 0.1), T thres = 20.0) {
    if (v < thres) {
        return log1p(pow(base, v)) / log(base);
    } else {
        return v;
    }
}

auto 
safediv0(auto a, auto b) {
    return (b != 0 ? a / b : 0);
}

auto 
mathsquare(auto x) { 
    return x * x; 
}

const double
prob2odds(double p) {
    assert(0.0 < p && p < 1.0);
    return p /  (1.0 - p);
}

const double
logit(double p) {
    return log(prob2odds(p));
}

const double
logit2(double a, double b) {
    return logit(a/(a+b));
}

template <bool TIsBiDirectional = false>
constexpr double
calc_binom_10log10_likeratio(double prob, double a, double b) {
    assert(prob > 0 && prob < 1);
    a += DBL_EPSILON;
    b += DBL_EPSILON;
    double A = (      prob) * (a + b);
    double B = (1.0 - prob) * (a + b);
    if (TIsBiDirectional || a > A) {
        // likeratio = pow(a / A, a) * pow(b / B, b);
        return 10.0 / log(10.0) * (a * log(a / A) + b * log(b / B));
    } else {
        return 0.0;
    }
}

// h0 and h1 are two explanations, while h2 is responsible for the absence of h1 
// (e.g., h0 is tumor var quality, h1 is germline polymorphism, and h2 is contamination)
// h0 and h1 are two possible causes of the outcome, while h2 inhibits h1
template <class T>
constexpr T
max_min01_sub02(T h0, T h1, T h2) {
    return MAX(MIN(h0, h1), h0 - h2);
}

#ifdef TEST_calc_binom
int main(int argc, char **argv) {
    double prob = atof(argv[1]);
    double a = atof(argv[2]);
    double b = atof(argv[3]);
    printf("calc_binom_10log10_likeratio(%f, %f, %f) = %f\n", prob, a, b, calc_binom_10log10_likeratio(prob, a, b));
}
#endif

static_assert(abs(calc_binom_10log10_likeratio(0.1, 10, 90)) < 1e-4);
static_assert(calc_binom_10log10_likeratio(0.1, 90, 10) > 763); // 10/log(10) * (90*log(9)+10*log(1/9))
static_assert(calc_binom_10log10_likeratio(0.1, 90, 10) < 764); // 10/log(10) * (90*log(9)+10*log(1/9))
static_assert(abs(calc_binom_10log10_likeratio(0.1, 1, 99)) < 1e-4); // 10/log(10) * (90*log(9)+10*log(1/9))

template <class T> 
void 
autoswap ( T& a, T& b ) {
    T c(a); a=b; b=c;
}

double invmax(double x) { return MAX(x, 1/x); }

template <bool TRefExcludeAlt=false, bool TIsPseudocountZero = false, class T>
double
_any4_to_biasfact(T dp0, T dp1, T ad0, T ad1, const bool is_inv, double pseudocount) {
    if (TRefExcludeAlt) {
        return _any4_to_biasfact<false>(dp0 + ad0, dp1 + ad1, ad0, ad1, is_inv, pseudocount);
    }
    if (is_inv) {
        // These assertions now result in errors when enabled due to extended possible usage
        //assert(dp0 <= ad0 || !(std::cerr << dp0 << " <= " << ad0 << " failed!" << std::endl));
        //assert(dp1 <= ad1 || !(std::cerr << dp1 << " <= " << ad1 << " failed!" << std::endl));
        assert(dp0 >= -pseudocount/2);
        assert(dp1 >= -pseudocount/2);
    } else {
        assert(dp0 >= ad0);
        assert(dp1 >= ad1);
        assert(ad0 >= -pseudocount/2);
        assert(ad1 >= -pseudocount/2);
    }
    double t00 = (TIsPseudocountZero ? dp0 : (dp0 + pseudocount));
    double t01 = (TIsPseudocountZero ? dp1 : (dp1 + pseudocount));
    double t10 = (TIsPseudocountZero ? ad0 : (ad0 + pseudocount));
    double t11 = (TIsPseudocountZero ? ad1 : (ad1 + pseudocount));
    
    double t0sum = t00 + t01;
    double t1sum = t10 + t11;
    
    if (!TIsPseudocountZero) { 
        if (t0sum > t1sum) {
            t00 += pseudocount * (t0sum/t1sum - 1);
            t01 += pseudocount * (t0sum/t1sum - 1);
        } else {
            t10 += pseudocount * (t1sum/t0sum - 1);
            t11 += pseudocount * (t1sum/t0sum - 1);
        }
    }
    double t00f = t00 / t01; // observed 
    double t10f = t10 / t11; // expected 
    double weight1over01 = t11 / (t10 + t11);
    
    double raw_biasfact = (t00f - t10f) * (weight1over01 / t10f) + 1;
    return raw_biasfact;    
}

template <bool TRefExcludeAlt=false, class T>
unsigned int
any4_to_biasfact100(T dp0, T dp1, T ad0, T ad1, const bool is_inv = false, double pseudocount = 1) {
    return floor(_any4_to_biasfact(dp0, dp1, ad0, ad1, is_inv, pseudocount) * (100 * (1 + DBL_EPSILON)));
}

double
biasfact100_to_imba(unsigned int biasfact100) {
    return MAX(SIGN2UNSIGN(100), biasfact100) / (double)100;
}

#ifdef TEST_any4_to_bias
int 
main(int argc, char **argv) {
    double t00 = atof(argv[1]);
    double t01 = atof(argv[2]);
    double t10 = atof(argv[3]);
    double t11 = atof(argv[4]);
    double p = atof(argv[5]);

    auto bf1 = any4_to_biasfact100<false>(t00, t01, t10, t11, false, p);
    auto bf2 = any4_to_biasfact100<false>(t01, t00, t11, t10, false, p);
    printf("any4_to_bias100(%f %f %f %f %f) == %d\n", t00, t01, t10, t11, p, bf1);
    printf("any4_to_bias100(%f %f %f %f %f) == %d\n", t01, t00, t11, t10, p, bf2);
}

#endif

double 
geometric_sum_to_nterms(double geosum, double term1, double ratio) {
    double ret = log(geosum * (ratio - 1) / term1 + 1) / log(ratio);
#ifdef TEST_h01_to_phredlike
    printf("geometric_sum_to_nterms(%f, %f, %f) = %f\n", geosum, term1, ratio, ret);
#endif
    return ret;
}

double 
dlog(double n, double r) {
    return log(n * (r-1) + 1) / log(r);
}

template <bool TIsConsensual, bool TIsPseudocountZero = false, bool TIsErrAmpRatioOne = false> 
double
h01_to_phredlike(double h0pos, double h0tot, double h1pos, double h1tot, 
        double pseudocount, double err_amp_ratio) {
    assert(h0pos <  h0tot || !fprintf(stderr, "%lf <  %lf failed", h0pos, h0tot));
    assert(h1pos <= h1tot || !fprintf(stderr, "%lf <= %lf failed", h1pos, h1tot));
    assert(h0pos >  0     || !fprintf(stderr, "%lf >  %lf failed", h0pos, (double)0));
    assert(h1pos >= 0     || !fprintf(stderr, "%lf >= %lf failed", h1pos, (double)0));
    if (h1pos <= (DBL_EPSILON * pseudocount)) {
        return 0;
    }
    double ratio = (h1tot) / h0tot;
    h0tot *= ratio;
    h0pos *= ratio;
    double bf;
    if (TIsPseudocountZero) {
        bf = _any4_to_biasfact<false, true>(h0tot, h1tot, h0pos, h1pos, true, 0);
    } else {
        bf = _any4_to_biasfact(h0tot + pseudocount, h1tot + pseudocount, h0pos + pseudocount, h1pos, true, 0);
    }
    double dlogval;
    if (TIsErrAmpRatioOne) {
        dlogval = bf;
    } else {
        dlogval = dlog(bf, err_amp_ratio);
    }
    return dlogval * log((h1pos * h0tot) / (h1tot * h0pos)) * (10.0 / log(10.0));
}

struct Any4Value {
    const double v1;
    const double v2;
    const double v3;
    const double v4;
    Any4Value(double v1, double v2, double v3, double v4) : v1(v1), v2(v2), v3(v3), v4(v4) {}
    const double to_phredlike(double d, double pc = 0, double err_amp_ratio = 1+1e-5) const {
        return (d + h01_to_phredlike<false, true, true>(v1 + d, (v2 + d) * (1 + DBL_EPSILON), v3 + d, (v4 + d) * (1 + DBL_EPSILON), pc, err_amp_ratio));
    }
};

// retrieved from https://en.wikipedia.org/wiki/Golden-section_search
Any4Value
gss(const auto & any4Value, double a, double b, double tol = 0.5) {
    Any4Value ret(0, 0, 0, 0);
    double invphi = (sqrt(5) - 1) / 2;
    double invphi2 = (3 - sqrt(5)) / 2;
    double h = b - a;
    assert(h > tol);
    /*
    if (h <= tol) {
        return Any4Value(a, b, 0, 0);
    }
    */
    int n = (int)(ceil(log(tol/h)/log(invphi)));
    double c = a + invphi2 * h;
    double d = a + invphi * h;
    double yc = any4Value.to_phredlike(c);
    double yd = any4Value.to_phredlike(d);
    for (int k = 0; k < n - 1; k++) {
        if (yc < yd) {
            b = d;
            d = c;
            yd = yc;
            h = invphi*h;
            c = a + invphi2 * h;
            yc = any4Value.to_phredlike(c);
        } else {
            a = c;
            c = d;
            yc = yd;
            h = invphi*h;
            d = a + invphi * h;
            yd = any4Value.to_phredlike(d);
        }
    }
    if (yc < yd) {
        double ya = any4Value.to_phredlike(a);
        double yb = any4Value.to_phredlike(b);
        return Any4Value(a, b, ya, yb);
    } else {
        return Any4Value(c, d, yc, yd);
    }
}

double
sumBQ4_to_phredlike(double & bestAddValue, 
        double normalAllBQsum, double normalAltBQsum, double tumorAllBQsum, double tumorAltBQsum) {
    Any4Value any4Value(normalAltBQsum, normalAllBQsum, tumorAltBQsum, tumorAllBQsum);
    Any4Value argmin2_min2 = gss(any4Value, 1, 3 + tumorAltBQsum);
    bestAddValue = argmin2_min2.v2;
    return argmin2_min2.v4;
}

// This function is here for regression test
// one-way conversion of information into other measures of information (based on information theory)
// consensual means in the same MIG, 
// homogeneity=0 means reads are completely independent (e.g. different DNA fragments), homogeneity=1 means reads are completely dependent (e.g. in the same MIG)
template <bool TIsConsensual = true> 
double
_old_h01_to_phredlike(double h0pos, double h0tot, double h1pos, double h1tot, 
        double pseudocount = 1, double homogeneity = (TIsConsensual ? 1 : 0), double err_amp_ratio = (TIsConsensual ? 2 : pow(2, 0.5))) {
    assert(h0pos <  h0tot || !fprintf(stderr, "%lf <  %lf failed", h0pos, h0tot));
    assert(h1pos <= h1tot || !fprintf(stderr, "%lf <= %lf failed", h1pos, h1tot));
    assert(h0pos >  0     || !fprintf(stderr, "%lf >  %lf failed", h0pos, 0));
    assert(h1pos >= 0     || !fprintf(stderr, "%lf >= %lf failed", h1pos, 0));
    
    double h1neg = h1tot - h1pos;
    double h0freq = ((double)(h0pos)) / ((double)(h0tot));
    double h1freq = ((double)(h1pos)) / ((double)(h1tot));
    h0freq = MIN(h1freq, h0freq);
    double kldiv;
    if (TIsConsensual) {
        kldiv = ((h1pos <= 0 ? 0 : h1pos * log(h1freq / h0freq)) + (h1neg <= 0 ? 0 : h1neg * log((1.0 - h1freq) / (1.0 - h0freq)))) / h1tot + DBL_EPSILON;
    } else {
        kldiv = ((h1pos <= 0 ? 0 : log(h1freq / h0freq))) + DBL_EPSILON;
    }
    assert(kldiv >= 0 || !fprintf(stderr, "kldiv value of %lf is found for (%lf , %lf , %lf , %lf ) and TIsConsensual=%d\n", kldiv, h0pos, h0tot, h1pos, h1tot, TIsConsensual));
    // a * (1 -r^n) / (1-r) = m, n = math.log(m*(r-1)/a + 1, r)
    double nreads = geometric_sum_to_nterms(h1pos, h0freq * h1tot + pseudocount, err_amp_ratio);
    return (kldiv) * nreads * (10 / log(10));
}

#ifdef TEST_h01_to_phredlike
int 
main(int argc, char **argv) {
    double result1 = _old_h01_to_phredlike<true>(10, 30000, 20, 30);
    double result2 = h01_to_phredlike<true>(10, 30000, 20, 30, 1, 1.5);
    printf("result12 = %f %f \n", result1, result2);
    if (argc > 4) {
        double a1 = atof(argv[1]);
        double a2 = atof(argv[2]);
        double b1 = atof(argv[3]); 
        double b2 = atof(argv[4]);
        double pc = atof(argv[5]);
        double ear = atof(argv[6]);
        double result1 = _old_h01_to_phredlike<false>(a1, a2, b1, b2, pc);
        double result2 = h01_to_phredlike<false>(a1, a2, b1, b2, pc, ear);
        printf("result12 = %f %f \n", result1, result2);
    }
}
#endif

double 
dp4_to_sratio(double all_fw0, double all_rv0, double alt_fw0, double alt_rv0, double pseudocount = 1) {
    assert(all_fw0 >= 0 - FLT_EPSILON);
    assert(all_rv0 >= 0 - FLT_EPSILON);
    assert(alt_fw0 >= 0 - FLT_EPSILON);
    assert(alt_rv0 >= 0 - FLT_EPSILON);
    double all_fw = all_fw0 + pseudocount;
    double all_rv = all_rv0 + pseudocount;
    double alt_fw = alt_fw0 + pseudocount;
    double alt_rv = alt_rv0 + pseudocount;
    double rawratio = (all_fw * alt_rv) / (all_rv * alt_fw);
    double sumratio = rawratio + 1.0 / rawratio;
    double allratio = MIN(all_fw, all_rv) / MAX(all_fw, all_rv);
    double altratio = MIN(alt_fw, alt_rv) / MAX(alt_fw, alt_rv);
    double ret = sumratio * allratio / altratio;
    assert(0 < ret);
    return ret;
}

// conversion between different defintions in bioinformatics

const unsigned int 
char2phred(const unsigned char charvalue) {
    return charvalue - 33;
}

const unsigned char 
phred2char(const unsigned int phredvalue) {
    return phredvalue + 33;
}

const double 
phred2prob(const unsigned int phredvalue) {
    return pow(10, -((float)phredvalue) / 10);
}

const unsigned int 
prob2phred(const double probvalue) {
    return floor(-10 * log(probvalue) / log(10));
}

const unsigned int 
phred2bucket(const unsigned int phredvalue) {
    assert(phredvalue < NUM_BUCKETS);
    return phredvalue;
    // With MIN(32-1, phredvalue / 2) we may some some resolution
    // 0 - 8 -> 1, 8 - 40 -> 2, 40 - .. -> 4 ; 
    //  0 - 10 : 1 -> 10
    // 10 - 20 : 2 -> 5
    // 20 - 30 : 1 -> 10
    // 30 - 40 : 2 -> 5
    // 40 - 60 : 4 -> 5
    // Another binning schme is : floor(20 * log(phredvalue + 1) / log(10))
}

const unsigned int 
bucket2phred(const unsigned int bucketvalue) {
    return bucketvalue;
}

struct _PhredToErrorProbability {
    const unsigned int CONST2POW16 = 256*256;
    double data[128];
    uint32_t over2pow16[128];
    _PhredToErrorProbability() {
        for (int i = 0; i < 128; i++) {
            this->data[i] = phred2prob(i);
            this->over2pow16[i] = ceil(this->data[i] * CONST2POW16);
        }
    }
};

const _PhredToErrorProbability THE_PHRED_TO_ERROR_PROBABILITY;

#endif

