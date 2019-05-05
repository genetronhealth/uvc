#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double lbeta(double a, double b) {
    return lgamma(a) + lgamma(b) - lgamma(a+b);
}

double lnchooser(double n, double r) {
    return lgamma(n+1) - lgamma(n-r+1) - lgamma(r+1);
}

double loglike(double count, double total) {
    return (count <= 0 ? 0 : count * log(count/total));
}

double lratio(double count, double total) {
    assert(count >= 0);
    assert(count <= total);
    return loglike(count, total) + loglike(total-count, total);
}

double lbbratio(double ad, double dp, double a, double b) {
    double h0loglike = lnchooser(dp, ad) + lbeta(ad+a ,   (dp-ad)+b) - lbeta(a , b); // lbeta(ad+a, dp-ad+b) - lbeta(a, b);
    double h1loglike = log(1/(dp+1)) ; // lbeta(ad+1,    (dp-ad)+1) - lbeta(1, 1); // lratio(ad, dp) ; //lbeta(ad+ad, dp) - lbeta(ad, ad);
    double h0loglike1 = lratio(a/(a+b)*dp, dp);
    fprintf(stderr, "h0： %f ; h1: %f ; h2: %f\n", h0loglike, h1loglike, h0loglike1);
    return h1loglike - h0loglike;
}

double lbratio(double ad, double dp, double a, double b) {
    double h0loglike = loglike(a, b); // lbeta(ad+a, dp-ad+b) - lbeta(a, b);
    double h1loglike = loglike(ad, dp) ; //lbeta(ad+ad, dp) - lbeta(ad, ad);
    //fprintf(stderr, "%f - %f\n", h0loglike, h1loglike);
    return h1loglike - h0loglike;
}

double nat2phred(double nat) {
    return 10/log(10) * nat;
}

int main(int argc, char **argv) {
    assert(argc > 4);
    int lbbtest = !strcmp("lbbratio", argv[1]);
    int lbtest  = !strcmp("lbratio",  argv[1]); 
    if (lbbtest || lbtest) {
        double result = 0;
        int niters = atoi(argv[2]);
        double eps = atof(argv[3]);
        double mult = atof(argv[4]); 
        for (int i = 0; i < niters; i++) {
            double mult2 = mult * (double) i;
            if (lbbtest) {
                result += lbbratio(1 + eps*mult, 3000 + eps*mult, 1 + eps*mult, 100*1000 + eps*mult);
            } else if (lbtest) {
                result += lbratio(1 + eps*mult, 3000 + eps*mult, 1 + eps*mult, 100*1000 + eps*mult); 
            }
        }
        printf("dummy-result-for-%d-iters = %f\n", niters, result);
    } else {
        double ad = atof(argv[1]);
        double dp = atof(argv[2]);
        double err_alpha = atof(argv[3]);
        double err_beta  = atof(argv[4]);
        double nat = lbbratio(ad, dp, err_alpha, err_beta);
        printf("%f\t%f\t%f\t%f\t%f\n", ad, dp, err_alpha, err_beta, nat2phred(nat));
    }
}

