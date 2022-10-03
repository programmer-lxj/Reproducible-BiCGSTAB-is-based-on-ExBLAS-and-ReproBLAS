
#ifndef COMMON_H
#define COMMON_H

#include <math.h>

double norm_inf(int n_dist, double *res_err) {
    double nrm = 0.0;

    for (int i = 0; i < n_dist; i++) {
        double tmp = fabs(res_err[i]);
        if (nrm < tmp)
            nrm = tmp;
    }

    return nrm;
}

#endif // COMMON_H
