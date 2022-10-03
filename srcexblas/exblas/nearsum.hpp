
// Courtesy of the OzBLAS by Daichi Mukunoki and Takeshi Ogita
//
// =========================================
// NearSum (and AccSum)
// from "ACCURATE FLOATING-POINT SUMMATION" by S.M.RUMP, T.OGITA, S.OISHI (2005)
// http://www.ti3.tu-harburg.de/paper/rump/RuOgOi06.pdf
// =========================================

#define MAX(a, b) ((a) > (b) ? (a) : (b))

__inline__
void Transform (const int32_t n, double *vec, const int32_t ld, double rho, double &tau1, double &tau2) {
	int32_t i, m;
	double tmp, mu, sigma, t, tau;

	mu = fabs(vec[0]);
	for (i = 1; i < n; i++)
		mu = MAX (mu, fabs(vec[i*ld])); 
	if ((n == 0) || (mu == 0.)) {
		tau1 = rho;
		tau2 = 0.;
		return;
	}

	m = ceil (log2((double)(n+2)));
	sigma = scalbn (1., m+ceil(log2(mu)));
	t = rho;
	while (1) {
		// ExtractVector
		tau = 0.;
		for (i = 0; i < n; i++) {
			tmp = (sigma + vec[i*ld]) - sigma;
			vec[i*ld] -= tmp;  // <- output
			tau += tmp;
		}
		// here, tau1 = t1
		tau1 = t + tau;
		if ((sigma <= scalbn (1., -1022)) || (fabs(tau1) >= scalbn (1., 2*m+1-53)*sigma)) {
			//FastTwoSum (t, tau, &tau1, &tau2);
			//tau1 = t + tau
			tau2 = tau - (tau1 - t);
			return;
		}
		sigma = scalbn (1., m-53) * sigma;
		t = tau1;
	} 
}

__inline__
void TransformK (const int32_t n, double *vec, const int32_t ld, double rho, double &res, double &r) { 
	int32_t i;
	double tmp, tau1, tau2;
	Transform (n, vec, ld, rho, tau1, tau2);
	tmp = 0.;
	for (i = 0; i < n; i++)
		tmp += vec[i*ld];
	res = tau1 + (tau2 + tmp);
	r = tau2 - (res - tau1);
}

__inline__
double Sign (double v) {
	return (v < 0) ? -1.:1.;
}

double NearSum (const int32_t n, double *vec, const int32_t ld) {
	double tmp, res, res2, r, r2, mu, delta, delta2;
	double eps = scalbn (1., -53);

	TransformK (n, vec, ld, 0, res, r);
	TransformK (n, vec, ld, r, delta, r2);
	if (delta == 0) 
		return res;
	res2 = res + Sign (delta) * eps * fabs(res);
	if (res2 == res) {
		mu = Sign (delta) * eps * fabs(res);
		res2 = res + 2. * Sign (delta) * eps * fabs(res);
	} else {
		mu = (res2 - res) / 2.;
	}
	if (fabs(delta) < fabs(mu)) 
		return res;
	if (fabs(delta) > fabs(mu)) 
		return res2;
	TransformK (n, vec, ld, r2, delta2, tmp);
	if (delta2 == 0) 
		return res + mu;
	if (Sign (delta2) == Sign (mu))
		return res2;
	return res;
}
