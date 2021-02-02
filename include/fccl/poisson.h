/* This namespace defines some functions
 * relative to 1D poissonian distributions
 *
 * Multi-D equivalent definitions are implementation specific
 */

#ifndef POISSON_H
#define POISSON_H

#include <cmath>

#include "fccl/rate.h"

namespace poisson {
	// n observed events
	// s is signal
	// b is background
	// function above has minimum value for n = s + b
	inline double llratio(double n, double s, double b) {
		// correct also if any term is zero
		return 2 * (s - std::max(n - b, 0.))
		     - 2 * (n > 0 ? n * (std::log(s + b)
				       - std::log(std::max(n, b))) : 0);
	}

	inline double llratio(double n, const rate &r) {
		return llratio(n, r.sig, r.bak);
	}

	inline double poisson(double n, double s) {
		if (n > 0) {
			double p = std::exp(- s / n) * s;
			long double ret = p / n;
			for (--n; n > 0; --n) {
				ret *= p / n;
			}

			return ret;
		}

		return std::exp(-s);
	}

	inline double poisson(double n, const rate &r) {
		return poisson(n, r.sig + r.bak);
	}

	/* poisson is exp(-s) * pow(s, n) / n!
	 * CL is reached by sum_n exp(-s) * pow(s, n) / n!
	 * exp(-s) * sum[ pow(s, n) / n! ]
	 * and n contiguos in a range around n0
	 * exp(-s) * sum[ pow(s, n0 + n - n0) / n0! * n0! / n! ]
	 * exp(-s) * pow(s, n0) / n0! * sum[ pow(s, n - n0) * n0! / n! ]
	*/

	// function computes a factorized term in the poissonian sum
	// partial is basicall poisson(n, r) / poisson(n0, r)
	inline double partial(double n0, double n, double s) {
		// n > n0 correct order
		long double ret = 1.;
		for (double k = std::max(n, n0); k > std::min(n, n0); --k)
			ret *= s / k;

		return n >= n0 ? ret : 1./ret;
	}

	inline double partial(double n0, double n, const rate &r) {
		return partial(n0, n, r.sig + r.bak);
	}
};

#endif
