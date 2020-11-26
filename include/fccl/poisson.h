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
	template<typename S>
	double llratio(S i, double s, double b) {
		double n = static_cast<double>(i);
		// correct also if any term is zero
		return 2 * (s - std::max(n - b, 0.))
		     - 2 * (n > 0 ? n * (std::log(s + b)
				       - std::log(std::max(n, b))) : 0);
	}

	template<typename S>
	double llratio(S i, const rate &r) {
		return llratio(i, r.sig, r.bak);
	}

	template<typename S>
	double poisson(S i, double s) {
		if (i > 0) {
			double n = static_cast<double>(i);
			double p = std::exp(- s / n) * s;
			long double ret = p / n;
			for (--n; n > 0; --n) {
				ret *= p / n;
			}

			return ret;
		}

		return std::exp(-s);
	}

	template<typename S>
	double poisson(S i, const rate &r) {
		double n = static_cast<double>(i);
		return poisson(i, r.sig + r.bak);
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
	template <typename S>
	double partial(S i0, S i, double s) {
		// n > n0 correct order
		double n0 = static_cast<double>(i0);
		double n = static_cast<double>(i);
		long double ret = 1.;
		for (double k = std::max(n, n0); k > std::min(n, n0); --k)
			ret *= s / k;

		return n >= n0 ? ret : 1./ret;
	}

	template <typename S>
	double partial(S i0, S i, const rate &r) {
		return partial(i0, i, r.sig + r.bak);
	}

};

#endif
