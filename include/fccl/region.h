/* Creates a confidence belt in a <N> dimensional space
 * given expected signal and expected background for 
 * each dimension as input.
 * A dimension can be a detection channel for example.
 * The channels are all uncorrelated with each other.
 *
 * The belt lives on a grid of integer values in N dimensions
 * and the points encompassed by the belt are determined using
 * the Feldman Cousin principle for a given confidence level (CL),
 * i.e. orderered by minimum log-likelihood ratio such that
 * the sum of probabilities inside the region is equal to the CL
 */

#ifndef REGION_H
#define REGION_H

#include <vector>

#include "fccl/point.h"
#include "fccl/rate.h"
#include "fccl/poisson.h"
#include "fccl/belt.h"

template <size_t N>
class region {
	public:
		region() = default;
		region(const point &pp) : bb(belt<N>(pp)) {};

		// set the start point in a variadic way
		// this will overwrite the existing belt
		template<typename... Args>
		void start(Args... args) {
			bb = belt<N>(args...);
		}

		void start(const point &pp) {
			bb = belt<N>(pp);
		}

		// routine will exapand belt until CL is reached
		template<typename... Args, typename F = double>
		belt<N> & expand(F CL, Args... args) {
			std::vector<rate> rr { { args...} };
			return expand(CL, rr);
		}

		template<typename F = double>
		belt<N> & expand(F CL, const std::vector<rate> &rr) {
			auto points = bb.points();
			if (points.size() == 0) // empty belt
				throw std::logic_error("belt is not initialized with a start point");
			if (points.size() > 1) // belt not fresh
				throw std::logic_error("belt does not contain only the start point");

			double best = prob(points.front(), rr);
			double sum = 1.;	 // only sum after the peak was reached

			while (sum < CL / best) {

				// next point with minimum likelihood
				point pp = next(rr);
				sum += part(points.front(), pp, rr);
			}

			return bb;
		}

		point next(const std::vector<rate> &rr) {
			if (rr.size() != N)
				throw std::invalid_argument("rates not compatible with this region");
			std::vector<point> points = bb.closest();
			std::vector<point>::iterator ip = points.begin(), imin = ip;

			double min_val = hood(*ip, rr);
			++ip;
			for (; ip != points.end(); ++ip) {
				if (hood(*ip, rr) < min_val) {
					min_val = hood(*ip, rr);
					imin = ip;
				}
			}

			// expand belt
			bb.add(*imin);

			return *imin;
		}

		/*
		std::vector<double> reject(const std::vector<double> &bak) {
			if (bak.size() != N)
				throw std::invalid_argument("background not compatible with this region");
			std::vector<rate> rr;
			rr.reserve(bak.size());
			for (const auto & b : bak)
				rr.emplace_back(sqrt(b), b);
			// starting point
			point sp(bak.begin(), bak.end());
			start(sp);	// <- set starting point
			while (bb.contains(sp)) {
				sig += 0.01;
				point pp = {size_t(std::round(bak + sig))};
				rate rr(sig, bak);

				reg.start(pp);
				myb = reg.expand(0.99, rr);

				// when background is not in belt anymore
				// this loop will end
			}

		}
		*/

		// computes poissonian probability of the full point
		double prob(const point &pp, const std::vector<rate> &rr) {
			if (rr.size() != pp.size())
				throw std::invalid_argument("rates not compatible with this point");

			auto ip = pp.begin();
			auto ir = rr.begin();
			double ret = poisson::poisson(*ip, *ir);
			++ip; ++ir;
			for ( ; ip != pp.end(); ++ip, ++ir)
				ret *= poisson::poisson(*ip, *ir);

			return ret;
		}

		// computes partial poissonian probability of the full point
		double part(const point &p0, const point &pp, const std::vector<rate> &rr) {
			if (rr.size() != pp.size())
				throw std::invalid_argument("rates not compatible with this point");

			auto ip = pp.begin(), i0 = p0.begin();
			auto ir = rr.begin();
			double ret = poisson::partial(*i0, *ip, *ir);
			++i0; ++ip; ++ir;
			for ( ; ip != pp.end(); ++ip, ++ir)
				ret *= poisson::partial(*i0, *ip, *ir);

			return ret;
		}

		// computes the log-likelihood ratio of the full point
		double hood(const point &pp, const std::vector<rate> &rr) {
			if (rr.size() != pp.size())
				throw std::invalid_argument("rates not compatible with this point");

			auto ip = pp.begin();
			auto ir = rr.begin();
			double ret = poisson::llratio(*ip, *ir);
			++ip; ++ir;
			for ( ; ip != pp.end(); ++ip, ++ir)
				ret += poisson::llratio(*ip, *ir);

			return ret;
		}

	private:
		belt<N> bb;
};

#endif
