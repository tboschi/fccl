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
#include <random>

#include "fccl/point.h"
#include "fccl/rate.h"
#include "fccl/poisson.h"
#include "fccl/belt.h"

template <size_t N>
class region {
	public:
		region() : gen(std::random_device()()) { };
		template <typename RNG>
		region(RNG rng) : gen(rng) { };
		//region() = default;
		//region(const point &pp) : bb(belt<N>(pp)) {};

		// set the start point in a variadic way
		// this will overwrite the existing belt
		//template<typename... Args>
		//void start(Args... args) {
		//	bb = belt<N>(args...);
		//}

		//void start(const point &pp) {
		//	bb = belt<N>(pp);
		//}

		// routine will exapand belt until CL is reached
		template<typename... Args, typename F = double>
		belt<N> expand(F CL, Args... args) {
			std::vector<rate> rr { args... };
			return expand(CL, rr);
		}

		template<typename F = double>
		belt<N> expand(F CL, const std::vector<rate> &rr) {
			if (rr.size() != N)
				throw std::invalid_argument("rates not compatible with this region");
			if (CL >= 1.0)
				throw std::invalid_argument("CL cannot be greater than 1");

			// initialize start point
			point start(rr.size());
			std::transform(rr.begin(), rr.end(), start.begin(),
					[](rate r0) { return size_t(std::round(r0.sig + r0.bak)); } );
			belt<N> bb(start);

			double best = prob(start, rr);
			double sum = 1.;	 // only sum after the peak was reached

			while (sum < CL / best) {

				// next point with minimum likelihood
				std::vector<point> points = bb.closest();
#ifdef DEBUG
				std::cout << "closest points\n";
				for (const auto &pp : points)
					std::cout << "\t" << pp;
				std::cout << "\n\n";
#endif
				point pp = next(points.begin(), points.end(), rr, start);

				sum += part(start, pp, rr);

				bb.add(pp);
			}

			return bb;
		}

		// return a vector of points to break any tie wit random choice
		point next(std::vector<point>::iterator first, std::vector<point>::iterator last,
			   const std::vector<rate> &rr, const point &ss) {
			if (rr.size() != N)
				throw std::invalid_argument("rates not compatible with this region");
			std::vector<point>::iterator imin = first;
			std::vector<point> pmins = {*first};

			double min_val = hood(*first++, rr);
			for (; first != last; ++first) {
				double val = hood(*first, rr);
				if (val < min_val) {
					min_val = val;
					pmins = {*first};
				}
				else if (val == min_val)
					pmins.push_back(*first);
			}

			if (pmins.size() == 1)
				return pmins[0];

			std::uniform_int_distribution<int> dis(0, pmins.size()-1);
			int pick = dis(gen);
			return pmins[dis(gen)];
		}

		// computes poissonian probability of the full point
		double prob(const point &pp, const std::vector<rate> &rr) const {
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
		double part(const point &p0, const point &pp, const std::vector<rate> &rr) const {
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
		double hood(const point &pp, const std::vector<rate> &rr) const {
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
		std::mt19937 gen;
};

#endif
