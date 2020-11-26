/* belt object describes a surface in a N-dimensional
 * discretized space. Instead of saving each point of
 * the surface individually, this structure saves on-
 * ly the minimum and maximum coordinates in each di-
 * mension. The memory needed is order
 *   2 * (1 + n * (1 + n * ( .. )))  =
 *   2 * (1 + n + n^2 + n^3 ... + n^{N-1})
 * which is conveniente especially for large dimensions
 *
 * The class uses metaprogramming to generate N nested
 * containers (deque) at compilation time.
 */

#ifndef BELT_H
#define BELT_H

#include <cmath>
#include <deque>
#include <vector>
#include <set>
#include <memory>
#include <numeric>
#include <algorithm>

#include "fccl/point.h"

// N D belt
template <int N, typename std::enable_if<(N > 0), bool>::type = true>
struct belt {
	belt() : _order(N) { };

	template<typename T, typename... Args>
		//typename std::enable_if<std::is_fundamental<T>::value>::type = true>
	belt(T t, Args... args) : _order(N), nn{t, t} {
		static_assert(sizeof ...(args) == N-1, "incorrect number of arguments");
		up.push_back(belt<N-1>(args...));
	}

	//template<>
	belt(const point &pp) : _order(N) {
		if (pp.size() != _order)
			throw std::invalid_argument("point not compatible with this belt");
		nn = {pp.front(), pp.front()};
		up.push_back(belt<N-1>(point(pp.begin()+1, pp.end())));
	}

	// returns number of 1D belts in memory
	size_t memory() const {
		return 1 + std::accumulate(up.begin(), up.end(), 0,
				[](size_t sz, const belt<N-1> &b) {
					return sz += b.memory(); } );
	}

	// returns size of belt, i.e. number of points
	size_t size() const {
		return  std::accumulate(up.begin(), up.end(), 0,
				[](size_t sz, const belt<N-1> &b) {
					return sz += b.size(); } );
	}

	bool contains(const point &pp) const {
		if (pp.size() != N)
			throw std::invalid_argument("point not compatible with this belt");
		if (nn[0] > pp[0] || pp[0] > nn[1])
			return false;

		// this should be a valid entry
		return up[pp[0] - nn[0]].contains(point(pp.begin()+1, pp.end()));
	}

	// true if this belt and the input ones share a point
	bool share(const belt &bb) const {
		if (_order != bb._order)
			throw std::invalid_argument("comparing belts with different orders");

		for (size_t at = std::max(nn[0], bb.nn[0]);
			    at < std::min(nn[1], bb.nn[1]); ++at)
			if (up[at-nn[0]].share(bb.up[at-nn[1]]))
				return true;

		return false;
	}

	void add(const point &pp) {
		if (pp.size() != N)
			throw std::invalid_argument("point not compatible with this belt");
		if (contains(pp)) // no op
			return;

		if (pp[0] == nn[0] - 1) {
			--nn[0];
			up.push_front(point(pp.begin()+1, pp.end()));
		}
		else if (pp[0] == nn[1] + 1) {
			++nn[1];
			up.push_back(point(pp.begin()+1, pp.end()));
		}
		else
			throw std::range_error("only adjacient points can be added");
	}

	// return a vector the points delineating this belt
	std::vector<point> points() const {
		std::vector<point> this_cord;
		// by implementation nn[1] >= nn[0]
		this_cord.reserve(2 * (nn[1] - nn[0]) + 1);
		size_t at = nn[0];
		// up size should be exactly as nn[1] - nn[0] + 1
		for (const auto &ib : up) {
			std::vector<point> up_cord = ib.points();
			for (auto &pp : up_cord) {
				pp.push_front(at);
				this_cord.push_back(std::move(pp));
			}
			++at;
		}

		this_cord.shrink_to_fit();
		return this_cord;
	}

	// large definition
	std::vector<point> neighbors() const {
		std::set<point> this_cord;
		std::vector<point> all_cord = points();
		
		size_t pos = (nn[0] > 0 ? nn[0] - 1 : nn[0]), end = nn[0]+1;
		for (const auto &ib : up) {
			std::vector<point> up_cord = ib.points();
			std::vector<point> up_next = ib.neighbors();
			up_next.insert(up_next.end(), up_cord.begin(), up_cord.end());
			for (size_t at = pos; at <= end; ++at) {
				for (auto pp : up_next) {
					pp.push_front(at);
					// do not add belt points
					if (std::find(all_cord.begin(), all_cord.end(), pp)
							== all_cord.end()) {
						std::cout << pp << std::endl;
						this_cord.insert(pp);
					}
				}
			}
			++end;
			pos = end - 2;
		}

		return std::vector<point>(this_cord.begin(), this_cord.end());
	}

	// smallest definition
	std::vector<point> closest() const {
		std::vector<point> this_cord;
		// by implementation nn[1] >= nn[0]
		this_cord.reserve(9 * (nn[1] - nn[0]) - 1);

		std::vector<point> up_cord;
		// first look before first point
		if (nn[0] > 0) {
			size_t at = nn[0] - 1;
			up_cord = up.front().points();
			for (auto &pp : up_cord) {
				pp.push_front(at);
				this_cord.push_back(std::move(pp));
			}
		}

		// next do middle points
		size_t at = nn[0];
		for (const auto &ib : up) {
			std::vector<point> up_cord = ib.closest();
			for (auto &pp : up_cord) {
				pp.push_front(at);
				this_cord.push_back(std::move(pp));
			}
			++at;
		}

		// finally after last point
		up_cord = up.back().points();
		for (auto &pp : up_cord) {
			pp.push_front(at);
			this_cord.push_back(std::move(pp));
		}

		this_cord.shrink_to_fit();
		return this_cord;
	}

	private:
		std::array<size_t, 2> nn;
		std::deque<belt<N-1> > up;
		size_t _order;
};

// 1 D belt
template<>
struct belt<1> {
	belt() : _order(1) { };

	template<typename T>
	belt(T t) : _order(1), nn({t, t}) { };

	//template<>
	belt(const point &pp) : _order(1) {
		if (pp.size() != _order)
			throw std::invalid_argument("point not compatible with this belt");
		nn = {pp.front(), pp.front()};
	}

	size_t memory() const {
		return 1;
	}

	size_t size() const {
		return (nn[1] == nn[0]) ? 1 : 2;
	}

	bool contains(const point &pp) const {
		if (pp.size() != _order)
			throw std::invalid_argument("point not compatible with this belt");
		return nn[0] <= pp[0] && pp[0] <= nn[1];
	}

	void add(const point &pp) {
		if (pp.size() != _order)
			throw std::invalid_argument("point not compatible with this belt");
		if (contains(pp)) // no op
			return;

		if (pp[0] == nn[0] - 1)
			--nn[0];
		else if (pp[0] == nn[1] + 1)
			++nn[1];
		else
			throw std::range_error("only adjacient points can be added");
	}

	std::vector<point> points() const {
		if (nn[1] == nn[0])
			return {point({nn[0]})};
		else
			return {point({nn[0], nn[1]})};
	}

	std::vector<point> neighbors() const {
		return closest();
	}

	// only poisitive numbers are accepted
	std::vector<point> closest() const {
		if (nn[0] > 0)
			return {point({nn[0]-1}), point({nn[1]+1})};
		else
			return {point({nn[1]+1})};
	}


	bool share(const belt &bb) const {
		return !(nn[0] > bb.nn[1] || nn[1] < bb.nn[0]);
	}

	private:
		std::array<size_t, 2> nn;
		size_t _order;
};

#endif
