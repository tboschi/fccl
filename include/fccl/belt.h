/* belt object describes a surface in a N-dimensional
 * discretized space. Instead of saving each point of
 * the surface individually, this structure saves on-
 * ly the minimum and maximum coordinates in each di-
 * mension. The memory needed is order
 *   2 * (1 + n * (1 + n * ( .. )))  =
 *   2 * (1 + n + n^2 + n^3 ... + n^{N-1})
 * which is convenient especially for large dimensions.
 * Everything else is lazily evaluated.
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

	/*    ************    
	 *    CONSTRUCTORS
	 *    ************    */

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



	/*    *****************    
	 *    MEMORY MANAGEMENT
	 *    *****************    */

	// number of 1D belts in memory
	size_t memory() const {
		return 1 + std::accumulate(up.begin(), up.end(), 0,
				[](size_t sz, const belt<N-1> &b) {
					return sz += b.memory(); } );
	}

	// size of belt, i.e. number of points delimiting belt
	size_t size() const {
		return  std::accumulate(up.begin(), up.end(), 0,
				[](size_t sz, const belt<N-1> &b) {
					return sz += b.size(); } );
	}

	// points contained by the belt
	size_t capacity() const {
		return  std::accumulate(up.begin(), up.end(), 0,
				[](size_t sz, const belt<N-1> &b) {
					return sz += b.capacity(); } );
	}


	/*    *********************    
	 *    POINT-WISE OPERATIONS
	 *    *********************    */


	bool contains(const point &pp) const {
		if (pp.size() != _order)
			throw std::invalid_argument("point not compatible with this belt");
		if (nn[0] > pp[0] || pp[0] > nn[1])
			return false;

		// this should be a valid entry
		return up[pp[0] - nn[0]].contains(point(pp.begin()+1, pp.end()));
	}

	// true if this belt and the input ones share a point
	bool share(const belt<N> &bb) const {
		if (_order != bb._order)
			throw std::invalid_argument("comparing belts with different orders");

		for (size_t at  = std::max(nn[0], bb.nn[0]);
			    at <= std::min(nn[1], bb.nn[1]); ++at) {
			if (up[at-nn[0]].share(bb.up[at-bb.nn[0]]))
				return true;
		}

		return false;
	}

	void add(const point &pp) {
		if (pp.size() != _order)
			throw std::invalid_argument("point not compatible with this belt");
		if (contains(pp)) // no op
			return;

		if (pp[0] == nn[0] - 1) {
			--nn[0];
			up.emplace_front(point(pp.begin()+1, pp.end()));
		}
		else if (pp[0] == nn[1] + 1) {
			++nn[1];
			up.emplace_back(point(pp.begin()+1, pp.end()));
		}
		else if (pp[0] >= nn[0] && pp[0] <= nn[1]) {
			up[pp[0] - nn[0]].add(point(pp.begin()+1, pp.end()));
		}
		else
			throw std::range_error("only adjacient points can be added");
	}


	/*    ***************    *
	 *    LIST OPERATIONS    *
	 *    ***************    */

	// return a vector the points delineating this belt
	std::vector<point> delim() const {
		std::vector<point> this_cord;
		this_cord.reserve(size());

		size_t at = nn[0];
		for (const auto &ib : up) {
			std::vector<point> up_cord = ib.delim();
			extrude_low(up_cord, this_cord, at++);
		}

		return this_cord;
	}

	// vector of all points belonging (contained) to belt
	std::vector<point> points() const {
		std::vector<point> this_cord;
		this_cord.reserve(capacity());

		size_t at = nn[0];
		for (const auto &ib : up) {
			std::vector<point> up_cord = ib.points();
			extrude_low(up_cord, this_cord, at++);
		}

		return this_cord;
	}

	/*
	// large definition of adjacient points
	std::vector<point> neighbors() const {
		std::set<point> this_cord;
		std::vector<point> all_cord = delim();
		
		size_t pos = (nn[0] > 0 ? nn[0] - 1 : nn[0]), end = nn[0]+1;
		for (const auto &ib : up) {
			std::vector<point> up_cord = ib.delim();
			std::vector<point> up_next = ib.neighbors();
			up_next.insert(up_next.end(), std::make_move_iterator(up_cord.begin()),
						      std::make_move_iterator(up_cord.end()));
			for (size_t at = pos; at <= end; ++at) {
				for (auto pp : up_next) {
					pp.push_front(at);
					// do not add belt points
					if (std::find(all_cord.begin(), all_cord.end(), pp)
							== all_cord.end())
						this_cord.insert(pp);
				}
			}
			++end;
			pos = end - 2;
		}

		return std::vector<point>(this_cord.begin(), this_cord.end());
	}
	*/

	// smallest definition of adjacient points (no diagonal)
	std::vector<point> closest() const {
		std::vector<point> this_cord;
		// by implementation nn[1] >= nn[0]
		this_cord.reserve(size());

		// first look before first point
		if (nn[0] > 0) {
			size_t at = nn[0] - 1;
			std::vector<point> up_cord = up.front().points();
			extrude_low(up_cord, this_cord, at);
		}

		// next do middle points
		size_t at = nn[0];
		for (const auto &ib : up) {
			std::vector<point> up_cord = ib.closest();
			extrude_low(up_cord, this_cord, at++);
		}

		std::vector<point> up_cord = up.back().points();
		extrude_low(up_cord, this_cord, at);

		return this_cord;
	}

	// add a dimension to an iterable list of points
	// from a upper dimension
	// extrude in this lower dimension at coordinate <coord>
	void extrude_low(std::vector<point> &up,
			std::vector<point> &lo,
			size_t coord) const {
		std::transform(up.begin(), up.end(), std::back_inserter(lo),
			[&coord](point &pp) -> point { pp.push_front(coord); return std::move(pp); });
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

	// number of pairs of coordinates
	size_t memory() const {
		return 1;
	}

	// number of points delimiting belt
	size_t size() const {
		return (nn[1] == nn[0]) ? 1 : 2;
	}

	// number of points contained in belt
	size_t capacity() const {
		return nn[1] - nn[0] + 1;
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

		if (nn[0] > 0 && pp[0] == nn[0] - 1)
			--nn[0];
		else if (pp[0] == nn[1] + 1)
			++nn[1];
		else
			throw std::range_error("only adjacient points can be added");
	}

	// vector of points delimiting the belt
	std::vector<point> delim() const {
		if (nn[1] == nn[0])
			return {point({nn[0]})};
		return {point({nn[0]}), point({nn[1]})};
	}

	// vector of all points belonging to belt
	std::vector<point> points() const {
		std::vector<point> ret;
		ret.reserve(capacity());	// reserve as many points 
		size_t at = nn[0];
		std::generate_n(std::back_inserter(ret), capacity(),
					[&at]() { return point({at++}); });
		return ret;
	}

	/*
	std::vector<point> neighbors() const {
		return closest();
	}
	*/

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
