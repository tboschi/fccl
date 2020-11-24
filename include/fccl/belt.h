#ifndef BELT_H
#define BELT_H

#include <cmath>
#include <deque>
#include <memory>
#include <numeric>

// N D belt
template <int N, typename std::enable_if<(N > 0), bool>::type = true>
struct belt {
	belt() = default;

	template<typename T, typename... Args>
	belt(T t, Args... args) {
		nn = {t, t};
		up.push_back(belt<N-1>(args...));
	}

	std::array<size_t, 2> nn;
	std::deque<belt<N-1> > up;

	size_t size() const {
		return 1 + std::accumulate(up.begin(), up.end(), 0,
				[](size_t sz, const belt<N-1> &b) {
					return sz += b.size(); } );
	}
};

// 1 D belt
template<>
struct belt<1> {
	belt() = default;

	template<typename T>
	belt(T t) {
		nn = {t, t};
	}

	std::array<size_t, 2> nn;

	size_t size() const {
		return 1;
	}
};

#endif
