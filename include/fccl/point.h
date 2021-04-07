/* A point in a N-dimensional space of the observed events
 * per channel type (with N channels)
 */

#ifndef POINT_H
#define POINT_H

#include <deque>
#include <ostream>

namespace fccl {
	using point = std::deque<size_t>;
} /*fccl*/

inline std::ostream & operator<<(std::ostream &os, const fccl::point& pp) {
	const auto separator = ", ";
	const auto* del = "";
	os << "<";
	for (const auto &ip : pp) {
		os << del << ip;
		del = separator;
	}
	return os << ">";
}

#endif
