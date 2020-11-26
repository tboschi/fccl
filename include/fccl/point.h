/* A point in a N-dimensional space of the observed events
 * per channel type (with N channels)
 */

#ifndef POINT_H
#define POINT_H

#include <deque>
#include <ostream>

using point = std::deque<size_t>;

std::ostream & operator<<(std::ostream &os, const point& pp) {
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
