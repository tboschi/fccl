/* A rate of event is defined by a signal and a background
 */

#ifndef RATE_H
#define RATE_H

#include <ostream>

struct rate {
	double sig;
	double bak;

	template<typename S, typename B>
	rate(S s, B b = 0) : sig(s), bak(b) { };
};

std::ostream & operator<<(std::ostream &os, const rate& rr) {
	return os << "{" << rr.sig << ", " << rr.bak << "}";
}

#endif
