/* A rate of event is defined by a signal and a background
 */

#ifndef RATE_H
#define RATE_H

#include <ostream>

namespace fccl {
	struct rate {
		double sig;
		double bak;

		template<typename S, typename B>
			rate(S s, B b = 0) : sig(s), bak(b) { };

		//rate(const rate &r) : sig(r.sig), bak(r.bak) { };
	};
} /*fccl*/

inline std::ostream & operator<<(std::ostream &os, const fccl::rate& rr) {
	return os << "{" << rr.sig << ", " << rr.bak << "}";
}

#endif
