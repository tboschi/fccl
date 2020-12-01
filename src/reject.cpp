/* given expected background, this code finds
 * the minimum expected signal to reject the
 * null hypothesis
 * 	H_0 = "there is no signal, just background"
 */

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <getopt.h>

#include "fccl/region.h"

int main(int argc, char** argv) {

	double CL = 0.90;
	int iarg = 0; 
	while((iarg = getopt(argc,argv, "C:h")) != -1)
	{
		switch(iarg)
		{
			case 'C':
				CL = strtod(optarg, NULL);
				CL /= 100.0;
				break;
			case 'h':
				return 1;
			default:
				break;
		}
	}
	if (argv[optind] == NULL) {
		std::cerr << "Missing mandatory argument\n";
		return 1;
	}

	region<1> reg;

	double bak = std::atof(argv[optind]);
	double sig = sqrt(bak);

	// background
	point bb = {size_t(bak)};
	belt<1> myb(bb);
	while (myb.contains(bb)) {
		sig += 0.001;
		rate ex(sig, bak);	// expected rate

		myb = reg.expand(CL, ex);

		// when background is not in belt anymore
		// this loop will end
	}

	auto plist = myb.points();

	std::cout << "For background of " << bak << " mean signal is " << sig << std::endl;
	std::cout << "Confidence belt (" << CL*100. << "\%) "
		  << "is contained between " << plist.front() << " and " << plist.back() << std::endl;

	return 0;
}
