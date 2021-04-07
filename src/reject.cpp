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
				if (CL >= 1.)
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

	double bak = std::atof(argv[optind]);

	// null hypothesis point, i.e. only observed background
	fccl::point bb{size_t(bak)};

	// expected rate of signal and background
	// starting signal is at least sqrt(bkg)
	fccl::rate ex(std::sqrt(bkg), bkg);

	// object to calculate F&C CL region
	fccl::region<1> reg;


	// expand region and find CL belt
	fccl::belt<1> myb = reg.expand(CL, ex);
	// continue expanding 1D CL region until background leaves region
	while (myb.contains(bb)) {
		// increment signal
		ex.sig += 0.001;

		myb = reg.expand(CL, ex);
	}

	auto plist = myb.points();

	std::cout << "For background of " << bak << " mean signal is " << sig << std::endl;
	std::cout << "Confidence belt (" << CL*100. << "\%) "
		  << "is contained between " << plist.front() << " and " << plist.back() << std::endl;

	return 0;
}
