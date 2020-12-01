#include <iostream>
#include <fstream>
#include <cstdlib>

#include "fccl/region.h"

int main(int argc, char** argv) {

	region<2> dirac, major;

	double bak = std::atof(argv[1]);
	double sig = sqrt(bak);

	double CL = 0.99;	//confidence level
	double scale = 0.2;	// LNV wrt to LNC channel

	std::ofstream outd("lnv_dirac.dat");
	std::ofstream outm("lnv_major.dat");
	std::ofstream adjd("adj_dirac.dat");
	std::ofstream adjm("adj_major.dat");
	int index = 0;
	for (double sig = 1; sig < 100; ++sig)
	{
		rate LNC_dirac(sig, bak);
		rate LNV_dirac(sig * scale, bak);

		rate LNC_major(sig * (1 + scale) / 2., bak);
		rate LNV_major(sig * (1 + scale) / 2., bak);

		auto dirac_belt = dirac.expand(CL, LNC_dirac, LNV_dirac);
		auto major_belt = major.expand(CL, LNC_major, LNV_major);

		std::cout << index++ << " Signal " << sig << " is good for LNV? "
			  << std::boolalpha << !dirac_belt.share(major_belt) << std::endl;

		auto pdirac = dirac_belt.points();
		for (const auto p : pdirac)
			outd << p[0] << "\t" << p[1] << "\n";
		outd << "\n\n";
		auto adirac = dirac_belt.closest();
		for (const auto p : adirac)
			adjd << p[0] << "\t" << p[1] << "\n";
		adjd << "\n\n";

		auto pmajor = major_belt.points();
		for (const auto p : pmajor)
			outm << p[0] << "\t" << p[1] << "\n";
		outm << "\n\n";
		auto amajor = major_belt.closest();
		for (const auto p : amajor)
			adjm << p[0] << "\t" << p[1] << "\n";
		adjm << "\n\n";

	}

	outd.close();
	outm.close();
	adjd.close();
	adjm.close();

	return 0;
}
