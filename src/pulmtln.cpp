
#include "Driver.h"

#include <filesystem>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main()
{
	po::

	std::string inputFilename;


	std::string folder{ 
		std::filesystem::path{ inputFilename }.parent_path().string()
	};

	auto driver{ pulmtln::Driver::loadFromFile(inputFilename) };
	driver.getMTLPUL();
	driver.setExportFolder(folder);
}
