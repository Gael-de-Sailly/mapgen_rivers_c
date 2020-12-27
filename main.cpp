#include <fstream>
#include <cstdint>
#include "map.cpp"
#include "rivermapper.h"

inline bool file_exists(const char* fname) {
	std::ifstream file(fname);
	return file.good();
}

int main(int argc, char** argv)
{
	const char* file_to_read = (argc > 1) ? argv[1] : "dem";
	if (file_exists(file_to_read)) {
		std::cout << "Reading file '" << file_to_read << "' ..." << std::endl;
		Map<double>* dem = read_file<double>(file_to_read);
		RiverMapper rm(dem);
		rm.flow();
		rm.erode(1.0f);
		rm.flow();

		std::cout << "Writing files..." << std::endl;
		const char* file_dem = (argc > 2) ? argv[2] : "dem_new";
		dem->write_file(file_dem);

		Map<uint8_t> dirs(dem->width, dem->height, rm.dirs);
		const char* file_dirs = (argc > 3) ? argv[3] : "dirs";
		dirs.write_file(file_dirs);

		Map<double> water(dem->width, dem->height, rm.water);
		const char* file_rivers = (argc > 4) ? argv[4] : "rivers";
		water.write_file(file_rivers);

		Map<double> lakes(dem->width, dem->height, rm.lakes);
		const char* file_lakes = (argc > 5) ? argv[5] : "lakes";
		lakes.write_file(file_lakes);
		std::cout << "Done." << std::endl;
	} else {
		std::cout << "No file named '" << file_to_read << "'. Run 'genmap.py' to generate it." << std::endl;
	}

	return 0;
}
