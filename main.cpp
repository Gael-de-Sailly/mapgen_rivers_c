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
		Map<double>* map = read_file<double>(file_to_read);
		RiverMapper rm(map);
		rm.flow();

		std::cout << "Writing files..." << std::endl;
		Map<uint8_t> dirs(map->width, map->height, rm.dirs);
		const char* file_dirs = (argc > 2) ? argv[2] : "dirs";
		dirs.write_file(file_dirs);

		Map<double> water(map->width, map->height, rm.water);
		const char* file_rivers = (argc > 3) ? argv[3] : "rivers";
		water.write_file(file_rivers);

		Map<double> lakes(map->width, map->height, rm.lakes);
		const char* file_lakes = (argc > 4) ? argv[4] : "lakes";
		lakes.write_file(file_lakes);
		std::cout << "Done." << std::endl;
	} else {
		std::cout << "No file named '" << file_to_read << "'. Run 'genmap.py' to generate it." << std::endl;
	}

	return 0;
}
