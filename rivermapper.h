#pragma once

#include <cstdint>
#include "map.cpp"
#include "random.h"

class RiverMapper
{
public:
	RiverMapper(Map<double>* dem_map);
	~RiverMapper();

	double* dem;
	double* water;
	uint8_t* dirs;
	double* lakes;

	void flow();
	void accumulate();
	void erode(const double time, const double K=1.0f, const double m=0.5f, const double sea_level=0.0f);

	const std::size_t width;
	const std::size_t height;
	const std::size_t size;

private:
	int flow_local(const double zdiffs[], const int n);
	size_t* dirs_ref;
	uint8_t* ndonors;
	uint32_t* basin_id;
	double* erosion_time;
	PcgRandom* pcg;
};
