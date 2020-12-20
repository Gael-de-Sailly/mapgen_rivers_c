#pragma once

#include <cstdint>

class PcgRandom
{
public:
	PcgRandom(uint64_t nseed);
	PcgRandom();
	uint32_t next();
	double rand();
	void seed(uint64_t s);
	void seed();

private:
	static const uint64_t multiplier = 6364136223846793005u;
	static const uint64_t increment = 1442695040888963407u;
	uint64_t state;
};
