#include <cstdint>
#include <ctime>
#include "random.h"

PcgRandom::PcgRandom()
{
	state = (uint64_t)time(0);
}

PcgRandom::PcgRandom(uint64_t nseed) : state(nseed) {}

uint32_t PcgRandom::next()
{
	state = state * multiplier + increment;

	uint32_t xorshifted = ((state >> 18u) ^ state) >> 27u;
	uint32_t rot = state >> 59u;
	return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

double PcgRandom::rand()
{
	return (double)next() / 4294967296.f;
}

void PcgRandom::seed(uint64_t s)
{
	state = s;
}

void PcgRandom::seed()
{
	state = (uint64_t)time(0);
}
