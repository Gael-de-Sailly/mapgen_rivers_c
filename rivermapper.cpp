#include <cstdint>
#include <vector>
#include <unordered_map>
#include <stack>
#include <ctime>
#include <limits>
#include <cmath>
#include <iostream>
#include "rivermapper.h"
#include "random.h"

using namespace std;

RiverMapper::RiverMapper(Map<double>* dem_map) : width(dem_map->width), height(dem_map->height), size(dem_map->size)
{
	dem = dem_map->data;
	water = new double[size];
	dirs_ref = new size_t[size];
	lakes = new double[size];
	dirs = new uint8_t[size];
	ndonors = new uint8_t[size];
	basin_id = new uint32_t[size];
	pcg = new PcgRandom();
}

RiverMapper::~RiverMapper()
{
	delete[] water;
	delete[] dirs_ref;
	delete[] lakes;
	delete[] dirs;
	delete[] ndonors;
	delete[] basin_id;
	delete[] pcg;
}

inline double diff_or_zero(double a, double b) {
	return (a > b) ? a-b : 0;
}

inline double max(double a, double b) {
	return (a > b) ? a : b;
}

inline int RiverMapper::flow_local(const double zdiffs[], const int nz)
{
	double sum = 0.0f;
	for (int i=0; i<nz; i++) {
		sum += zdiffs[i];
	}
	if (sum <= 0.0f)
		return 0;

	double rn = pcg->rand() * sum;
	for (int i=0; i<nz; i++) {
		if (zdiffs[i] > rn) {
			return i+1;
		}
		rn -= zdiffs[i];
	}

	return 0;
}

inline uint64_t get_key(uint32_t k1, uint32_t k2) {
	return ((uint64_t)k1 << 32) | (uint64_t)k2;
}

struct Link {uint32_t basins[2]; size_t x; size_t y; double elev; bool isY;};

struct LinkManager {
	Link* linklist;
	unordered_map<uint64_t, Link*> linkmap;
	size_t n = 0;

	LinkManager(size_t num_basins) {
		size_t max_num_links = 1;
		if (num_basins > 2)
			max_num_links = (size_t)num_basins*3 - 6; // From Euler's characteristics, a planar graph of n vertices can have at most 3n-6 edges.
		linklist = new Link[max_num_links];
		linkmap.reserve(max_num_links);
	};

	~LinkManager() {
		delete[] linklist;
	};

	void add_link(uint32_t b1, uint32_t b2, size_t x, size_t y, double elev, bool isY) {
		if (b1 > b2) {
			uint32_t b3 = b1;
			b1 = b2; b2 = b3;
		}

		uint64_t key = get_key(b1, b2);
		auto item = linkmap.find(get_key(b1, b2));

		if (item == linkmap.end()) {
			Link* ln = &linklist[n++];
			*ln = {b1, b2, x, y, elev, isY};
			linkmap.insert({key, ln});

		} else {
			Link* ln = item->second;
			if (elev < ln->elev) {
				ln->x = x;
				ln->y = y;
				ln->elev = elev;
				ln->isY = isY;
			}
		}
	}
};

struct Adjacency{size_t next; size_t nbasin;};
struct Basin{size_t begin; size_t size;};

void planar_boruvka(Link* linklist, Link** graph, const size_t nlinks, const size_t nbasins)
{
	//cout << "Building adjacency list" << endl;
	/*
	adjacencylist: List that stores all links of linklist, in the same order, but duplicated (one for both directions). The '.next' parameters points to another link coming from the same node.
	basinlist: first edge of every node, and number of edges.
	*/
	Adjacency* adjacencylist = new Adjacency[nlinks*2];
	Basin* basinlist = new Basin[nbasins];
	bool* active_edges = new bool[nlinks];
	size_t* bucket = new size_t[nbasins];
	for (int i=0; i<nbasins; i++)
		basinlist[i] = {(size_t)-1, 0};

	for (size_t i=0; i<nlinks; i++) { // Fill the adjacency list
		active_edges[i] = true;
		Link& link = linklist[i];

		// Add link from basin 1 to basin 2
		Basin* basin = &basinlist[link.basins[0]];
		adjacencylist[i*2] = {basin->begin, link.basins[1]};
		basin->begin = i*2;
		basin->size++;

		// Add link from basin 2 to basin 1
		basin = &basinlist[link.basins[1]];
		adjacencylist[i*2+1] = {basin->begin, link.basins[0]};
		basin->begin = i*2+1;
		basin->size++;
	}


	vector<size_t> lowlevel;
	vector<size_t> highlevel;
	for (size_t i=0; i<nbasins; i++) // Select basins based on their number of links
		if (basinlist[i].size <= 8)
			lowlevel.push_back(i);
		else {
			highlevel.push_back(i);
		}

	for (size_t i=0; i<nbasins; i++)
		bucket[i] = (size_t)-1;
	vector<size_t> in_bucket;
	size_t g = 0;
	while (!lowlevel.empty()) {
		for (const size_t b1 : lowlevel) { // Process all basins with 8 edges or less
			if (basinlist[b1].size > 8) { // Check that the basin has less than 8 edges. It may have changed.
				highlevel.push_back(b1);
				continue;
			}

			size_t lowest_last_n;
			size_t lowest_n;
			double lowest_elev = numeric_limits<double>::infinity();
			size_t n = basinlist[b1].begin;
			size_t nlast = (size_t)-1;
			size_t count = 0;

			while (n != (size_t)-1) { // Loop over b1's edges
				const size_t nlink = n >> 1; // Corresponding index in linklist

				if (!active_edges[nlink]) {
					// If edge is inactive, remove it and skip
					// Removing the nth edge: making the n-1th edge pointing to the n+1th one.
					if (nlast == (size_t)-1)
						basinlist[b1].begin = adjacencylist[n].next;
					else
						adjacencylist[nlast].next = adjacencylist[n].next;
					n = adjacencylist[n].next; // Move to next link
					continue;
				}

				// If it is the lowest edge yet
				if (linklist[nlink].elev < lowest_elev) {
					lowest_elev = linklist[nlink].elev;
					lowest_last_n = nlast;
					lowest_n = n;
				}

				count++;
				nlast = n;
				n = adjacencylist[n].next; // Move to next link
			}
			if (count != basinlist[b1].size)
				cout << "Inconsistent size for node B" << b1 << ";: " << count << " links instead of " << basinlist[b1].size << endl;


			if (!isfinite(lowest_elev)) {
				// If node is empty! (Should happen only for last node)
				continue;
			}

			// Add link to the graph
			graph[g++] = &linklist[lowest_n>>1];

			const size_t b2 = adjacencylist[lowest_n].nbasin;

			// We are going to merge basin b1 into b2.
			// First, edges that were pointing to basin b1 should be redirected to b2
			n = basinlist[b1].begin;
			while (n != (size_t)-1) { // Loop again over the b1's edges
				if (adjacencylist[n].nbasin == b2) { // Links that point to b2 should be removed
					active_edges[n>>1] = false; // Mark the link as disused, it will be removed when encountered
					basinlist[b1].size--; // Decrement basin sizes (applies to both basins)
					basinlist[b2].size--;
				}

				// Redirect reciprocal edges to b2
				adjacencylist[n^(size_t)1].nbasin = b2; // n^1 corresponds to the reciprocal edge (because pairs of edges are stored next to each other)
				n = adjacencylist[n].next;
			}


			adjacencylist[nlast].next = basinlist[b2].begin; // Append b2 after b1 by referencing last b1 member on first b2 one
			basinlist[b2].begin = basinlist[b1].begin; // Merge b2 and b1
			basinlist[b2].size += basinlist[b1].size; // Increment b2's size

			basinlist[b1] = {(size_t)-1, 0}; // Delete b1
		}

		// Clean the graph, to remove many duplicate or self-loop edges that may appear
		size_t cur_highlevel = 0;
		lowlevel.clear();

		if (highlevel.size() == 0)
			break;

		for (const size_t b1 : highlevel) { // High level nodes are the only ones remaining
			in_bucket.clear();

			size_t n = basinlist[b1].begin;
			size_t j = 0;
			while (n != (size_t)-1) { // Loop over edges
				const size_t b2 = adjacencylist[n].nbasin;

				if (!active_edges[n>>1]) { // Skip disused edges
					n = adjacencylist[n].next;
					continue;
				}
				size_t jref = bucket[b2]; // Check the content of the bucket for the opposite basin. If that basin has already been seen, it contains an index to the corresponding item in 'in_bucket'.

				if (jref == (size_t)-1) { // First time a link to this basin is encountered, keep it
					bucket[b2] = j++;
					in_bucket.push_back(n);
				} else { // There is already a link to this basin
					if (linklist[n >> 1].elev < linklist[in_bucket[jref] >> 1].elev) { // Check elevation: if the new link found is lower, replace the older, otherwise, ignore it.
						active_edges[in_bucket[jref]>>1] = false;
						in_bucket[jref] = n;
					} else {
						active_edges[n>>1] = false;
					}
				}
				n = adjacencylist[n].next; // Move to next edge
			}


			size_t nlast = (size_t)-1;
			if (in_bucket.size() == 0) {
				basinlist[b1] = {(size_t)-1, 0};
				continue;
			}

			// Reconstruct adjacency of b1
			for (const size_t na : in_bucket) {// na to avoid confusion with n
				bucket[adjacencylist[na].nbasin] = (size_t)-1;
				if (nlast == (size_t)-1) {
					basinlist[b1].begin = na;
				} else
					adjacencylist[nlast].next = na;
				nlast = na;
			}
			adjacencylist[nlast].next = (size_t)-1;

			// Recompute size and push into appropriate list
			basinlist[b1].size = in_bucket.size();
			if ((basinlist[b1].size <= 8) && (basinlist[b1].size > 0))
				lowlevel.push_back(b1);
			else {
				highlevel[cur_highlevel++] = b1;
			}
		}

		highlevel.resize(cur_highlevel);
	}


	delete[] adjacencylist;
	delete[] basinlist;
	delete[] active_edges;
	delete[] bucket;
}

inline double disptime(clock_t delta) {
	double time = (double)delta / (double)CLOCKS_PER_SEC;
	return time;
}

void RiverMapper::flow()
{
	clock_t t0, t1;
	clock_t t_total = 0;
	cout << "Calculating local flow directions..." << endl;
	t0 = clock();

	// Reinitialize array
	for (size_t i=0; i<size; i++) {
		dirs[i] = 0u;
	}

	const size_t Xmax = width-1;
	const size_t Ymax = height-1;

	vector<size_t> singular;

	// Determine local flow dirs
	size_t i = 0;
	double zdiffs[4];
	for (size_t y=0; y<height; y++) {
		for (size_t x=0; x<width; x++, i++) {
			double z = dem[i];
			zdiffs[0] = (x < Xmax) ? diff_or_zero(z, dem[i+1]) : 0;
			zdiffs[1] = (y < Ymax) ? diff_or_zero(z, dem[i+width]) : 0;
			zdiffs[2] = (x > 0) ? diff_or_zero(z, dem[i-1]) : 0;
			zdiffs[3] = (y > 0) ? diff_or_zero(z, dem[i-width]) : 0;

			uint8_t d = flow_local(zdiffs, 4u);

			switch (d) {
				case 1: {
					dirs[i+1] |= (uint8_t)1u;
					dirs_ref[i] = i+1;
					break;
				}
				case 2: {
					dirs[i+width] |= (uint8_t)2u;
					dirs_ref[i] = i+width;
					break;
				}
				case 3: {
					dirs[i-1] |= (uint8_t)4u;
					dirs_ref[i] = i-1;
					break;
				}
				case 4: {
					dirs[i-width] |= (uint8_t)8u;
					dirs_ref[i] = i-width;
					break;
				}
				default: {
					singular.push_back(i);
					dirs_ref[i] = SIZE_MAX;
				}
			}
		}
	}

	t1 = clock();
	t_total += (t1-t0);
	cout << "\tCompleted in " << disptime(t1-t0) << " s" << endl;
	cout << "Percolating drainage basins (" << singular.size() << " basins)..." << endl;

	t0 = clock();

	uint32_t nb = 0;
	std::stack<size_t> s;
	for (size_t i : singular) {
		s.push(i);
		while (!s.empty()) {
			i = s.top();
			s.pop();
			uint8_t d = dirs[i];
			if (d & (uint8_t)1u)
				s.push(i-1);
			if (d & (uint8_t)2u)
				s.push(i-width);
			if (d & (uint8_t)4u)
				s.push(i+1);
			if (d & (uint8_t)8u)
				s.push(i+width);
			basin_id[i] = nb;
		}
		nb++;
	}

	const size_t basin_max = nb;

	t1 = clock();
	t_total += (t1-t0);
	cout << "\tCompleted in " << disptime(t1-t0) << " s" << endl;
	cout << "Finding lowest passes between basins..." << endl;

	t0 = clock();

	LinkManager linkmgr(basin_max+1);

	// Link basins
	i = 0;
	for (size_t y=0; y<height; y++) {
		uint32_t b1 = basin_id[i];
		linkmgr.add_link(b1, basin_max, 0, y, dem[i], false);
		i++;
		for (size_t x=1; x<width; x++, i++) {
			const uint32_t b2 = basin_id[i];
			if (b1 != b2)
				linkmgr.add_link(b1, b2, x, y, max(dem[i], dem[i-1]), false);
			b1 = b2;
		}
		linkmgr.add_link(b1, basin_max, width, y, dem[i-1], false);
	}

	for (size_t x=0; x<width; x++) {
		i = x;
		uint32_t b1 = basin_id[i];
		linkmgr.add_link(b1, basin_max, x, 0, dem[i], true);
		i+=width;
		for (size_t y=1; y<height; y++, i+=width) {
			const uint32_t b2 = basin_id[i];
			if (b1 != b2)
				linkmgr.add_link(b1, b2, x, y, max(dem[i], dem[i-width]), true);
			b1 = b2;
		}
		linkmgr.add_link(b1, basin_max, x, height, dem[i-width], true);
	}


	t1 = clock();
	t_total += (t1-t0);
	cout << "\tCompleted in " << disptime(t1-t0) << " s" << endl;
	cout << "Building basin flow tree..." << endl;

	t0 = clock();

	Link** basin_tree = new Link*[basin_max];
	planar_boruvka(linkmgr.linklist, basin_tree, linkmgr.n, basin_max+1); // This computes basin tree

	t1 = clock();
	t_total += (t1-t0);
	cout << "\tCompleted in " << disptime(t1-t0) << " s" << endl;
	cout << "Orienting basin tree and updating flow directions..." << endl;

	t0 = clock();

	// Compute adjacency list of the basin tree
	size_t* adjacencylist = new size_t[basin_max*2];
	size_t* basinlist = new size_t[basin_max+1];
	for (size_t i=0; i<basin_max+1; i++)
		basinlist[i] = (size_t)-1;

	for (size_t i=0; i<basin_max; i++) {
		Link& link = *basin_tree[i];

		adjacencylist[i*2] = basinlist[link.basins[0]];
		basinlist[link.basins[0]] = i*2;

		adjacencylist[i*2+1] = basinlist[link.basins[1]];
		basinlist[link.basins[1]] = i*2+1;
	}


	// Resolve links
	double* lake_elev = new double[basin_max+1];
	lake_elev[basin_max] = -numeric_limits<double>::infinity();
	stack<size_t> basin_stack;
	basin_stack.push(basin_max);
	while (!basin_stack.empty()) {
		const size_t b1 = basin_stack.top();
		const double b1_lake_elev = lake_elev[b1];
		basin_stack.pop();

		size_t n = basinlist[b1];
		while (n != (size_t)-1) {
			if (basin_tree[n>>1] == nullptr) {
				n = adjacencylist[n];
				continue;
			}
			Link& link = *basin_tree[n>>1];
			const size_t b2 = (link.basins[0]==b1) ? link.basins[1] : link.basins[0];

			basin_stack.push(b2);
			lake_elev[b2] = max(b1_lake_elev, link.elev);

			// Resolve link
			size_t offset = link.isY ? width : 1;
			size_t i1 = link.y*width + link.x;
			bool forward = (link.x == width) || (link.y == height) || (basin_id[i1] != b2);
			bool border = link.isY ? ((link.y==0) || (link.y==height)) : ((link.x==0) || (link.x==width));
			size_t i2;
			if (forward) {
				i2 = i1;
				i1 -= offset;
			} else {
				i2 = i1 - offset;
			}

			if (border) {
				i2 = SIZE_MAX;
			}

			while (i1 != SIZE_MAX) {
				if (basin_id[i1] != b2)
					cout << "Basin " << basin_id[i1] << " instead of " << b2 << "!" << endl;
				const size_t temp = dirs_ref[i1];
				dirs_ref[i1] = i2;
				i2 = i1;
				i1 = temp;
			}

			basin_tree[n>>1] = nullptr; // To avoid looping twice on the same link
			n = adjacencylist[n];
		}
	}

	for (size_t i=0; i<size; i++) {
		lakes[i] = lake_elev[basin_id[i]];

		const int dir_diff = dirs_ref[i] - i;
		uint8_t d = 0u;
		if (dir_diff==-width)
			d = 1u;
		else if (dir_diff==1)
			d = 2u;
		else if (dir_diff==width)
			d = 3u;
		else if (dir_diff==-1)
			d = 4u;
		dirs[i] = d;
	}
				


	delete[] basin_tree;
	delete[] adjacencylist;
	delete[] basinlist;
	delete[] lake_elev;

	t1 = clock();
	t_total += (t1-t0);
	cout << "\tCompleted in " << disptime(t1-t0) << " s" << endl;
	cout << "Accumulating water..." << endl;

	t0 = clock();

	for (size_t i=0; i<size; i++)
		water[i] = 1;
	accumulate();

	t1 = clock();
	t_total += (t1-t0);
	cout << "\tCompleted in " << disptime(t1-t0) << " s" << endl;
	cout << "Full flow calculation completed in " << disptime(t_total) << " s" << endl;
}


void RiverMapper::accumulate()
{
	for (size_t i=0; i<size; i++)
		ndonors[i] = 0;
	for (size_t i=0; i<size; i++) {
		size_t d = dirs_ref[i];
		if (d == SIZE_MAX)
			continue;
		ndonors[d]++;
	}

	for (size_t i=0; i<size; i++) {
		if (ndonors[i] > 0)
			continue;

		double* w1 = &water[i];
		size_t iw = i;
		while(true) {
			iw = dirs_ref[iw];

			if (iw == SIZE_MAX) {
				break;
			}

			double* w2 = &water[iw];
			*w2 += *w1; // Flow w1 into w2
			w1 = w2; // Now we will follow w2
			if (ndonors[iw] > 1) {
				ndonors[iw]--;
				break;
			}
		}
	}
}
