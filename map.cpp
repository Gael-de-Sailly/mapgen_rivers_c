#pragma once

#include <iostream>
#include <fstream>
#include <cstdint>

using namespace std;

template<typename T>
class Map
{
public:
	Map(size_t widthp, size_t heightp);
	Map(size_t widthp, size_t heightp, T* content);
	Map(size_t widthp, size_t heightp, ifstream& in_file);
	~Map();
	T* data;
	const size_t width;
	const size_t height;
	const size_t size;

	void write_file(const char* filename);
	void print();
private:
	const bool original_data;
};

template<typename T>
Map<T>::Map(size_t widthp, size_t heightp) : width(widthp), height(heightp), size(widthp*heightp), original_data(true)
{
	data = new T [size];
}

template <typename T>
Map<T>::Map(size_t widthp, size_t heightp, T* content) : width(widthp), height(heightp), size(widthp*heightp), original_data(false), data(content) {}

template<typename T>
Map<T>::Map(size_t widthp, size_t heightp, ifstream& in_file) : width(widthp), height(heightp), size(widthp*heightp), original_data(true)
{
	data = new T [size];
	in_file.read(reinterpret_cast<char*>(data), size*sizeof(T));
}

template<typename T>
Map<T>::~Map()
{
	if (original_data)
		delete[] data;
}

template<typename T>
Map<T>* read_file(const char* filename)
{
	ifstream in_file(filename, ios::in|ios::binary);
	in_file.seekg(0, ios::beg);

	uint16_t width, height;
	in_file.read(reinterpret_cast<char*>(&width), sizeof(uint16_t));
	in_file.read(reinterpret_cast<char*>(&height), sizeof(uint16_t));

	Map<T>* output_map = new Map<T>((size_t)width, (size_t)height, in_file);

	in_file.close();
	return output_map;
}

template<typename T>
void Map<T>::write_file(const char* filename)
{
	ofstream out_file(filename, ios::out|ios::binary);

	size_t widthp = width;
	size_t heightp = height;
	out_file.write(reinterpret_cast<char*>(&widthp), sizeof(uint16_t));
	out_file.write(reinterpret_cast<char*>(&heightp), sizeof(uint16_t));

	size_t write_size = size*sizeof(T);
	out_file.write(reinterpret_cast<char*>(data), size*sizeof(T));

	out_file.close();
}

template<typename T>
void Map<T>::print()
{
	cout << width << endl;
	cout << height << endl;
	size_t i=0;
	for (size_t y=0; y<height; y++) {
		for (size_t x=0; x<width; x++, i++) {
			cout << data[i] << "\t";
		}
		cout << endl;
	}
}
