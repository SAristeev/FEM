#pragma once
#ifndef __FEM_HPP__
#define __FEM_HPP__
#include "external/json/json.hpp"
#include "external/base64/base64.h"
#include "fp.hpp"
#include <set>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mkl.h>
#include <tbb/tick_count.h>
#include <tbb/parallel_for.h>
#include <cmath>
#include <span>
using json = nlohmann::json;

struct material_t
{
	unsigned char id;
	int threshold;

	fp::fp_t E;
	fp::fp_t nu;

	static material_t make_materal(unsigned char id, int threshold, fp::fp_t E, fp::fp_t nu)
	{
		return material_t{ id, threshold, E, nu };
	}
};

struct UnstructedMesh {
	std::vector<vec3<double>> nodes;
	std::vector<int>          nids;
	std::vector<int>          elems;
	std::vector<int>          elemids;
	std::vector<int>          nodes_per_elem;
	std::vector<uint8_t>      elem_type;
	std::map<int, int>        map_node_numeration;
	std::map<int, int>        map_element_numeration;
};

void read_dimensions(const json& fc, int& dim);
void read_materials(const json& fc, std::vector<material_t>& materials, std::map<int, unsigned char>& matid_threshold_map);
void read_blocks(const json& fc, std::map<int, unsigned char>& matid_threshold_map, std::vector<int>& blocks, std::vector<unsigned char>& thresholds, std::map<int, unsigned char>& block_threshold_map);
void read_mesh(const json& fc, UnstructedMesh& mesh);

void createLoads(const int& dim, const json& fc, std::span<double>& F, const UnstructedMesh& mesh);
void applyconstraints(const json& fc, std::span<double>& K, const std::span<MKL_INT>& rows, const std::span<MKL_INT>& cols, std::span<double>& F, const UnstructedMesh& mesh);
void solve(const MKL_INT& blocksize, const std::span<double>& A, const std::span<MKL_INT>& rows, const std::span<MKL_INT>& cols, const MKL_INT& nrhs, const std::span<double>& b, std::span<double>& x);
void buildFullGlobalMatrixStruct(const UnstructedMesh& mesh, std::span<MKL_INT>& rows, std::span<MKL_INT>& cols);
void buildFullGlobalMatrix(const int& dim, std::span<double>& K, material_t mat, const UnstructedMesh& mesh, const std::span<MKL_INT>& rows, const std::span<MKL_INT>& cols);
void resultants(const int& dim, material_t material, std::span<double>& eps, std::span<double>& sigma, const std::span<double>& x, const UnstructedMesh& mesh, const std::span<MKL_INT>& rows, const std::span<MKL_INT>& cols);

void print_matrix(const int& dim, std::span<double>& K, const std::span<MKL_INT>& rows, const std::span<MKL_INT>& cols);
#endif // !__FEM_HPP__