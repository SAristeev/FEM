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
#include "mkl.h"
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
};

void read_dimensions(const json& fc, int& dim);
void read_materials(const json& fc, std::vector<material_t>& materials, std::map<int, unsigned char>& matid_threshold_map);
void read_blocks(const json& fc, std::map<int, unsigned char>& matid_threshold_map, std::vector<int>& blocks, std::vector<unsigned char>& thresholds, std::map<int, unsigned char>& block_threshold_map);
void read_mesh(const json& fc, UnstructedMesh& mesh);

void createLoads(const int& dim, const json& fc, std::vector<double>& F, const UnstructedMesh& mesh);
void applyconstraints(const json& fc, std::vector<double>& K, const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols, std::vector<double>& F, const UnstructedMesh& mesh);
void solve(const int& dim, const std::vector<double>& K, const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols, const std::vector<double>& F, std::vector<double>& x);
void buildFullGlobalMatrixStruct(const UnstructedMesh& mesh, std::vector<MKL_INT>& rows, std::vector<MKL_INT>& cols);
void buildFullGlobalMatrix(const int& dim, std::vector<double>& K, material_t mat, const UnstructedMesh& mesh, const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols);
void resultants(const int& dim, material_t material, std::vector<double>& eps, std::vector<double>& sigma, const std::vector<double>& x, const UnstructedMesh& mesh, const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols);


#endif // !__FEM_HPP__