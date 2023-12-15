#include "FEM.hpp"

void read_dimensions(const json& fc, int& dim) {
	if (fc["settings"]["dimensions"] == "2D") {
		dim = 2;
	}
	if (fc["settings"]["dimensions"] == "3D") {
		dim = 3;
	}
}

void read_materials(const json& fc, std::vector<material_t>& materials, std::map<int, unsigned char>& matid_threshold_map) {
	unsigned int mat_id = 1;
	for (auto& mat : fc["materials"]) {
		std::string e_b64 = mat["elasticity"][0]["constants"][0];
		std::string nu_b64 = mat["elasticity"][0]["constants"][1];

		double E, nu;
		b64decode_host(e_b64.data(), e_b64.size(), reinterpret_cast<char*>(&E));
		b64decode_host(nu_b64.data(), nu_b64.size(), reinterpret_cast<char*>(&nu));

		materials.emplace_back(material_t::make_materal(mat_id, mat_id, E, nu));

		matid_threshold_map[mat["id"]] = mat_id;
		mat_id++;

	}
}

void read_blocks(const json& fc, std::map<int, unsigned char>& matid_threshold_map, std::vector<int>& blocks, std::vector<unsigned char>& thresholds, std::map<int, unsigned char>& block_threshold_map) {
	// create block-threshold map
	for (auto& block : fc["blocks"]) {
		block_threshold_map[block["id"]] = matid_threshold_map[block["material_id"]];
	}

	for (auto const& imap : block_threshold_map) {
		blocks.push_back(imap.first);
		thresholds.push_back(imap.second);
	}
}

void read_mesh(const json& fc, UnstructedMesh& mesh) {
	// read base64 data
	auto get_mesh_data = [&fc](std::string const& key, size_t* out_size = nullptr) -> char*
		{
			std::string data64 = fc["mesh"][key];
			char* data = (char*)malloc(data64.size() * sizeof(char));
			b64decode_host(data64.data(), data64.size(), data);
			size_t size = b64dsize_host(data64.size(), data64.data());

			if (out_size) *out_size = size;

			return data;
		};
	size_t elems_count = fc["mesh"]["elems_count"];
	size_t nodes_count = fc["mesh"]["nodes_count"];


	std::vector<int> nodes_per_elem(elems_count + 1);
	std::vector<uint8_t> elem_type(elems_count + 1);

	// Only 2d Tri mesh
	int offset = 0;
	unsigned char* etypes = reinterpret_cast<unsigned char*>(get_mesh_data("elem_types"));
	for (int elem_ID = 0; elem_ID < elems_count; elem_ID++) {
		switch (etypes[elem_ID]) {
		case '\n':
			nodes_per_elem[elem_ID] = offset;
			elem_type[elem_ID] = 5;
			offset += 3;
			break;
		default:
			throw std::runtime_error("non tri mesh.\n");
			break;
		}
	}
	nodes_per_elem[elems_count] = offset;
	free(etypes);




	int* nids_raw = reinterpret_cast<int*>(get_mesh_data("nids"));
	std::vector<int> nids(nids_raw, nids_raw + nodes_count);

	for (int node_id = 0; node_id < nodes_count; node_id++) {
		if (!mesh.map_node_numeration.insert_or_assign(nids_raw[node_id], node_id).second) {
			throw std::runtime_error("Some nodes with the same ID in the mesh.\nToDo: Verify mesh nodes");
		}
	}
	free(nids_raw);

	std::map<int, int> map_element_numeration;
	int* elemids_raw = reinterpret_cast<int*>(get_mesh_data("elemids"));
	std::vector<int>  elemids(elemids_raw, elemids_raw + elems_count);
	for (int elem_id = 0; elem_id < elems_count; elem_id++) {
		if (!map_element_numeration.insert_or_assign(elemids_raw[elem_id], elem_id).second) {
			throw std::runtime_error("Some elements with the same ID in the mesh.\nToDo: Verify mesh elements");
		}
	}
	free(elemids_raw);


	vec3<double>* nodes_raw = reinterpret_cast<vec3<double>*>(get_mesh_data("nodes"));
	std::vector<vec3<double>> nodes(nodes_raw, nodes_raw + nodes_count);
	free(nodes_raw);

	size_t elems_size;
	int* elems_raw = reinterpret_cast<int*>(get_mesh_data("elems", &elems_size));
	elems_size = elems_size / sizeof(int);
	std::vector<int> elems(elems_raw, elems_raw + elems_size);

	free(elems_raw);

	mesh.nodes = nodes;
	mesh.elems = elems;
	mesh.nids = nids;
	mesh.elemids = elemids;
	mesh.elem_type = elem_type;
	mesh.nodes_per_elem = nodes_per_elem;
}


void buildFullGlobalMatrixStruct(const UnstructedMesh& mesh, std::vector<MKL_INT>& rows, std::vector<MKL_INT>& cols) {
	if (!rows.empty() || !cols.empty()) {
		throw std::runtime_error("buildFullGlobalMatrixStruct: try to fill non-empty matrix");
	}
	size_t elems_size = mesh.elemids.size();
	size_t n = mesh.nids.size();

	std::set<std::pair<int, int>> coostruct;
	int offset = 0;
	for (int elem_id = 0; elem_id < elems_size; elem_id++) {
		for (int nodei = 0; nodei < mesh.nodes_per_elem[elem_id + 1] - mesh.nodes_per_elem[elem_id]; nodei++) {
			for (int nodej = 0; nodej < mesh.nodes_per_elem[elem_id + 1] - mesh.nodes_per_elem[elem_id]; nodej++)
				coostruct.insert({ mesh.elems[offset + nodei], mesh.elems[offset + nodej] });
		}
		offset += mesh.nodes_per_elem[elem_id + 1] - mesh.nodes_per_elem[elem_id];
	}
	cols.resize(coostruct.size());
	rows.resize(n + 1);

	int currow = 0;
	int curnnz = rows[0] = 0;

	for (auto x : coostruct) {
		// zero based indexing
		if (x.first == currow + 2) {
			rows[++currow] = curnnz;
		}
		cols[curnnz] = x.second - 1;
		curnnz++;
	}
	rows[++currow] = curnnz;
}


void buildFullGlobalMatrix(const int& dim, std::vector<double>& K, material_t material, const UnstructedMesh& mesh, const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols) {
	if (dim != 2) {
		throw std::runtime_error("buildFullGlobalMatrix: try to solve non 2D task");
	}
	int blocksize = dim;
	size_t n = rows.size() - 1;
	size_t nnz = rows[n];

	double* raw_K = reinterpret_cast<double*>(malloc(blocksize * blocksize * nnz * sizeof(double)));
	for (int i = 0; i < nnz * blocksize * blocksize; i++) {
		raw_K[i] = 0.0;
	}
	
	// D matrix	
	int Ddim = 3;

	double* D = reinterpret_cast<double*>(malloc(Ddim * Ddim * sizeof(double)));

	// filling D
	fp::fp_t E = material.E;
	fp::fp_t nu = material.nu;

	D[0] = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
	D[1] = E * nu / ((1 + nu) * (1 - 2 * nu));
	D[2] = 0;
	D[3] = E * nu / ((1 + nu) * (1 - 2 * nu));
	D[4] = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
	D[5] = 0;
	D[6] = 0;
	D[7] = 0;
	D[8] = E / (2 * (1 + nu));

	// only if Tri
	int Brows = 3;
	int Bcols = 6;
	double* B = reinterpret_cast<double*>(malloc(Brows * Bcols * sizeof(double)));
	double* A = reinterpret_cast<double*>(malloc(Bcols * Bcols * sizeof(double)));
	double* Z = reinterpret_cast<double*>(malloc(Ddim * Bcols * sizeof(double)));

	for (int e = 0; e < mesh.elemids.size(); e++) {
		std::vector<int> ijk = { mesh.elems[mesh.nodes_per_elem[e] + 0] - 1,
									mesh.elems[mesh.nodes_per_elem[e] + 1] - 1,
									mesh.elems[mesh.nodes_per_elem[e] + 2] - 1 };
		//std::sort(ijk.begin(), ijk.end());

		int i = ijk[0];
		int j = ijk[1];
		int k = ijk[2];

		double S = std::abs((mesh.nodes[j].x - mesh.nodes[i].x) * (mesh.nodes[k].y - mesh.nodes[i].y) -
			(mesh.nodes[k].x - mesh.nodes[i].x) * (mesh.nodes[j].y - mesh.nodes[i].y)) / 2;
		
		B[0] = mesh.nodes[j].y - mesh.nodes[k].y;
		B[1] = 0;
		B[2] = mesh.nodes[k].x - mesh.nodes[j].x;

		B[3] = 0;
		B[4] = mesh.nodes[k].x - mesh.nodes[j].x;
		B[5] = mesh.nodes[j].y - mesh.nodes[k].y;

		B[6] = mesh.nodes[k].y - mesh.nodes[i].y;
		B[7] = 0;
		B[8] = mesh.nodes[i].x - mesh.nodes[k].x;

		B[9] = 0;
		B[10] = mesh.nodes[i].x - mesh.nodes[k].x;
		B[11] = mesh.nodes[k].y - mesh.nodes[i].y;

		B[12] = mesh.nodes[i].y - mesh.nodes[j].y;
		B[13] = 0;
		B[14] = mesh.nodes[j].x - mesh.nodes[i].x;

		B[15] = 0;
		B[16] = mesh.nodes[j].x - mesh.nodes[i].x;
		B[17] = mesh.nodes[i].y - mesh.nodes[j].y;
	
		CBLAS_LAYOUT layout = CblasColMajor;
		CBLAS_TRANSPOSE nontrans = CblasNoTrans;
		CBLAS_TRANSPOSE trans = CblasTrans;
		const double beta = 0.0;
		cblas_dgemm(layout, nontrans, nontrans, Ddim, Bcols, Ddim, 1.0 / (2 * S), D, Ddim, B, Brows, beta, Z, Ddim);
		cblas_dgemm(layout, trans, nontrans, Bcols, Bcols, Ddim, 0.5, B, Ddim, Z, Brows, beta, A, Bcols);
		
		int i1 = -1;
		int j1 = -1;
		int k1 = -1;
		for (int q = rows[i]; q < rows[i + 1]; q++) {
			//zero or one indexing
			if (cols[q] == i) {
				i1 = q;
			}
			if (cols[q] == j) {
				j1 = q;
			}
			if (cols[q] == k) {
				k1 = q;
			}
		}
		raw_K[4 * i1 + 0] += A[0 * Bcols + 0];
		raw_K[4 * i1 + 1] += A[1 * Bcols + 0];
		raw_K[4 * i1 + 2] += A[0 * Bcols + 1];
		raw_K[4 * i1 + 3] += A[1 * Bcols + 1];

		raw_K[4 * j1 + 0] += A[2 * Bcols + 0];
		raw_K[4 * j1 + 1] += A[3 * Bcols + 0];
		raw_K[4 * j1 + 2] += A[2 * Bcols + 1];
		raw_K[4 * j1 + 3] += A[3 * Bcols + 1];

		raw_K[4 * k1 + 0] += A[4 * Bcols + 0];
		raw_K[4 * k1 + 1] += A[5 * Bcols + 0];
		raw_K[4 * k1 + 2] += A[4 * Bcols + 1];
		raw_K[4 * k1 + 3] += A[5 * Bcols + 1];

		int i2 = -1;
		int j2 = -1;
		int k2 = -1;
		for (int q = rows[j]; q < rows[j + 1]; q++) {
			//zero or one indexing
			if (cols[q] == i) {
				i2 = q;
			}
			if (cols[q] == j) {
				j2 = q;
			}
			if (cols[q] == k) {
				k2 = q;
			}
		}
		raw_K[4 * i2 + 0] += A[0 * Bcols + 2];
		raw_K[4 * i2 + 1] += A[1 * Bcols + 2];
		raw_K[4 * i2 + 2] += A[0 * Bcols + 3];
		raw_K[4 * i2 + 3] += A[1 * Bcols + 3];
		
		raw_K[4 * j2 + 0] += A[2 * Bcols + 2];
		raw_K[4 * j2 + 1] += A[3 * Bcols + 2];
		raw_K[4 * j2 + 2] += A[2 * Bcols + 3];
		raw_K[4 * j2 + 3] += A[3 * Bcols + 3];
		
		raw_K[4 * k2 + 0] += A[4 * Bcols + 2];
		raw_K[4 * k2 + 1] += A[5 * Bcols + 2];
		raw_K[4 * k2 + 2] += A[4 * Bcols + 3];
		raw_K[4 * k2 + 3] += A[5 * Bcols + 3];

		int i3 = -1;
		int j3 = -1;
		int k3 = -1;
		for (int q = rows[k]; q < rows[k + 1]; q++) {
			//zero or one indexing
			if (cols[q] == i) {
				i3 = q;
			}
			if (cols[q] == j) {
				j3 = q;
			}
			if (cols[q] == k) {
				k3 = q;
			}
		}
		raw_K[4 * i3 + 0] += A[0 * Bcols + 4];
		raw_K[4 * i3 + 1] += A[1 * Bcols + 4];
		raw_K[4 * i3 + 2] += A[0 * Bcols + 5];
		raw_K[4 * i3 + 3] += A[1 * Bcols + 5];
						 
		raw_K[4 * j3 + 0] += A[2 * Bcols + 4];
		raw_K[4 * j3 + 1] += A[3 * Bcols + 4];
		raw_K[4 * j3 + 2] += A[2 * Bcols + 5];
		raw_K[4 * j3 + 3] += A[3 * Bcols + 5];
						 
		raw_K[4 * k3 + 0] += A[4 * Bcols + 4];
		raw_K[4 * k3 + 1] += A[5 * Bcols + 4];
		raw_K[4 * k3 + 2] += A[4 * Bcols + 5];
		raw_K[4 * k3 + 3] += A[5 * Bcols + 5];

	}
	K = std::vector<double>(raw_K, raw_K + blocksize * blocksize * nnz);
	free(raw_K);
	free(B);
	free(A);
	free(Z);
	free(D);
}

void createLoads(const int& dim, const json& fc, std::vector<double>& F, const UnstructedMesh& mesh) {
	F.resize(dim * mesh.nodes.size());
	for (auto& load : fc["loads"]) {
		auto get_load_data = [&load](std::string const& key, size_t* out_size = nullptr) -> char*
			{
				std::string data64 = load[key];
				char* data = (char*)malloc(data64.size() * sizeof(char));
				b64decode_host(data64.data(), data64.size(), data);
				size_t size = b64dsize_host(data64.size(), data64.data());

				if (out_size) *out_size = size;

				return data;
			};
		if (load["name"] == "Force" && load["type"] == 5) {
			
			size_t apply_to_size = load["apply_to_size"];

			int* apply_to_raw = reinterpret_cast<int*>(get_load_data("apply_to"));
			std::vector<int> apply_to(apply_to_raw, apply_to_raw + apply_to_size);
			free(apply_to_raw);
			std::vector<double> data(6);
			for (int i = 0; i < 6; i++) {
				std::string data_b64 = load["data"][i];
				b64decode_host(data_b64.data(), data_b64.size(), reinterpret_cast<char*>(&(data[i])));
			}

			for (int& node : apply_to) {
				for (int i = 0; i < 2; i++) {
					F[2 * mesh.map_node_numeration.at(node) + i] += data[i];
				}
			}
		}		

		else if (load["name"] == "Pressure" && load["type"] == 4) {
			size_t apply_to_size = load["apply_to_size"];

			int* apply_to_raw = reinterpret_cast<int*>(get_load_data("apply_to"));
			std::vector<int> apply_to(apply_to_raw, apply_to_raw + 2 * apply_to_size);
			free(apply_to_raw);

			std::string data_b64 = load["data"][0];
			double data;
			b64decode_host(data_b64.data(), data_b64.size(), reinterpret_cast<char*>(&data));
			for (int p = 0; p < apply_to_size; p++) {
				int elem = 2 * p;
				int edge = 2 * p + 1;
				int shift1 = apply_to[edge];
				int shift2 = (apply_to[edge] + 1) % 3;
				int shift3 = (apply_to[edge] + 2) % 3;

				int i = mesh.elems[mesh.nodes_per_elem[apply_to[elem] - 1] + shift1] - 1;
				int j = mesh.elems[mesh.nodes_per_elem[apply_to[elem] - 1] + shift2] - 1;
				int k = mesh.elems[mesh.nodes_per_elem[apply_to[elem] - 1] + shift3] - 1;

				double nx = -(mesh.nodes[j].y - mesh.nodes[i].y);
				double ny = mesh.nodes[j].x - mesh.nodes[i].x;

				if (nx * (mesh.nodes[k].x - mesh.nodes[i].x) + ny * (mesh.nodes[k].y - mesh.nodes[i].y) < 0) {
					double tmp = nx;
					nx = -ny;
					ny = tmp;
				}
				double len = std::sqrt(nx * nx + ny * ny);
				nx /= len;
				ny /= len;

				F[2 * i + 0] += data * nx * len / 2;
				F[2 * i + 1] += data * ny * len / 2;

				F[2 * j + 0] += data * nx * len / 2;
				F[2 * j + 1] += data * ny * len / 2;
			}
		}
		else {
			throw std::runtime_error("not Force or Pressure load. Not supported yet");
		}
		

		// Force to pressure
		/*if (std::abs(data[0]) > 1e-8 ) {
			std::map<double, int> x_side;
			for (int p = 0; p < apply_to_size; p++) {
				x_side[mesh.nodes[mesh.map_node_numeration.at(apply_to[p])].y] = p;
			}

			auto find_by_value = [&x_side](int target) {
				for (const auto& [key, value] : x_side)
					if (value == target)
						return key;
				};
			
			auto begin0 = find_by_value(0);
			auto begin1 = find_by_value(1);
			F[2 * mesh.map_node_numeration.at(apply_to[0])] += data[0] * std::abs((begin0 - begin1)) / 2 ;

			for(int p = 1; p < apply_to_size - 1; p++) {
				auto prev = find_by_value(p - 1);
				auto cur  = find_by_value(p);
				auto next = find_by_value(p + 1);
				F[2 * mesh.map_node_numeration.at(apply_to[p])] += data[0] * std::abs((next - cur)) / 2;
				F[2 * mesh.map_node_numeration.at(apply_to[p])] += data[0] * std::abs((prev - cur)) / 2;
			}
			auto end0 = find_by_value(apply_to_size - 2);
			auto end1 = find_by_value(apply_to_size - 1);
			F[2 * mesh.map_node_numeration.at(apply_to[apply_to_size - 1])] += data[0] * std::abs((end0 - end1)) / 2;

		}
		if (std::abs(data[1]) > 1e-8) {
			std::map<double, int> y_side;
			for (int p = 0; p < apply_to_size; p++) {
				y_side[mesh.nodes[mesh.map_node_numeration.at(apply_to[p])].x] = p;
			}
			auto find_by_value = [&y_side](int target) {
				for (const auto& [key, value] : y_side)
					if (value == target)	
						return key;
				};


			auto begin0 = find_by_value(0);
			auto begin1 = find_by_value(1);
			F[2 * mesh.map_node_numeration.at(apply_to[0]) + 1] += data[1] * std::abs((begin0 - begin1)) / 2;

			for (int p = 1; p < apply_to_size - 1; p++) {
				auto prev = find_by_value(p - 1);
				auto cur = find_by_value(p);
				auto next = find_by_value(p + 1);
				F[2 * mesh.map_node_numeration.at(apply_to[p]) + 1] += data[1] * std::abs((next - cur)) / 2;
				F[2 * mesh.map_node_numeration.at(apply_to[p]) + 1] += data[1] * std::abs((prev - cur)) / 2;
			}
			auto end0 = find_by_value(apply_to_size - 2);
			auto end1 = find_by_value(apply_to_size - 1);
			F[2 * mesh.map_node_numeration.at(apply_to[apply_to_size - 1]) + 1] += data[1] * std::abs((end0 - end1)) / 2;
		}
	*/}
}

void applyconstraints(const json& fc, std::vector<double>& K, const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols, std::vector<double>& F, const UnstructedMesh& mesh) {

	auto get_mesh_data = [&fc](std::string const& key, size_t* out_size = nullptr) -> char*
		{
			std::string data64 = fc["mesh"][key];
			char* data = (char*)malloc(data64.size() * sizeof(char));
			b64decode_host(data64.data(), data64.size(), data);
			size_t size = b64dsize_host(data64.size(), data64.data());

			if (out_size) *out_size = size;

			return data;
		};


	for (auto& constraint : fc["restraints"]) {
		int flag[6];
		int coordinate = -1;
		for (int i = 0; i < 6; i++) {
			flag[i] = constraint["flag"][i];
			if (i >= 3 && flag[i] != 0) {
				throw std::runtime_error("constraints at pivot. Not supported yet");
			}
			else if (coordinate == -1 && flag[i] == 1) {
				coordinate = i;
			}
			else if (coordinate != -1 && flag[i] != 0) {
				throw std::runtime_error("constraints at 2 or many coord. Not supported yet");
			}
		}
		size_t apply_to_size = constraint["apply_to_size"];
		auto get_constraint_data = [&constraint](std::string const& key, size_t* out_size = nullptr) -> char*
			{
				std::string data64 = constraint[key];
				char* data = (char*)malloc(data64.size() * sizeof(char));
				b64decode_host(data64.data(), data64.size(), data);
				size_t size = b64dsize_host(data64.size(), data64.data());

				if (out_size) *out_size = size;

				return data;
			};

		int* apply_to_raw = reinterpret_cast<int*>(get_constraint_data("apply_to"));
		std::vector<int> apply_to(apply_to_raw, apply_to_raw + apply_to_size);
		free(apply_to_raw);


		std::vector<double> data(6);
		for (int i = 0; i < 6; i++) {
			std::string data_b64 = constraint["data"][i];
			b64decode_host(data_b64.data(), data_b64.size(), reinterpret_cast<char*>(&(data[i])));
		}

		for (int& node : apply_to) {
			int row = mesh.map_node_numeration.at(node);
			//mesh.nodes[row];
			int n = mesh.nodes.size();
			for (int j = 0; j < n; j++) {
				for (int i = rows[j]; i < rows[j + 1]; i++) {
					if (flag[0] && cols[i] == row) {
						K[4 * i + 0] = 0;
						K[4 * i + 2] = 0;
					}
					if (flag[1] && cols[i] == row) {
						K[4 * i + 1] = 0;
						K[4 * i + 3] = 0;
					}
				}
			}

			for (int i = rows[row]; i < rows[row + 1]; i++) {
				if (flag[0]) {
					if (cols[i] == row) {
						K[4 * i + 0] = 1;
					}
					else {
						K[4 * i + 0] = 0;
					}
					K[4 * i + 1] = 0;
				}
				if (flag[1]) {
					K[4 * i + 2] = 0;
					if (cols[i] == row) {
						K[4 * i + 3] = 1;
					}
					else {
						K[4 * i + 3] = 0;
					}

				}
			}
			if (flag[0]) {
				F[2 * row + 0] = data[0];
			}
			if (flag[1]) {
				F[2 * row + 1] = data[1];
			}
		}
	}
}

void solve(const int& dim, const std::vector<double>& K, const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols, const std::vector<double>& F, std::vector<double>& x) {
	MKL_INT _iparm[64];
	void* _pt[64];

	// Setup Pardiso control parameters
	for (int i = 0; i < 64; i++) {
		_iparm[i] = 0;
	}

	_iparm[0] = 1;  /* No solver default */
	_iparm[1] = 3; /* Fill-in reordering from METIS */ // !!! = 0
	/* Numbers of processors, value of OMP_NUM_THREADS */
	_iparm[2] = 1;
	_iparm[3] = 0; /* No iterative-direct algorithm */
	_iparm[4] = 0; /* No user fill-in reducing permutation */
	_iparm[5] = 0; /* If =0 then write solution only into x. If =1 then the RightHandSide-array will replaced by the solution*/
	_iparm[6] = 0; /* Not in use */
	_iparm[7] = 2; /* Max numbers of iterative refinement steps */
	_iparm[8] = 0; /* Not in use */
	_iparm[11] = 0; /* Not in use */
	if (1) { // sym by default
		_iparm[9] = 8;
		_iparm[10] = 0; /* Disable scaling. Default for symmetric indefinite matrices. */
		_iparm[12] = 0; /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	}
	else {
		_iparm[9] = 13;
		_iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
		_iparm[12] = 1; /*1!!! Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
	}
	//_iparm[12] = 1; /*1!!! Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
	_iparm[13] = 0; /* Output: Number of perturbed pivots */
	_iparm[14] = 0; /* Not in use */
	_iparm[15] = 0; /* Not in use */
	_iparm[16] = 0; /* Not in use */
	_iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	_iparm[18] = -1; /* Output: Mflops for LU factorization */
	_iparm[19] = 0; /* Output: Numbers of CG Iterations */
	_iparm[26] = 1; /* Matrix Checker */
	if (1) // double by default
		_iparm[27] = 0;
	else
		_iparm[27] = 1;

	_iparm[59] = 0;

	_iparm[34] = 1; // zero-indexing
	_iparm[36] = dim; // bsr: block size
	for (int i = 0; i < 64; i++) {
		_pt[i] = 0;
	}



	MKL_INT n = rows.size() - 1;
	MKL_INT nnz = rows[n];

	const MKL_INT* h_RowsA = rows.data();
	const MKL_INT* h_ColsA = cols.data();
	const double* h_ValsA = K.data();
	x.resize(2*n);
	const double* h_b = F.data();
	double* h_x = x.data();
	MKL_INT nrhs = 1;
	double ddum = 0.0;
	MKL_INT maxfct = 1;
	MKL_INT msglvl = 0;
	MKL_INT mnum = 1;
	MKL_INT mtype = 11;
	MKL_INT idum = 0;
	MKL_INT phase = 11;
	MKL_INT error = 0;

	//phase11

	pardiso(&_pt[0], (MKL_INT*)&maxfct, (MKL_INT*)&mnum, (MKL_INT*)&mtype, (MKL_INT*)&phase, (MKL_INT*)&n,
		(void*)h_ValsA, (MKL_INT*)h_RowsA, (MKL_INT*)h_ColsA,
		(MKL_INT*)&idum, (MKL_INT*)&nrhs, (MKL_INT*)&_iparm[0], (MKL_INT*)&msglvl, &ddum, &ddum, (MKL_INT*)&error);


	//phase22
	phase = 22;
	pardiso((MKL_INT*)&_pt[0], (MKL_INT*)&maxfct, (MKL_INT*)&mnum, (MKL_INT*)&mtype, (MKL_INT*)&phase, (MKL_INT*)&n,
		(void*)h_ValsA, (MKL_INT*)h_RowsA, (MKL_INT*)h_ColsA,
		(MKL_INT*)&idum, (MKL_INT*)&nrhs, (MKL_INT*)&_iparm[0], (MKL_INT*)&msglvl, &ddum, &ddum, (MKL_INT*)&error);

	//phase33
	phase = 33;
	pardiso((MKL_INT*)&_pt[0], (MKL_INT*)&maxfct, (MKL_INT*)&mnum, (MKL_INT*)&mtype, (MKL_INT*)&phase, (MKL_INT*)&n,
		(void*)h_ValsA, (MKL_INT*)h_RowsA, (MKL_INT*)h_ColsA,
		(MKL_INT*)&idum, (MKL_INT*)&nrhs, (MKL_INT*)&_iparm[0], (MKL_INT*)&msglvl, (void*)h_b, (void*)h_x, (MKL_INT*)&error);

	//phase -1
	phase = -1;
	pardiso(&_pt[0], (MKL_INT*)&maxfct, (MKL_INT*)&mnum, (MKL_INT*)&mtype, (MKL_INT*)&phase, (MKL_INT*)&n,
		(void*)h_ValsA, (MKL_INT*)h_RowsA, (MKL_INT*)h_ColsA,
		(MKL_INT*)&idum, (MKL_INT*)&nrhs, (MKL_INT*)&_iparm[0], (MKL_INT*)&msglvl, (void*)h_b, (void*)h_x, (MKL_INT*)&error);
	mkl_free_buffers();
}

void resultants(const int& dim, material_t material, std::vector<double>& eps, std::vector<double>& sigma, const std::vector<double>& x, const UnstructedMesh& mesh, const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols) {
	if (dim != 2) {
		throw std::runtime_error("resultants: try to solve non 2D task");
	}
	int blocksize = dim;
	size_t n = rows.size() - 1;
	size_t nnz = rows[n];
	size_t elems_size = mesh.elemids.size();
	
	std::span<double> C(reinterpret_cast<double*>(malloc(nnz * sizeof(double))), nnz);
	for (int i = 0; i < nnz; i++) {
		C[i] = 0.0;
	}
	std::span<double> b(reinterpret_cast<double*>(malloc(6 * n * sizeof(double))), 6 * n);
	for (int i = 0; i < 6 * n; i++) {
		b[i] = 0.0;
	}

	// D matrix	
	int Ddim = 3;

	std::span<double> D(reinterpret_cast<double*>(malloc(Ddim * Ddim * sizeof(double))), Ddim * Ddim);

	// filling D
	fp::fp_t E = material.E;
	fp::fp_t nu = material.nu;

	D[0] = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
	D[1] = E * nu / ((1 + nu) * (1 - 2 * nu));
	D[2] = 0;
	D[3] = E * nu / ((1 + nu) * (1 - 2 * nu));
	D[4] = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
	D[5] = 0;
	D[6] = 0;
	D[7] = 0;
	D[8] = E / (2 * (1 + nu));

	// only if Tri
	int Brows = 3;
	int Bcols = 6;
	std::span<double> B(reinterpret_cast<double*>(malloc(Brows * Bcols * sizeof(double))), Brows * Bcols);
	std::span<double> u(reinterpret_cast<double*>(malloc(Bcols * sizeof(double))), Bcols);
	std::span<double> sigma_loc(reinterpret_cast<double*>(malloc(Brows * sizeof(double))), Brows);
	std::span<double> eps_loc(reinterpret_cast<double*>(malloc(Brows * sizeof(double))), Brows);
	for (int e = 0; e < elems_size; e++) {
		std::fill(sigma_loc.begin(), sigma_loc.end(), 0.0);
		std::fill(eps_loc.begin(), eps_loc.end(), 0.0);
		std::vector<int> ijk = { mesh.elems[mesh.nodes_per_elem[e] + 0] - 1,
									mesh.elems[mesh.nodes_per_elem[e] + 1] - 1,
									mesh.elems[mesh.nodes_per_elem[e] + 2] - 1 };
		//std::sort(ijk.begin(), ijk.end());

		int i = ijk[0];
		int j = ijk[1];
		int k = ijk[2];
		double J = std::abs((mesh.nodes[j].x - mesh.nodes[i].x) * (mesh.nodes[k].y - mesh.nodes[i].y) -
			(mesh.nodes[k].x - mesh.nodes[i].x) * (mesh.nodes[j].y - mesh.nodes[i].y));

		double S = J / 2.0;

		B[0] = mesh.nodes[j].y - mesh.nodes[k].y;
		B[1] = 0;
		B[2] = mesh.nodes[k].x - mesh.nodes[j].x;

		B[3] = 0;
		B[4] = mesh.nodes[k].x - mesh.nodes[j].x;
		B[5] = mesh.nodes[j].y - mesh.nodes[k].y;

		B[6] = mesh.nodes[k].y - mesh.nodes[i].y;
		B[7] = 0;
		B[8] = mesh.nodes[i].x - mesh.nodes[k].x;

		B[9] = 0;
		B[10] = mesh.nodes[i].x - mesh.nodes[k].x;
		B[11] = mesh.nodes[k].y - mesh.nodes[i].y;

		B[12] = mesh.nodes[i].y - mesh.nodes[j].y;
		B[13] = 0;
		B[14] = mesh.nodes[j].x - mesh.nodes[i].x;

		B[15] = 0;
		B[16] = mesh.nodes[j].x - mesh.nodes[i].x;
		B[17] = mesh.nodes[i].y - mesh.nodes[j].y;

		u[0] = x[2 * i + 0];
		u[1] = x[2 * i + 1];
		u[2] = x[2 * j + 0];
		u[3] = x[2 * j + 1];
		u[4] = x[2 * k + 0];
		u[5] = x[2 * k + 1];

		CBLAS_LAYOUT layout = CblasColMajor;
		CBLAS_TRANSPOSE nontrans = CblasNoTrans;
		CBLAS_TRANSPOSE trans = CblasTrans;
		const double alpha = 1.0;
		const double beta = 0.0;
		cblas_dgemv(layout, nontrans, Brows, Bcols, 1 / (2 * S), B.data(), Brows, u.data(), 1, beta, eps_loc.data(), 1);
		eps_loc[2] /= 2;
		cblas_dgemv(layout, nontrans, Ddim, Ddim, alpha, D.data(), Ddim, eps_loc.data(), 1, beta, sigma_loc.data(), 1);
		sigma_loc[2] *= 2;
		int i1 = -1;
		int j1 = -1;
		int k1 = -1;
		for (int q = rows[i]; q < rows[i + 1]; q++) {
			//zero or one indexing
			if (cols[q] == i) {
				i1 = q;
			}
			if (cols[q] == j) {
				j1 = q;
			}
			if (cols[q] == k) {
				k1 = q;
			}
		}
		C[i1] += J / 12.0;
		C[j1] += J / 24.0;
		C[k1] += J / 24.0;
		
		b[i + 0 * n] += J * eps_loc[0] / 6.0;
		b[i + 1 * n] += J * eps_loc[1] / 6.0;
		b[i + 2 * n] += J * eps_loc[2] / 6.0;

		b[i + 3 * n] += J * sigma_loc[0] / 6.0;
		b[i + 4 * n] += J * sigma_loc[1] / 6.0;
		b[i + 5 * n] += J * sigma_loc[2] / 6.0;
		
		int i2 = -1;
		int j2 = -1;
		int k2 = -1;
		for (int q = rows[j]; q < rows[j + 1]; q++) {
			//zero or one indexing
			if (cols[q] == i) {
				i2 = q;
			}
			if (cols[q] == j) {
				j2 = q;
			}
			if (cols[q] == k) {
				k2 = q;
			}
		}
		C[i2] += J / 24.0;
		C[j2] += J / 12.0;
		C[k2] += J / 24.0;
		
		b[j + 0 * n] += J * eps_loc[0] / 6.0;
		b[j + 1 * n] += J * eps_loc[1] / 6.0;
		b[j + 2 * n] += J * eps_loc[2] / 6.0;

		b[j + 3 * n] += J * sigma_loc[0] / 6.0;
		b[j + 4 * n] += J * sigma_loc[1] / 6.0;
		b[j + 5 * n] += J * sigma_loc[2] / 6.0;
		
		int i3 = -1;
		int j3 = -1;
		int k3 = -1;
		for (int q = rows[k]; q < rows[k + 1]; q++) {
			//zero or one indexing
			if (cols[q] == i) {
				i3 = q;
			}
			if (cols[q] == j) {
				j3 = q;
			}
			if (cols[q] == k) {
				k3 = q;
			}
		}
		C[i3] += J / 24.0;
		C[j3] += J / 24.0;
		C[k3] += J / 12.0;

		b[k + 0 * n] += J * eps_loc[0] / 6.0;
		b[k + 1 * n] += J * eps_loc[1] / 6.0;
		b[k + 2 * n] += J * eps_loc[2] / 6.0;
		
		b[k + 3 * n] += J * sigma_loc[0] / 6.0;
		b[k + 4 * n] += J * sigma_loc[1] / 6.0;
		b[k + 5 * n] += J * sigma_loc[2] / 6.0;

	}
	

	
	MKL_INT _iparm[64];
	void* _pt[64];

	// Setup Pardiso control parameters
	for (int i = 0; i < 64; i++) {
		_iparm[i] = 0;
	}

	_iparm[0] = 1;  /* No solver default */
	_iparm[1] = 3; /* Fill-in reordering from METIS */ // !!! = 0
	/* Numbers of processors, value of OMP_NUM_THREADS */
	_iparm[2] = 1;
	_iparm[3] = 0; /* No iterative-direct algorithm */
	_iparm[4] = 0; /* No user fill-in reducing permutation */
	_iparm[5] = 0; /* If =0 then write solution only into x. If =1 then the RightHandSide-array will replaced by the solution*/
	_iparm[6] = 0; /* Not in use */
	_iparm[7] = 2; /* Max numbers of iterative refinement steps */
	_iparm[8] = 0; /* Not in use */
	_iparm[11] = 0; /* Not in use */
	if (0) { // sym by default
		_iparm[9] = 8;
		_iparm[10] = 0; /* Disable scaling. Default for symmetric indefinite matrices. */
		_iparm[12] = 0; /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	}
	else {
		_iparm[9] = 13;
		_iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
		_iparm[12] = 1; /*1!!! Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
	}
	//_iparm[12] = 1; /*1!!! Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
	_iparm[13] = 0; /* Output: Number of perturbed pivots */
	_iparm[14] = 0; /* Not in use */
	_iparm[15] = 0; /* Not in use */
	_iparm[16] = 0; /* Not in use */
	_iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	_iparm[18] = -1; /* Output: Mflops for LU factorization */
	_iparm[19] = 0; /* Output: Numbers of CG Iterations */
	_iparm[26] = 1; /* Matrix Checker */
	if (1) // double by default
		_iparm[27] = 0;
	else
		_iparm[27] = 1;

	_iparm[59] = 0;

	_iparm[34] = 1; // zero-indexing
	_iparm[36] = 0; // bsr: block size
	for (int i = 0; i < 64; i++) {
		_pt[i] = 0;
	}



	const MKL_INT* h_RowsA = rows.data();
	const MKL_INT* h_ColsA = cols.data();
	const double* h_ValsA = C.data();

	const double* h_b = b.data();
	std::span<double> x_(reinterpret_cast<double*>(malloc(6 * n * sizeof(double))), 6 * n);
	double* h_x = x_.data();
	MKL_INT nrhs = 6;
	double ddum = 0.0;
	MKL_INT maxfct = 1;
	MKL_INT msglvl = 0;
	MKL_INT mnum = 1;
	MKL_INT mtype = 11;
	MKL_INT idum = 0;
	MKL_INT phase = 11;
	MKL_INT error = 0;

	//phase11

	pardiso(&_pt[0], (MKL_INT*)&maxfct, (MKL_INT*)&mnum, (MKL_INT*)&mtype, (MKL_INT*)&phase, (MKL_INT*)&n,
		(void*)h_ValsA, (MKL_INT*)h_RowsA, (MKL_INT*)h_ColsA,
		(MKL_INT*)&idum, (MKL_INT*)&nrhs, (MKL_INT*)&_iparm[0], (MKL_INT*)&msglvl, &ddum, &ddum, (MKL_INT*)&error);


	//phase22
	phase = 22;
	pardiso((MKL_INT*)&_pt[0], (MKL_INT*)&maxfct, (MKL_INT*)&mnum, (MKL_INT*)&mtype, (MKL_INT*)&phase, (MKL_INT*)&n,
		(void*)h_ValsA, (MKL_INT*)h_RowsA, (MKL_INT*)h_ColsA,
		(MKL_INT*)&idum, (MKL_INT*)&nrhs, (MKL_INT*)&_iparm[0], (MKL_INT*)&msglvl, &ddum, &ddum, (MKL_INT*)&error);

	//phase33
	phase = 33;
	pardiso((MKL_INT*)&_pt[0], (MKL_INT*)&maxfct, (MKL_INT*)&mnum, (MKL_INT*)&mtype, (MKL_INT*)&phase, (MKL_INT*)&n,
		(void*)h_ValsA, (MKL_INT*)h_RowsA, (MKL_INT*)h_ColsA,
		(MKL_INT*)&idum, (MKL_INT*)&nrhs, (MKL_INT*)&_iparm[0], (MKL_INT*)&msglvl, (void*)h_b, (void*)h_x, (MKL_INT*)&error);

	//phase -1
	phase = -1;
	pardiso(&_pt[0], (MKL_INT*)&maxfct, (MKL_INT*)&mnum, (MKL_INT*)&mtype, (MKL_INT*)&phase, (MKL_INT*)&n,
		(void*)h_ValsA, (MKL_INT*)h_RowsA, (MKL_INT*)h_ColsA,
		(MKL_INT*)&idum, (MKL_INT*)&nrhs, (MKL_INT*)&_iparm[0], (MKL_INT*)&msglvl, (void*)h_b, (void*)h_x, (MKL_INT*)&error);
	mkl_free_buffers();


	eps.resize(3 * n);
	sigma.resize(3 * n);
	for (int i = 0; i < n; i++) {
		eps[3 * i + 0] = h_x[i + 0 * n];
		eps[3 * i + 1] = h_x[i + 1 * n];
		eps[3 * i + 2] = h_x[i + 2 * n];

		sigma[3 * i + 0] = h_x[i + 3 * n];
		sigma[3 * i + 1] = h_x[i + 4 * n];
		sigma[3 * i + 2] = h_x[i + 5 * n];
	}
	//eps = std::vector(h_x, h_x + 1 * n);
	//sigma = std::vector(h_x + 3 * n, h_x + 6 * n);
	free(C.data());
	free(b.data());
	free(D.data());
	free(B.data());
	free(u.data());
	free(sigma_loc.data());

}