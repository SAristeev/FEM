#include "FEM.hpp"
#include "vtu_writer/include/vtu_writer.h"
int main(int argc, char* argv[]) {

    std::unordered_map<std::string, std::string>  parsed_params;//in the pair {key,param} param may be empty

    for (int pos = 1; pos < argc; ++pos) {
        if (argv[pos][0] == '-') {//key is found
            std::string key(argv[pos]);

            if (pos + 1 < argc && argv[pos + 1][0] != '-') {//value is found
                std::string value(argv[pos + 1]);

                parsed_params.insert(std::make_pair(key, value));
                ++pos;
                continue;
            }
            parsed_params.insert(std::make_pair(key, std::string()));
        }
    }
    std::string FileName;
    {
        const auto it = parsed_params.find("--fc");
        if (it != parsed_params.end() && !it->second.empty()) {
            FileName = it->second;
            parsed_params.erase(it);
        }
        else {
            printf("ERROR: Path to fc is not provided!\n");
        }

    }
    
    std::ifstream fc_file(FileName, std::ios::in);
    if (!fc_file) {
        throw std::runtime_error("cannot open fc file: " + FileName);
    }
    auto fc = nlohmann::json::parse(fc_file);
    fc_file.close();

    int dim;
    read_dimensions(fc, dim);

    // read materials
    std::map<int, unsigned char> matid_threshold_map;
    std::map<int, unsigned char> block_threshold_map;
    std::vector<material_t> materials;


    read_materials(fc, materials, matid_threshold_map);

    std::vector<int> blocks;
    std::vector<unsigned char> thresholds;

    read_blocks(fc, matid_threshold_map, blocks, thresholds, block_threshold_map);

    UnstructedMesh mesh;
    read_mesh(fc, mesh);

    std::vector<MKL_INT> rows;
    std::vector<MKL_INT> cols;
    buildFullGlobalMatrixStruct(mesh, rows, cols);

    std::vector<double> K;
    buildFullGlobalMatrix(dim, K, materials[0], mesh, rows, cols);

    std::vector<double> F;
    createLoads(dim, fc, F, mesh);
    applyconstraints(fc, K, rows, cols, F, mesh);
    std::vector<double> x;
    solve(dim, K, rows, cols, F, x);
    if (0){
        std::cout << "K" << std::endl;
        for (int i = 0; i < K.size(); ++i) {
        
            std::cout << std::scientific << K[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "cols" << std::endl;    
        for (int i = 0; i < cols.size(); ++i) {
            std::cout << cols[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "rows" << std::endl;
        for (int i = 0; i < rows.size(); ++i) {
            std::cout << rows[i] << ", ";
        }
        std::cout << std::endl;
    }
    if (0) {
        for (int i = 0; i < F.size(); ++i) {
            std::cout << "f(" << i << ") = " << F[i] << ", x(" << i << ") = "  << x[i] << "\n";
        }
    }
    std::vector<double> eps;
    std::vector<double> sigma;
    std::vector<double> eps_e;
    std::vector<double> sigma_e;
    std::vector<double> test;
    int e = mesh.elemids.size();
    test.resize(2 * e);
    for (int i = 0; i < 2 * e; i++) {
        if (i < e) {
            test[i] = 1;
        }
        else {
            test[i] = 2;
        }
    }
    resultants(dim, materials[0], eps, sigma, x, mesh, rows, cols);
    {
        std::vector<vtu::spatial_data_t> p_data;
        p_data.emplace_back(x.data(), "Displacement", 2);
        p_data.emplace_back(F.data(), "External Force", 2);
        std::vector<std::string> components = { "XX", "YY", "XY" };
        p_data.emplace_back(sigma.data(), "Stress", components);
        p_data.emplace_back(eps.data(), "Strain", components);

        std::vector<vtu::spatial_data_t> c_data;

        vtu::writer_t writer("abc.vtu", true);
        auto p_cords = std::span<double>{ (double*)mesh.nodes.data(), mesh.nodes.size() * 3 };
        std::vector<int> conn(mesh.elems.size());
        for (int i = 0; i < mesh.elems.size(); i++) {
            conn[i] = mesh.map_node_numeration[mesh.elems[i]];
        }
        writer.write(p_cords, conn, mesh.elem_type, p_data, c_data);
    }
    return 0;
}