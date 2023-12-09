#include "FEM.hpp"
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
            printf("ERROR: Path to mesh is not provided!\n");
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
    std::vector<double> sigma;
    std::cout << "******************************************************88";
    resultants(dim, materials[0], sigma, x, mesh, rows, cols);
    
}