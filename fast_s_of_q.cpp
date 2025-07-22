#include <string>
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cassert>
#include <ctime>
#include <sstream>
#include <cctype>
#include <chrono>
#include <ranges>

#include "s_of_q.hpp"


/**
 * @brief Data structure for storing particle system information extracted from .vtf files.
 *
 * This structure holds the essential information about a particle system that has been
 * parsed from a Visual Trajectory Format (.vtf) file. It serves as a container for
 * system metadata and particle coordinates, providing a convenient way to pass this
 * data between functions for further analysis.
 *
 * @note This structure is designed to work with the extractPartPositions() function.
 *
 * @see extractPartPositions() for the function that populates this structure
 * @see calculate_s_of_q() for typical usage of the position data
 */
struct Extraction {
    /**
     * @brief Number of particles in the system.
     */
    int N_part;

    /**
     * @brief Number of bonds in the system.
     */
    int N_bonds;

    /**
     * @brief 2D vector containing particle positions.
     * 
     * Each inner vector contains the [x, y, z] coordinates of a single particle.
     * The outer vector index corresponds to the particle ID, and positions[i]
     * contains the 3D coordinates for particle i.
     * 
     * Structure: positions[particle_id][coordinate] where coordinate is 0=x, 1=y, 2=z
     */
    std::vector<std::vector<double>> positions;

    /**
     * @brief Vector containing simulation box dimensions.
     */
    std::vector<double> box_l;
};

/**
 * @brief Extracts particle positions and system information from a .vtf file.
 *
 * This function parses a Visual Trajectory Format (.vtf) file to extract particle
 * positions, system box dimensions, and counts of particles and bonds. The function
 * processes the file in two phases: first reading the header information (unitcell,
 * atoms, bonds) and then reading the timestep data containing particle coordinates.
 *
 * The .vtf file format structure expected:
 * - Header section with 'unitcell', 'atom', and 'bond' declarations
 * - Timestep section starting with 'timestep' keyword
 * - Position data with particle ID followed by x, y, z coordinates
 *
 * @param vtfFilePath The path to the .vtf file to be parsed.
 *
 * @return Extraction A struct containing:
 *         - N_part: Number of particles in the system
 *         - N_bonds: Number of bonds in the system  
 *         - positions: 2D vector containing [x, y, z] coordinates for each particle
 *         - box_l: Vector containing the box dimensions
 *
 * @throws std::runtime_error If the file cannot be opened.
 *
 * @note The function assumes:
 *       - Lines starting with 'u' contain unitcell information
 *       - Lines starting with 'a' declare atoms (counted for N_part)
 *       - Lines starting with 'b' declare bonds (counted for N_bonds)
 *       - Lines starting with 't' indicate the start of timestep data
 *       - Position lines contain particle ID as first value, followed by coordinates
 *       - Case-insensitive parsing (converts to lowercase)
 *
 * @warning The function performs an assertion check to ensure consistency between
 *          the number of declared particles and the number of position entries read.
 *          The program will terminate if this check fails.
 *
 * @example
 * try {
 *     Extraction data = extractPartPositions("trajectory.vtf");
 *     std::cout << "Particles: " << data.N_part << std::endl;
 *     std::cout << "Box length: " << data.box_l[0] << std::endl;
 * } catch (const std::runtime_error& e) {
 *     std::cerr << "Error: " << e.what() << std::endl;
 * }
 */
Extraction extractPartPositions(const std::string& vtfFilePath) {
    int N_part = 0;
    int N_bonds = 0;
    bool pos_bool = false;
    std::vector<std::vector<double>> positions;
    std::vector<double> box_l;
    int part_id = -1;

    std::ifstream file(vtfFilePath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::string line;
    while (std::getline(file, line)) {
        if (!pos_bool) {
            // Trim and convert to lowercase
            std::transform(line.begin(), line.end(), line.begin(), ::tolower);
            line.erase(0, line.find_first_not_of(" \t"));
            
            if (line[0] == 'u') {
                std::istringstream iss(line.substr(9));
                double value;
                while (iss >> value) {
                    box_l.push_back(value);
                }
            } else if (line[0] == 'a') {
                N_part++;
            } else if (line[0] == 'b') {
                N_bonds++;
            } else if (line[0] == 't') {
                pos_bool = true;
            }
        } else {
            std::istringstream iss(line);
            std::vector<double> pos;
            double value;
            
            while (iss >> value) {
                pos.push_back(value);
            }
            
            if (!pos.empty()) {
                part_id = static_cast<int>(pos[0]);
                pos.erase(pos.begin()); // Remove the first element (particle ID)
                positions.push_back(pos);
            }
        }
    }
    
    assert(part_id == N_part - 1 && "Inconsistency between number of particles and number of positions.");

    return {N_part, N_bonds, positions, box_l};
}

/**
 * @brief Generates linearly spaced wave vector magnitudes based on the fundamental wave vector.
 *
 * Creates q values that are integer multiples of the fundamental wave vector q₀ = 2π/L,
 * where L is the box length. This approach ensures that the wave vectors are compatible
 * with the periodic boundary conditions of the simulation box.
 *
 * @param num_q Number of q values to generate
 * @param q_min Minimum q value threshold - used to determine the starting multiple of q₀
 * @param box_length Length of the cubic simulation box
 * @return std::vector<float> Linearly spaced q values as multiples of q₀
 *
 * @note The function:
 *       - Calculates q₀ = 2π/box_length
 *       - Determines n_min = max(2, floor(q_min/q₀)) as the starting multiplier
 *       - Generates q values: q[i] = q₀ * (n_min + i) for i = 0, 1, ..., num_q-1
 *       - Uses C++20 ranges with std::views::iota and std::views::transform
 *
 * @example For box_length=10, q_min=1.0, num_q=5:
 *          q₀ ≈ 0.628, n_min=2, generates q ≈ [1.26, 1.88, 2.51, 3.14, 3.77]
 */
std::vector<float> calculate_qs(int& num_q, float& q_min, float& box_length) {
    std::vector<float> q(num_q);
    float q_0 = 2 * M_PI / box_length; // q_0 = 2 * pi / box_length
    int val = static_cast<int>(std::floor(q_min / q_0));
    int n_min = val >= 2 ? val : 2; // Ensure at least 2
    
    // Using std::views::transform with ranges
    auto indices = std::views::iota(0, num_q);
    auto q_values = indices | std::views::transform([q_0, n_min](int i) {
        return q_0 * (n_min + i);
    });
    
    // Copy the transformed values to the output vector
    std::ranges::copy(q_values, q.begin());
    
    return q;
}

/**
 * @brief Main program for calculating structure factor S(q) from VTF trajectory files.
 *
 * Reads particle positions from a .vtf file, computes the structure factor using
 * the Debye formula, and outputs results to specified files.
 *
 * @param argc Number of command line arguments (5-7 expected)
 * @param argv Command line arguments:
 *             [1] path_to_traj - Path to .vtf trajectory file
 *             [2] num_q - Number of q values to calculate
 *             [3] q_min - Minimum q value
 *             [4] out_file - Output file for S(q) data (optional, default: "s_of_q.txt")
 *             [5] info_file - System info file (optional, default: "info.txt")
 *
 * @return 0 on success, 1 on file I/O error
 *
 * @note q values are generated as multiples of the fundamental wave vector q₀ = 2π/L,
 *       where L is the box length. The minimum q value is used to determine the starting
 *       multiple for generating the q values. If q_min is less than 2q₀, it will start
 *       from the first valid multiple, namely 2q₀.
 * @note Outputs two files:
 *       - S(q) data: Two-column format with q values and corresponding S(q)
 *       - System info: Metadata including particle count, timing, and parameters
 *
 * @example ./fast_s_of_q trajectory.vtf 50 0.1 10.0 results.txt info.txt
 */
int main(int argc, char *argv[]) {

    // Read command line arguments
    std::cout << "Reading command line arguments..." << std::endl;
    // --path_to_traj: Path to the .vtf file
    // --num_q: Number of q values
    // --q_min: Minimum q value
    std::string path_to_traj(argv[1]);
    int num_q = std::stoi(argv[2]);
    float q_min = std::stof(argv[3]);

    const std::string out_file = argc >= 5 ? argv[4] : "s_of_q.txt";
    const std::string info_file = argc == 6 ? argv[5] : "info.txt";

    std::cout << "Reading data from vtf file: " << path_to_traj << std::endl;
    
    std::cout << "Extracting particle positions..." << std::endl;
    Extraction data = extractPartPositions(path_to_traj);
    
    std::cout << "Defining variables..." << std::endl;
    float box_length = data.box_l[0];
    std::vector<float> q_vec = calculate_qs(num_q, q_min, box_length);
    std::vector<float> pos_x;
    std::vector<float> pos_y;
    std::vector<float> pos_z;
    
    std::cout << "Initialize positions.." << std::endl;
    for (const auto& pos : data.positions) {
        pos_x.push_back(pos[0]);
        pos_y.push_back(pos[1]);
        pos_z.push_back(pos[2]);
    }
    
    std::cout << "Calculate S(q)..." << std::endl;
    auto const start = std::chrono::high_resolution_clock::now();
    std::vector<float> s_of_q =  calculate_s_of_q(pos_x, pos_y, pos_z, q_vec, box_length);
    auto const end = std::chrono::high_resolution_clock::now();
    
    std::cout << "S(q) calculated." << std::endl;
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    
    std::cout << "Writing S(q) to file..." << std::endl;
    std::ofstream output_file(out_file);
    if (!output_file) {
        std::cerr << "Error opening output file." << std::endl;
        return 1;
    }
    output_file << "# q S(q)\n";
    for (size_t i = 0; i < q_vec.size(); ++i) {
        output_file << q_vec[i] << " " << s_of_q[i] << "\n";
    }
    output_file.close();

    std::cout << "Writing system information to file..." << std::endl;

    std::ofstream system_info_file(info_file);
    if (!system_info_file) {
        std::cerr << "Error opening system info file." << std::endl;
        return 1;
    }
    system_info_file << "vft file: " << path_to_traj << "\n";
    system_info_file << "Number of particles: " << data.N_part << "\n";
    system_info_file << "Number of bonds: " << data.N_bonds << "\n";
    system_info_file << "Box length: " << box_length << "\n";
    system_info_file << "Number of q values: " << num_q << "\n";
    system_info_file << "Minimum q value: " << q_vec[0] << "\n";
    system_info_file << "Maximum q value: " << q_vec.back() << "\n";
    system_info_file <<  "calculation time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s   OR   " \
	    <<  std::chrono::duration_cast<std::chrono::minutes>(end - start).count() << " min" << "\n";
    system_info_file.close();
    return 0;
}
