#pragma once
#include <ranges>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <cmath>
#include <numeric>
#include <omp.h> // Include OpenMP header

/**
 * @brief Applies the minimum image convention to a distance.
 *
 * This function adjusts a given distance to ensure it lies within the range
 * [-box_length/2, box_length/2) by applying the minimum image convention.
 *
 * @tparam T The type of the distance and box length (e.g., float, double).
 * @param dist The distance to be adjusted.
 * @param box_length The length of the simulation box.
 * @return The adjusted distance within the range [-box_length/2, box_length/2).
 */
template <typename T>
T minimum_image_convention(T dist, T box_length) {
    return dist - box_length * std::round(dist / box_length);
}   


/**
 * @brief Applies the minimum image convention to a distance and squares the result.
 *
 * @tparam T The type of the distance and box length (e.g., float, double).
 * @param dist The distance to be adjusted and squared.
 * @param box_length The length of the simulation box.
 * @return The squared adjusted distance within the range [-box_length/2, box_length/2).
 */
template <typename T>
T mic_and_sqr(T dist, T box_length) {
    auto const mdist = minimum_image_convention(dist, box_length);
    return mdist * mdist;
}

/**
 * @brief Computes the outer product of two vectors.
 *
 * This function calculates the outer product (also known as tensor product) of two 
 * vectors, resulting in a matrix where each element (i,j) is the product of the 
 * i-th element of the first vector and the j-th element of the second vector.
 *
 * @tparam T The type of elements in both input vectors (e.g., float, double, int).
 * @param vec1 The first input vector (left operand).
 * @param vec2 The second input vector (right operand).
 *
 * @return std::vector<std::vector<T>> A 2D matrix represented as a vector of vectors,
 *         where result[i][j] = vec1[i] * vec2[j].
 *         The resulting matrix has dimensions vec1.size() × vec2.size().
 *
 * @note The function creates a new matrix and does not modify the input vectors.
 *       For large vectors, consider memory usage as the result matrix size is
 *       O(vec1.size() × vec2.size()).
 *
 * @example
 * std::vector<float> a = {1, 2, 3};
 * std::vector<float> b = {4, 5};
 * auto result = outer_product(a, b);
 * // result = {{4, 5}, {8, 10}, {12, 15}}
 *
 * @complexity O(N × M) where N = vec1.size() and M = vec2.size().
 */
template <typename T>
std::vector<std::vector<T>> outer_product(const std::vector<T>& vec1, const std::vector<T>& vec2) {
    std::vector<std::vector<T>> result(vec1.size(), std::vector<T>(vec2.size()));

    // Use std::transform to fill the result matrix
    // Each row corresponds to an element of vec1, and each column corresponds to an element of vec2
    std::transform(vec1.cbegin(), vec1.cend(), result.begin(), [&vec2](const T val1) {
        std::vector<T> row(vec2.size());
        std::transform(vec2.cbegin(), vec2.cend(), row.begin(), [val1](const T val2) {
            return val1 * val2;
        });
        return row;
    });
    return result;
}

template <typename T>
T iso_wave(T q, T r) {
    T qr = q * r;
    return std::sin(qr) / qr;
}

/**
 * @brief Computes the isotropic wave function for structure factor calculations.
 *
 * Calculates 2*sin(qr)/qr, commonly used in the Debye formula for computing
 * structure factors in scattering analysis.
 *
 * @tparam T Numeric type (float, double)
 * @param qr Product of wave vector magnitude and distance
 * @return The isotropic wave function value
 *
 * @note Factor of 2 accounts for symmetry in pair summation
 */
template <typename T>
T iso_wave(T qr) {
    return 2 * std::sin(qr) / qr;
}


/**
 * @brief Computes the outer product of two vectors with isotropic wave function applied.
 *
 * This function calculates a specialized outer product where instead of simple multiplication,
 * it applies the isotropic wave function iso_wave to the product of each element pair and
 * then sums the results for each row.
 *
 * The implementation uses C++20 ranges and views with std::views::transform and 
 * std::ranges::fold_left for efficient functional programming style computation.
 * The result is a lazy-evaluated view that can be materialized when needed.
 *
 * @tparam T The type of elements in both input vectors (e.g., float, double).
 * @param vec1 The first input vector (typically wave vector magnitudes q).
 * @param vec2 The second input vector (typically distance values r).
 *
 * @return auto A ranges view that yields T values, where each element corresponds to
 *         the sum of iso_wave(vec1[i] * vec2[j]) for all j, for each i.
 *         The resulting view has the same size as vec1.
 *
 * @note This function:
 *       - Uses lazy evaluation through std::views::transform
 *       - Applies iso_wave function to each product vec1[i] * vec2[j]
 *       - Sums all iso_wave results for each element of vec1
 *       - Requires C++20 ranges support
 *
 * @example
 * std::vector<float> q = {1.0, 2.0};
 * std::vector<float> r = {0.5, 1.0, 1.5};
 * auto result = outer_product_iso_wave(q, r);
 * // result[0] = iso_wave(1.0*0.5) + iso_wave(1.0*1.0) + iso_wave(1.0*1.5)
 * // result[1] = iso_wave(2.0*0.5) + iso_wave(2.0*1.0) + iso_wave(2.0*1.5)
 *
 * @complexity O(N × M) where N = vec1.size() and M = vec2.size().
 *
 * @see iso_wave() for the isotropic wave function implementation
 * @see calculate_s_of_q() for usage in structure factor calculations
 */
template <typename T>
auto outer_product_iso_wave(const std::vector<T>& vec1, const std::vector<T>& vec2) {
    return vec1 | std::views::transform([&vec2](const T val1) {
        return std::ranges::fold_left(vec2, T{}, [val1](const T sum, const T val2) {
            return sum + iso_wave<T>(val1 * val2);
        });
    });
}

// template <std::ranges::range Range1, std::ranges::range Range2>
// auto outer_product_iso_wave(Range1& vec1, const Range2& vec2) {
//     // print length of vec1 and vec2
//     std::cout << "vec1 size: " << std::ranges::size(vec1) << ", vec2 size: " << std::ranges::size(vec2) << std::endl;
//     return vec1 | std::views::transform([&vec2](const auto val1) {
//         return std::ranges::fold_left(vec2, std::ranges::range_value_t<Range2>{}, [val1](const auto sum, const auto val2) {
//             return sum + iso_wave(val1 * val2);
//         });
//     });
//     // return std::ranges::transform(vec1, std::begin(vec1), [&vec2](const auto val1) {
//     //     return std::ranges::fold_left(vec2, std::ranges::range_value_t<Range2>{}, [val1](const auto sum, const auto val2) {
//     //         return sum + iso_wave(val1 * val2);
//     //     });
//     // });
// }


template <typename T>
void print_vector(std::vector<T> v) {
    for (T i : v)
        std::cout << i << ' ';
    std::cout << '\n';
}

template <typename T>
void initialize_dest(const T &pos, T &dest){
    std::copy(pos.begin(), pos.end(), dest.begin());
    std::copy(pos.begin(), pos.end(), dest.begin()+pos.size());
}

#pragma omp declare reduction(vec_float_plus : std::vector<float> : \
    std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<float>())) \
initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))


/**
 * @brief Calculates the structure factor S(q) for a given set of particle positions.
 *
 * This function computes the structure factor S(q) using the Debye formula for a system
 * of particles in a periodic box. The structure factor is calculated by summing over all
 * pairs of particles and applying the isotropic wave vector function to the product q*r,
 * where q is the wave vector and r is the inter-particle distance.
 *
 * The calculation is parallelized using OpenMP with a custom reduction for vector addition.
 * Each thread processes different iterations of the main loop independently, and the
 * partial results are combined using the vec_float_plus reduction.
 *
 * @param pos_x Vector containing the x-coordinates of all particles.
 * @param pos_y Vector containing the y-coordinates of all particles.
 * @param pos_z Vector containing the z-coordinates of all particles.
 * @param q Vector containing the wave vector magnitudes for which S(q) is calculated.
 * @param box_length The length of the cubic simulation box (assumes cubic box).
 *
 * @return std::vector<float> Vector containing the structure factor values S(q) 
 *         corresponding to each q value in the input vector.
 *
 * @note The function assumes:
 *       - The simulation box is cubic with equal dimensions
 *       - OpenMP is available for parallelization
 *
 * @warning This function requires the vec_float_plus custom reduction to be declared
 *          before use. The reduction handles the parallel accumulation of partial
 *          structure factor contributions.
 *
 * @complexity O(N² * M) where N is the number of particles and M is the number of q values.
 */
std::vector<float> calculate_s_of_q(const std::vector<float> &pos_x,
                                    const std::vector<float> &pos_y,
                                    const std::vector<float> &pos_z,
                                    const std::vector<float> &q,
                                    const float &box_length) {

    size_t N_part = pos_x.size();

    std::vector<float> s_of_q(q.size(), 0.0f);

    #pragma omp parallel reduction(vec_float_plus : s_of_q)
    {
        std::vector<float> diff(pos_x.size()-1);
        std::vector<float> r_diff(diff.size());

        #pragma omp for schedule(dynamic)
        for (size_t i = 1; i<pos_x.size(); i++) {
            auto const it_rdiff_end = r_diff.end()-i+1;
            std::fill(r_diff.begin(), it_rdiff_end, 0.0f);
            
            for (const auto &pos: {pos_x, pos_y, pos_z}) {

                // Difference between the vectors  dx = x2 - x1
                std::transform(pos.cbegin(), pos.cend()-i, pos.cbegin()+i, diff.begin(), std::minus<float>());
                
                auto const it_diff_end = diff.end()-i+1;
                // Apply the minimum image convention and square the result:  dx^2
                // and sum the squared differences: vector of dr^2 = dx^2 + dy^2 + dz^2
                std::transform(diff.begin(), it_diff_end, r_diff.begin(), r_diff.begin(), [box_length](float x, float sum) {
                    return sum + x * x;
                });
            }

            // Root the sum of the squared differences: vector of dr = sqrt(dx^2 + dy^2 + dz^2)
            std::transform(r_diff.begin(), it_rdiff_end, r_diff.begin(), [](float x) { return std::sqrt(x); });

            // Compute the product q*r; apply the isotropic wave vector
            // and sum the isotropic wave vector
            std::vector<float> r_diff_subset(r_diff.begin(), it_rdiff_end);
            auto const s_of_q_partial = outer_product_iso_wave(q, r_diff_subset);

            // Sum of partial contributions
            std::transform(s_of_q.begin(), s_of_q.end(), s_of_q_partial.begin(), s_of_q.begin(), std::plus<float>());
        }
    }

    std::for_each(s_of_q.begin(), s_of_q.end(), [N_part](float &x) { x = (x + N_part) / N_part; });

    return s_of_q;
}
