/*
!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
*/

// Parallel STL implementation
//
// sorts coordinate arrays to create global addressing array ibool/iglob
//
// infos for setup:
//  gcc:
//    https://solarianprogrammer.com/2019/05/09/cpp-17-stl-parallel-algorithms-gcc-intel-tbb-linux-macos/
//  intel oneAPI Threading Building Blocks:
//    https://software.intel.com/content/www/us/en/develop/tools/oneapi.html
//    https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onetbb.html
//
// example compilation with gcc:
// > g++ -std=c++17 -c sort_array_coordinates_c.cpp -I./setup
//
// original contribution by Dmitry Alexeev

#if defined(USE_PARALLEL_STL_SORTING)

#include "config.h"

#include <algorithm>
#include <execution>
#include <numeric>
#include <cmath>

extern "C"
{
    void FC_FUNC_(sort_array_coordinates_c,SORT_ARRAY_COORDINATES_C)
        (int *_n, const double *x, const double *y, const double *z, int *iglob, double *_epsilon, int *nglob);
}

struct Point
{
    int id;
    double x,y,z;
};

void FC_FUNC_(sort_array_coordinates_c,SORT_ARRAY_COORDINATES_C)
    (int *_n, const double *x, const double *y, const double *z, int *iglob, double *_epsilon, int *nglob)
{
    const int n = *_n;
    const double epsilon = *_epsilon;
    std::vector<Point> points(n);

    std::vector<int> sequence(n), segments(n+1);

    for (int i=0; i<n; i++)
        points[i] = {i, x[i], y[i], z[i]};

    auto cmp_x = [] (const auto& el1, const auto el2) {
        return el1.x < el2.x;
    };
    auto cmp_y = [] (const auto& el1, const auto el2) {
        return el1.y < el2.y;
    };
    auto cmp_z = [] (const auto& el1, const auto el2) {
        return el1.z < el2.z;
    };

    // Sort by x
    std::sort(std::execution::par, points.begin(), points.end(), cmp_x);

    // Find continuous
    std::iota(sequence.begin(), sequence.end(), 0);
    segments[0] = 0;
    auto last_seg = std::copy_if(std::execution::par_unseq,
        std::next(sequence.begin()), sequence.end(), std::next(segments.begin()),
        [epsilon, &points] (int id) {
            return std::fabs(points[id].x-points[id-1].x) > epsilon;
        }
    );
    *last_seg = n;

    // Sort by y
    #pragma omp parallel for
    for (int i=0; i < last_seg-segments.begin(); i++)
        if (segments[i+1] - segments[i] > 1)
            std::sort(points.begin() + segments[i], points.begin() + segments[i+1], cmp_y);

    // Again find continuous
    last_seg = std::copy_if(std::execution::par_unseq,
        std::next(sequence.begin()), sequence.end(), std::next(segments.begin()),
        [epsilon, &points] (int id) {
            return std::fabs(points[id].x-points[id-1].x) > epsilon ||
                   std::fabs(points[id].y-points[id-1].y) > epsilon;
        }
    );
    *last_seg = n;

    // Sort by z
    #pragma omp parallel for
    for (int i=0; i < last_seg-segments.begin(); i++)
        if (segments[i+1] - segments[i] > 1)
            std::sort(points.begin() + segments[i], points.begin() + segments[i+1], cmp_z);

    // Create iglob mapping
    int ig = 1;
    iglob[ points[0].id ] = ig;
    for (int i=1; i<n; i++)
    {
        if (std::fabs(points[i].x-points[i-1].x) > epsilon ||
            std::fabs(points[i].y-points[i-1].y) > epsilon ||
            std::fabs(points[i].z-points[i-1].z) > epsilon)
        {
            ig++;
        }
        iglob[ points[i].id ] = ig;
    }
    *nglob = ig;
}

#endif
