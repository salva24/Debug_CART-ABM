/**
 * @file Macros.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Utility macros for the cancer simulation system
 *
 * @details This file defines several utility macros used throughout the cancer simulation
 * codebase. These macros provide efficient implementations of commonly used operations,
 * particularly for spatial indexing and debugging.
 * 
 * In cancer research, efficient spatial indexing is critical for:
 * - Rapid lookup of cells in 3D tissue environments
 * - Efficient calculation of concentration gradients across tumor regions
 * - Performance optimization in large-scale tumor simulations
 * - Accurate representation of spatial heterogeneity in tumors
 * 
 * Inbound dependencies vector, cmath, string
 *
 * Outbound dependencies Used throughout the codebase for spatial indexing
 *
 * @usage Included in files that need efficient spatial indexing or debugging
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef _MACROS_H
#define _MACROS_H

#include <vector>
#include <cmath>
#include <string>

/**
 * @brief Macro for calculating linear voxel index from 3D coordinates
 * @details Converts 3D coordinates (i,j,k) to a linear index for accessing voxel data
 * in a 1D array. Uses the global coordinate vectors to determine dimensions.
 * 
 * In cancer research, this spatial indexing is essential for:
 * - Mapping biochemical gradients across tumor microenvironments
 * - Locating cells in specific tumor regions (e.g., hypoxic core vs. proliferative rim)
 * - Analyzing spatial patterns in heterogeneous tumors
 * - Efficient substrate diffusion calculations
 */
#define g_indice_de_voxel(i,j,k) (k*coordenadas_y.size() + j )*coordenadas_x.size() + i


/**
 * @brief Macro for calculating linear voxel index from 3D coordinates using a grid reference
 * @details Similar to g_indice_de_voxel but uses the grid (mgrilla) member variable's
 * coordinate vectors to determine dimensions. Used when the grid is a class member.
 * 
 * In cancer modeling, this enables:
 * - Object-oriented access to spatial data in tumor microenvironments
 * - Consistent indexing across different grid-based components
 * - Efficient spatial queries in complex tissue architectures
 */
#define fg_indice_de_voxel(i,j,k) (k*mgrilla.coordenadas_y.size() + j )*mgrilla.coordenadas_x.size() + i

///*
//Para debuggear el código
/**
 * @brief Debug information macro
 * @details Outputs file name, line number, and a custom message to stderr.
 * Used for debugging the simulation code.
 * 
 * In cancer research software development, this supports:
 * - Validation of model behavior against experimental data
 * - Identification of issues in complex multi-scale cancer simulations
 * - Tracking execution flow in parallel processing of tumor regions
 * - Debugging of rare events in stochastic cancer models
 */
#define INFO(msg) \
    fprintf(stderr, "info: %s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, "%s\n", msg);
//*/

      
#endif
