/**
 * @file Voxel.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Defines the Voxel class, which represents discrete spatial elements in the microenvironment
 *
 * @details
 * This file contains the definition of the Voxel class, which represents a single volume element 
 * in the discretized spatial domain of the simulation. Voxels are fundamental building blocks for
 * modeling spatial heterogeneity in cancer, allowing for the representation of localized
 * microenvironmental conditions (such as hypoxia, nutrient gradients, and acidosis) that strongly
 * influence tumor development and therapeutic response.
 * 
 * In cancer research, spatial discretization via voxels enables:
 * - Modeling of localized regions with distinct microenvironmental properties
 * - Tracking spatial distribution of cancer cells and their heterogeneous phenotypes
 * - Simulating diffusion of nutrients, growth factors, and therapeutic agents
 * - Analyzing tumor spatial architecture and invasion patterns
 * 
 * Inbound Dependencies: Vector.h, iostream
 *
 * Outbound Dependencies: Microambiente.h, Grillado.h
 *
 * @licence: Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __VOXEL_H__
#define __VOXEL_H__

#include "Vector.h"
#include <iostream>

/**
 * @brief Basic spatial unit for microenvironment discretization in cancer modeling
 * @details
 * The Voxel class represents a discrete cubic volume in 3D space, forming the basic unit of
 * spatial discretization for the microenvironment in cancer simulations. Each voxel has a defined
 * position (center), volume, and unique index for identification and accessing associated data.
 * 
 * In cancer research, voxels enable modeling of the tumor microenvironment's spatial heterogeneity,
 * which is critical for understanding cancer progression, as different regions within a tumor
 * experience varying conditions (oxygen, nutrients, pH) that drive diverse cellular behaviors
 * and therapeutic responses.
 * 
 * The voxel implementation allows for:
 * - Spatial tracking of biochemical substrate concentrations relevant to cancer
 * - Correlation of cell positions with local environmental conditions
 * - Computation of gradients that influence cell migration and polarization
 * - Mapping of therapy distribution for treatment response modeling
 */
class Voxel
{

 private:
	/**
	 * @brief Friend operator for stream output of voxel data
	 * @param os Output stream
	 * @param v Voxel to output
	 * @return Modified output stream
	 * @details
	 * Enables formatted output of voxel data for visualization and debugging.
	 * In cancer research, this facilitates data analysis of spatial aspects of
	 * the tumor microenvironment.
	 */
	friend std::ostream& operator<<(std::ostream& os, const Voxel& v); 

 public:
	/**
	 * @brief Default constructor
	 * @details
	 * Initializes a voxel with default settings: index 0, volume 10^3, and not a Dirichlet boundary.
	 * In cancer modeling, new voxels are typically initialized before being configured with
	 * specific parameters for the simulation domain.
	 */
	Voxel(); 
	
	/**
	 * @brief Unique identifier for the voxel within the simulation domain
	 * @details
	 * Index is used to access voxel-specific data in substrate and gradient vectors.
	 * In cancer research, this indexing system enables efficient lookup of microenvironmental
	 * conditions at specific spatial locations within the tumor or surrounding tissue.
	 */
	int indice;

	/**
	 * @brief Volume of the voxel in cubic micrometers
	 * @details
	 * Determines the spatial resolution of the microenvironment discretization.
	 * In cancer modeling, voxel volume affects the granularity of spatial heterogeneity
	 * representation, with smaller volumes providing higher resolution but at increased
	 * computational cost.
	 */
	double volumen;
	
	/**
	 * @brief 3D coordinates of the voxel's center position
	 * @details
	 * Defines the spatial location of the voxel in the simulation domain.
	 * In cancer research, the position enables spatial correlation between cells
	 * and their local microenvironment, allowing for analysis of position-dependent
	 * phenomena like invasion front dynamics and therapy penetration.
	 */
	Vector centro;
	
	/**
	 * @brief Flag indicating if this voxel has fixed boundary conditions
	 * @details
	 * When true, indicates this voxel maintains constant substrate concentrations
	 * (Dirichlet boundary condition).
	 * In cancer modeling, Dirichlet voxels often represent blood vessels or
	 * tissue boundaries with fixed nutrient/oxygen concentrations, which are
	 * critical determinants of tumor growth patterns and therapeutic response.
	 */
	bool es_dirichlet;
	
	/**
	 * @brief Outputs voxel data with specified units
	 * @param os Output stream
	 * @param unidades Units string (e.g., "micron")
	 * @details
	 * Formats voxel data with specified units for structured output.
	 * In cancer research, consistent units are crucial for meaningful
	 * interpretation of spatial data across different analysis tools.
	 */
	void stream_output_con_unidades( std::ostream& os , std::string unidades ) const;
};


#endif
