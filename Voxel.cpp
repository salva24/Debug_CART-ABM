/**
 * @file Voxel.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of the Voxel class for spatial discretization in cancer modeling
 *
 * @details
 * This file implements the functionality of the Voxel class, which represents the fundamental
 * spatial unit in the discretized tumor microenvironment. Voxels are essential for modeling
 * the heterogeneous conditions within the tumor microenvironment, allowing for realistic
 * representation of spatially variable phenomena that influence cancer progression.
 * 
 * The implementation includes:
 * - Constructor for initial voxel setup (defaults optimized for cancer tissue simulation)
 * - Stream output operators for data visualization and analysis
 * - Formatted output with specified units for interdisciplinary research compatibility
 * 
 * In cancer research, this spatial discretization approach enables:
 * - Representation of localized hypoxic regions that drive tumor evolution
 * - Modeling of treatment penetration barriers in solid tumors
 * - Analysis of spatially distinct cellular populations within heterogeneous tumors
 * 
 * Inbound Dependencies: Voxel.h
 *
 * Outbound Dependencies: None
 *
 * @licence: Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#include "Voxel.h"

/**
 * @brief Default constructor for Voxel objects
 * @details
 * Initializes a voxel with standard default values suitable for cancer simulations:
 * - Index 0 (will be assigned by the microenvironment system)
 * - Volume of 1000 cubic micrometers (10x10x10), a common resolution for single-cell modeling
 * - Center position at origin (to be configured based on spatial location)
 * - Not a Dirichlet boundary node by default
 * 
 * In cancer modeling, these default values provide a reasonable starting point for
 * configuring the spatial domain, with dimensions appropriate for representing
 * cellular-scale phenomena while balancing computational efficiency.
 */
Voxel::Voxel()
{
	indice = 0; 
	volumen = 10*10*10;
//	centro; 
	es_dirichlet = false;
}

/**
 * @brief Stream output operator for voxel visualization
 * @param os Output stream to write to
 * @param v Voxel to output
 * @return Reference to the output stream
 * @details
 * Formats voxel data for structured output, including ID, position, and volume.
 * In cancer research, visualization of spatial data is critical for analyzing
 * tumor architecture, identifying regions of interest (such as hypoxic cores or
 * invasive margins), and correlating cellular behaviors with local microenvironmental
 * conditions.
 */
std::ostream& operator<<(std::ostream& os, const Voxel& v)  
{
	static std::string tabbing = "\t\t\t"; 
	static std::string tabbing2 = "\t\t\t\t"; 
	os	<< tabbing << "<voxel ID=\"" << v.indice << "\">"  << std::endl
		<< tabbing2 << "centro: " << v.centro << std::endl  
		<< tabbing2 << "volumen: " << v.volumen << std::endl; 

 return os; 
}

/**
 * @brief Outputs voxel data with specified measurement units
 * @param os Output stream to write to
 * @param unidades Units string (e.g., "micron")
 * @details
 * Generates structured XML-like output of voxel data with explicit units for each measurement.
 * In cancer research, consistent units are essential when sharing data across different
 * analysis tools or integrating simulation results with experimental measurements.
 * 
 * This method enables:
 * - Clear documentation of spatial scales in multiscale cancer models
 * - Compatibility with data analysis pipelines
 * - Integration with visualization tools for tumor architecture analysis
 * - Standardized formats for reporting spatial aspects of cancer simulations
 */
void Voxel::stream_output_con_unidades( std::ostream& os , std::string unidades ) const 
{
	static std::string tabbing = "\t\t\t\t"; 
	static std::string tabbing2 = "\t\t\t\t\t"; 
	os	<< tabbing << "<voxel ID=\"" << indice << "\">"  << std::endl
		<< tabbing2 << "<centro " << centro << " unidades=\"" << unidades << "\" />" << std::endl 
		<< tabbing2 << "<volumen= " << volumen << " unidades= " << unidades << " cubicos\">" << std::endl		
//		<< tabbing2 << "<unidades de volumen=\ " << unidades << " cubicos\">" << volumen << "</volumen>" << std::endl
		<< tabbing  << "</voxel>"; 
	return; 
}
