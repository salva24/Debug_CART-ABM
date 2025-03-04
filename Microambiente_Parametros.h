/**
 * @file Microambiente_Parametros.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Microenvironment parameter definitions for tumor microenvironment simulation
 *
 * @details Defines the parameters that configure the biochemical microenvironment in which
 * cancer cells exist. These parameters control spatial discretization, biochemical substrate
 * distributions, boundary conditions, and gradients - all critical aspects of modeling the
 * heterogeneous tumor microenvironment. The class enables realistic simulation of spatial
 * variations in oxygen, nutrients, and signaling molecules that profoundly affect cancer
 * cell behavior.
 * 
 * Inbound Dependencies: <vector>, <string>
 *
 * Outbound Dependencies: Used by Microambiente.h for environment configuration
 * 
 * Usage: Instantiated to define microenvironment properties in cancer simulations
 *
 * @license: Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef _MICROAMBIENTE_PARAMETROS_H
#define _MICROAMBIENTE_PARAMETROS_H


#include <vector>
#include <string>

/**
 * @brief Parameter class defining the biochemical tumor microenvironment
 * @details Encapsulates the parameters that define how the biochemical environment
 * surrounding cancer cells is simulated. This includes:
 * 1. Spatial discretization properties (dx, dy, dz)
 * 2. Boundary conditions that model interactions with surrounding tissues
 * 3. Initial and Dirichlet conditions for substrate concentrations
 * 4. Gradient calculation options
 * 
 * Cancer research applications include:
 * - Modeling hypoxic gradients that drive tumor progression and treatment resistance
 * - Simulating nutrient limitations that affect cancer cell metabolism
 * - Representing spatial heterogeneity in the tumor microenvironment
 * - Creating realistic boundary conditions that model tumor-tissue interfaces
 */
class Microambiente_Parametros
{
 private:
 
 public: 
	/**
	 * @brief Name identifier for this microenvironment
	 * @details Provides a descriptive label for the microenvironment, useful for
	 * tracking different microenvironment configurations in complex simulations.
	 */
	std::string nombre; 
 
	/**
	 * @brief Spatial units for the microenvironment (typically "micron")
	 * @details Defines the physical interpretation of spatial dimensions, important
	 * for matching simulation scales to real tumor measurements.
	 */
	std::string unidades_espaciales; 
	
	/**
	 * @brief Temporal units for the microenvironment (typically "min")
	 * @details Defines the physical interpretation of time, important for matching
	 * simulation kinetics to experimental measurements of cancer cell behavior.
	 */
	std::string unidades_temporales;
	
	/**
	 * @brief Spatial discretization in the x-dimension (microns)
	 * @details Controls the voxel size in the x-direction, affecting spatial resolution
	 * of biochemical gradients. Smaller values provide finer resolution but increase
	 * computational cost. Typically set to match the scale of individual cells.
	 */
	double dx;
	
	/**
	 * @brief Spatial discretization in the y-dimension (microns)
	 * @details Controls the voxel size in the y-direction, affecting spatial resolution
	 * of biochemical gradients across the tumor.
	 */
	double dy; 
	
	/**
	 * @brief Spatial discretization in the z-dimension (microns)
	 * @details Controls the voxel size in the z-direction, affecting spatial resolution
	 * of biochemical gradients in 3D tumor models.
	 */
	double dz; 
	
	/**
	 * @brief Flag for applying external Dirichlet conditions
	 * @details When true, enables fixed-value boundary conditions that represent
	 * constant sources/sinks, such as blood vessels supplying oxygen. Critical for
	 * modeling the tumor-vasculature interface in cancer research.
	 */
	bool condiciones_de_Dirichlet_externas; 
	
	/**
	 * @brief Vector of Dirichlet condition values for each substrate
	 * @details Fixed concentration values used at Dirichlet boundaries, typically
	 * representing physiological levels in well-perfused tissues. For oxygen, often
	 * set to 38 mmHg (5% oxygen) to model normal tissue oxygenation.
	 */
	std::vector<double> vector_condicion_de_dirichlet; 
	
	/**
	 * @brief Vector of flags controlling activation of Dirichlet conditions
	 * @details Enables or disables Dirichlet conditions for each substrate, allowing
	 * selective application of boundary conditions. Useful for modeling different
	 * tissue interfaces around tumors.
	 */
	std::vector<bool> vector_activacion_dirichlet;

	/**
	 * @brief Vector of flags for applying Dirichlet conditions to all boundaries
	 * @details When true for a substrate, applies the corresponding Dirichlet
	 * condition to all domain boundaries. Used to model tumors surrounded by
	 * well-vascularized tissues on all sides.
	 */
	std::vector<bool> dirichlet_todo; 
	
	/**
	 * @brief Vector of flags for Dirichlet conditions at minimum x boundary
	 * @details Controls whether fixed-value conditions are applied at the
	 * minimum x boundary for each substrate. Important for modeling directional
	 * gradients across tumors.
	 */
	std::vector<bool> dirichlet_xmin; 
	
	/**
	 * @brief Vector of flags for Dirichlet conditions at maximum x boundary
	 * @details Controls whether fixed-value conditions are applied at the
	 * maximum x boundary for each substrate.
	 */
	std::vector<bool> dirichlet_xmax; 
	
	/**
	 * @brief Vector of flags for Dirichlet conditions at minimum y boundary
	 * @details Controls whether fixed-value conditions are applied at the
	 * minimum y boundary for each substrate.
	 */
	std::vector<bool> dirichlet_ymin; 
	
	/**
	 * @brief Vector of flags for Dirichlet conditions at maximum y boundary
	 * @details Controls whether fixed-value conditions are applied at the
	 * maximum y boundary for each substrate.
	 */
	std::vector<bool> dirichlet_ymax; 
	
	/**
	 * @brief Vector of flags for Dirichlet conditions at minimum z boundary
	 * @details Controls whether fixed-value conditions are applied at the
	 * minimum z boundary for each substrate.
	 */
	std::vector<bool> dirichlet_zmin; 
	
	/**
	 * @brief Vector of flags for Dirichlet conditions at maximum z boundary
	 * @details Controls whether fixed-value conditions are applied at the
	 * maximum z boundary for each substrate.
	 */
	std::vector<bool> dirichlet_zmax; 

	/**
	 * @brief Vector of Dirichlet values for minimum x boundary
	 * @details Fixed concentration values at the minimum x boundary. Can be used
	 * to model concentration gradients across the tumor, such as oxygen declining
	 * from well-vascularized to poorly-perfused regions.
	 */
	std::vector<double> dirichlet_xmin_valores; 
	
	/**
	 * @brief Vector of Dirichlet values for maximum x boundary
	 * @details Fixed concentration values at the maximum x boundary.
	 */
	std::vector<double> dirichlet_xmax_valores; 
	
	/**
	 * @brief Vector of Dirichlet values for minimum y boundary
	 * @details Fixed concentration values at the minimum y boundary.
	 */
	std::vector<double> dirichlet_ymin_valores; 
	
	/**
	 * @brief Vector of Dirichlet values for maximum y boundary
	 * @details Fixed concentration values at the maximum y boundary.
	 */
	std::vector<double> dirichlet_ymax_valores; 
	
	/**
	 * @brief Vector of Dirichlet values for minimum z boundary
	 * @details Fixed concentration values at the minimum z boundary.
	 */
	std::vector<double> dirichlet_zmin_valores; 
	
	/**
	 * @brief Vector of Dirichlet values for maximum z boundary
	 * @details Fixed concentration values at the maximum z boundary.
	 */
	std::vector<double> dirichlet_zmax_valores; 

	/**
	 * @brief Vector of initial condition values for each substrate
	 * @details Starting concentration values throughout the domain at simulation
	 * initialization. Important for modeling initial tumor conditions before
	 * the cells begin to alter their environment.
	 */
	std::vector<double> vector_condiciones_iniciales; 
	
	/**
	 * @brief Vector defining x-dimension domain bounds (microns)
	 * @details Two-element vector [min, max] defining the x-extent of the simulation
	 * domain. Typically set to cover the entire tumor with sufficient margin.
	 */
	std::vector<double> rango_en_X; 
	
	/**
	 * @brief Vector defining y-dimension domain bounds (microns)
	 * @details Two-element vector [min, max] defining the y-extent of the simulation domain.
	 */
	std::vector<double> rango_en_Y; 
	
	/**
	 * @brief Vector defining z-dimension domain bounds (microns)
	 * @details Two-element vector [min, max] defining the z-extent of the simulation domain.
	 */
	std::vector<double> rango_en_Z; 
	
	/**
	 * @brief Flag for gradient calculation
	 * @details When true, enables calculation of concentration gradients across
	 * the microenvironment. Gradients are critical for modeling chemotaxis and
	 * directed cell migration in cancer, especially in invasion and metastasis.
	 */
	bool calcular_gradientes; 
	
	/**
	 * @brief Flag for using oxygen as the first substrate
	 * @details When true, ensures oxygen is the first substrate in the microenvironment,
	 * which aligns with the oxygen-specific parameters in other components. Critical
	 * for cancer modeling as hypoxia is a key driver of tumor progression.
	 */
	bool usar_oxigeno_como_primer_sustrato;
	
};

#endif
