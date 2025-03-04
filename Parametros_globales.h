/**
 * @file Parametros_globales.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Global configuration parameters for the cell-based cancer simulation model
 *
 * @details
 * This file defines the central parameter repository that configures all aspects of the simulation,
 * including microenvironment properties, boundary conditions, cell properties, and immune response settings.
 * 
 * Inbound dependencies "Ciclo_Modelo.h", <vector>, <string>
 *
 * Outbound dependencies Various model components that access these global parameters
 * 
 * @usage Include this file to access or modify global simulation parameters
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __PARAMETROS_GLOBALES_H__
#define __PARAMETROS_GLOBALES_H__

#include "Ciclo_Modelo.h"
#include <vector>
#include <string>

/**
 * @brief Global configuration parameters for cancer cell simulation
 * 
 * @details This class serves as a central repository for all configuration parameters
 *          that govern the behavior of the simulation. It includes settings for the
 *          microenvironment, boundary conditions, cell properties, and immune response.
 *          The global instance 'pg' provides access to these parameters throughout the codebase.
 * 
 * @cancer_research Comprehensive parameter configuration is critical in cancer modeling to:
 *                 - Simulate diverse tumor microenvironments representing different cancer types
 *                 - Configure boundary conditions that model tissue-specific constraints
 *                 - Set cell properties that reflect cancer-specific behaviors
 *                 - Enable immunotherapy simulation through immune response parameters
 *                 - Control temporal and spatial resolution for different cancer processes
 */
class Parametros_globales{

	public:
    /**
     * @brief Random number generator seed
     * @details Controls reproducibility of stochastic processes in the simulation
     * @cancer_research Reproducible simulations are essential for methodical investigation
     *                 of stochastic cancer processes like mutation and treatment response
     */
    int seed;

	/**
	 * @brief Name of the microenvironment
	 * @details Descriptive identifier for the simulation environment
	 * @cancer_research Allows identification of specific tumor microenvironment configurations
	 */
	std::string m_nombre;

	/**
	 * @brief Spatial units label
	 * @details Units used for spatial measurements (typically microns)
	 * @cancer_research Ensures consistent spatial scaling for comparing simulation
	 *                 results with experimental cancer data
	 */
	std::string unidades_espaciales;
	
	/**
	 * @brief Temporal units label
	 * @details Units used for time measurements (typically hours)
	 * @cancer_research Ensures consistent temporal scaling for comparing simulation
	 *                 dynamics with experimental cancer progression data
	 */
	std::string unidades_temporales;
	
	/**
	 * @brief X-direction voxel size in microenvironment
	 * @details Spatial resolution in the X direction
	 * @cancer_research Controls the resolution of diffusion gradients that influence
	 *                 cancer cell behavior and migration
	 */
	double m_dx;
	
	/**
	 * @brief Y-direction voxel size in microenvironment
	 * @details Spatial resolution in the Y direction
	 * @cancer_research Controls the resolution of diffusion gradients that influence
	 *                 cancer cell behavior and migration
	 */
	double m_dy;
	
	/**
	 * @brief Z-direction voxel size in microenvironment
	 * @details Spatial resolution in the Z direction
	 * @cancer_research Controls the resolution of diffusion gradients that influence
	 *                 cancer cell behavior and migration
	 */
	double m_dz;

    /**
     * @brief Master periodicity flag
     * @details Controls whether domain boundaries repeat (true) or terminate (false)
     * @cancer_research Models boundary conditions that best represent the tumor environment,
     *                 from isolated tumors (non-periodic) to tissue samples (periodic)
     */
	bool condiciones_de_periodicidad;
	
	/**
	 * @brief X-direction periodicity flag
	 * @details Controls periodicity specifically in the X direction
	 * @cancer_research Directional periodicity can model tumors constrained by
	 *                 specific anatomical boundaries in one dimension
	 */
	bool condiciones_de_periodicidad_x;
	
	/**
	 * @brief Y-direction periodicity flag
	 * @details Controls periodicity specifically in the Y direction
	 * @cancer_research Directional periodicity can model tumors constrained by
	 *                 specific anatomical boundaries in one dimension
	 */
	bool condiciones_de_periodicidad_y;
	
	/**
	 * @brief Z-direction periodicity flag
	 * @details Controls periodicity specifically in the Z direction
	 * @cancer_research Directional periodicity can model tumors constrained by
	 *                 specific anatomical boundaries in one dimension
	 */
	bool condiciones_de_periodicidad_z;

	/**
	 * @brief External Dirichlet boundary conditions flag
	 * @details Controls whether fixed-value boundary conditions are applied at domain boundaries
	 * @cancer_research Models constant-concentration sources like blood vessels or
	 *                 tissue boundaries that influence tumor growth and metabolism
	 */
	bool condiciones_de_Dirichlet_externas;
	
	/**
	 * @brief Vector of Dirichlet condition values
	 * @details Concentration values for fixed-value boundaries for each substrate
	 * @cancer_research Defines oxygen, nutrient, or drug concentrations at blood vessels
	 *                 or tissue boundaries surrounding the tumor
	 */
	std::vector<double> vector_condicion_de_dirichlet;
	
	/**
	 * @brief Vector of Dirichlet condition activation flags
	 * @details Controls which substrates use Dirichlet boundary conditions
	 * @cancer_research Allows selective application of fixed boundaries for
	 *                 different factors, modeling varied tissue interfaces
	 */
	std::vector<bool> vector_activacion_dirichlet;

	/**
	 * @brief Vector of flags for all-boundary Dirichlet conditions
	 * @details Controls whether each substrate has Dirichlet conditions on all boundaries
	 * @cancer_research Models tumors completely surrounded by a tissue with fixed
	 *                 biochemical concentrations
	 */
	std::vector<bool> dirichlet_todo;

	/**
	 * @brief Vector of flags for Dirichlet conditions at minimum X boundary
	 * @details Controls which substrates have fixed values at the X minimum boundary
	 * @cancer_research Models directional substrate sources like blood vessels at
	 *                 specific tumor boundaries
	 */
	std::vector<bool> dirichlet_xmin;
	
	/**
	 * @brief Vector of flags for Dirichlet conditions at maximum X boundary
	 * @details Controls which substrates have fixed values at the X maximum boundary
	 * @cancer_research Models directional substrate sources like blood vessels at
	 *                 specific tumor boundaries
	 */
	std::vector<bool> dirichlet_xmax;
	
	/**
	 * @brief Vector of flags for Dirichlet conditions at minimum Y boundary
	 * @details Controls which substrates have fixed values at the Y minimum boundary
	 * @cancer_research Models directional substrate sources like blood vessels at
	 *                 specific tumor boundaries
	 */
	std::vector<bool> dirichlet_ymin;
	
	/**
	 * @brief Vector of flags for Dirichlet conditions at maximum Y boundary
	 * @details Controls which substrates have fixed values at the Y maximum boundary
	 * @cancer_research Models directional substrate sources like blood vessels at
	 *                 specific tumor boundaries
	 */
	std::vector<bool> dirichlet_ymax;
	
	/**
	 * @brief Vector of flags for Dirichlet conditions at minimum Z boundary
	 * @details Controls which substrates have fixed values at the Z minimum boundary
	 * @cancer_research Models directional substrate sources like blood vessels at
	 *                 specific tumor boundaries
	 */
	std::vector<bool> dirichlet_zmin;
	
	/**
	 * @brief Vector of flags for Dirichlet conditions at maximum Z boundary
	 * @details Controls which substrates have fixed values at the Z maximum boundary
	 * @cancer_research Models directional substrate sources like blood vessels at
	 *                 specific tumor boundaries
	 */
	std::vector<bool> dirichlet_zmax;
	
	/**
	 * @brief Vector of flags for Dirichlet conditions at blood vessel locations
	 * @details Controls which substrates have fixed values at blood vessel nodes
	 * @cancer_research Models substrate delivery through vasculature within tumors,
	 *                 essential for studying angiogenesis and drug delivery
	 */
	std::vector<bool> dirichlet_vs;

	/**
	 * @brief Values for Dirichlet conditions at minimum X boundary
	 * @details Fixed concentration values for each substrate at X minimum boundary
	 * @cancer_research Defines specific oxygen, nutrient, or drug concentrations
	 *                 at directional tissue boundaries
	 */
	std::vector<double> dirichlet_xmin_valores;
	
	/**
	 * @brief Values for Dirichlet conditions at maximum X boundary
	 * @details Fixed concentration values for each substrate at X maximum boundary
	 * @cancer_research Defines specific oxygen, nutrient, or drug concentrations
	 *                 at directional tissue boundaries
	 */
	std::vector<double> dirichlet_xmax_valores;
	
	/**
	 * @brief Values for Dirichlet conditions at minimum Y boundary
	 * @details Fixed concentration values for each substrate at Y minimum boundary
	 * @cancer_research Defines specific oxygen, nutrient, or drug concentrations
	 *                 at directional tissue boundaries
	 */
	std::vector<double> dirichlet_ymin_valores;
	
	/**
	 * @brief Values for Dirichlet conditions at maximum Y boundary
	 * @details Fixed concentration values for each substrate at Y maximum boundary
	 * @cancer_research Defines specific oxygen, nutrient, or drug concentrations
	 *                 at directional tissue boundaries
	 */
	std::vector<double> dirichlet_ymax_valores;
	
	/**
	 * @brief Values for Dirichlet conditions at minimum Z boundary
	 * @details Fixed concentration values for each substrate at Z minimum boundary
	 * @cancer_research Defines specific oxygen, nutrient, or drug concentrations
	 *                 at directional tissue boundaries
	 */
	std::vector<double> dirichlet_zmin_valores;
	
	/**
	 * @brief Values for Dirichlet conditions at maximum Z boundary
	 * @details Fixed concentration values for each substrate at Z maximum boundary
	 * @cancer_research Defines specific oxygen, nutrient, or drug concentrations
	 *                 at directional tissue boundaries
	 */
	std::vector<double> dirichlet_zmax_valores;

	/**
	 * @brief Initial condition values for substrates
	 * @details Starting concentration values for each substrate throughout the domain
	 * @cancer_research Sets baseline biochemical conditions that represent the
	 *                 pre-cancer or early tumor microenvironment
	 */
	std::vector<double> vector_condiciones_iniciales;

	/**
	 * @brief X-dimension range for microenvironment
	 * @details Minimum and maximum X coordinates of the simulation domain
	 * @cancer_research Defines the physical size of the tumor model in X direction
	 */
	std::vector<double> rango_en_X;
	
	/**
	 * @brief Y-dimension range for microenvironment
	 * @details Minimum and maximum Y coordinates of the simulation domain
	 * @cancer_research Defines the physical size of the tumor model in Y direction
	 */
	std::vector<double> rango_en_Y;
	
	/**
	 * @brief Z-dimension range for microenvironment
	 * @details Minimum and maximum Z coordinates of the simulation domain
	 * @cancer_research Defines the physical size of the tumor model in Z direction
	 */
	std::vector<double> rango_en_Z;

	/**
	 * @brief Gradient calculation flag
	 * @details Controls whether chemical gradients are calculated during simulation
	 * @cancer_research Gradients drive chemotaxis and influence cancer cell migration,
	 *                 crucial for modeling invasive behavior
	 */
	bool calcular_gradientes;

	/**
	 * @brief Oxygen as first substrate flag
	 * @details Controls whether oxygen is the first substrate in the microenvironment
	 * @cancer_research Standardizes oxygen position for hypoxia-dependent models,
	 *                 as oxygen is the most critical substrate in tumor growth models
	 */
	bool usar_oxigeno_como_primer_sustrato;


	/**
	 * @brief X-direction mechanical container voxel size
	 * @details Spatial resolution of the cell container in X direction
	 * @cancer_research Controls the granularity of mechanical interactions between
	 *                 cancer cells, affecting tumor density and morphology
	 */
	double c_dx;
	
	/**
	 * @brief Y-direction mechanical container voxel size
	 * @details Spatial resolution of the cell container in Y direction
	 * @cancer_research Controls the granularity of mechanical interactions between
	 *                 cancer cells, affecting tumor density and morphology
	 */
	double c_dy;
	
	/**
	 * @brief Z-direction mechanical container voxel size
	 * @details Spatial resolution of the cell container in Z direction
	 * @cancer_research Controls the granularity of mechanical interactions between
	 *                 cancer cells, affecting tumor density and morphology
	 */
	double c_dz;


	/**
	 * @brief Cell cycle model
	 * @details Model governing cell division and death decisions
	 * @cancer_research Defines proliferation characteristics of the specific
	 *                 cancer type being modeled
	 */
	Ciclo_Modelo ciclo;
	
	/**
	 * @brief Secretion rates
	 * @details Rate at which cells secrete biochemical factors
	 * @cancer_research Models how cancer cells modify their microenvironment
	 *                 through signaling molecules and metabolites
	 */
	double tasas_de_secrecion;
	
	/**
	 * @brief Consumption rates
	 * @details Rate at which cells consume biochemical factors
	 * @cancer_research Models cancer cells' metabolic demands, often elevated
	 *                 compared to normal cells (Warburg effect)
	 */
	double tasas_de_consumo;
	
	/**
	 * @brief Saturation densities
	 * @details Concentration thresholds for maximum secretion/consumption
	 * @cancer_research Models how cancer cells respond to varying levels of
	 *                 nutrients, oxygen, and other factors
	 */
	double densidades_de_saturacion;
	
	/**
	 * @brief Cell type name
	 * @details Descriptive identifier for the modeled cell type
	 * @cancer_research Identifies specific cancer cell types (e.g., MCF-7, A549)
	 *                 being simulated
	 */
	std::string c_nombre;
	
	/**
	 * @brief Cell type identifier
	 * @details Numeric identifier for the cell type
	 * @cancer_research Distinguishes between different cell populations in
	 *                 heterogeneous tumor models
	 */
	int tipo;
	
	/**
	 * @brief Oxygen saturation threshold for proliferation
	 * @details Oxygen concentration above which proliferation occurs normally
	 * @cancer_research Models how hypoxia affects cancer cell division,
	 *                 a critical factor in tumor growth kinetics
	 */
	double o2_saturacion_para_la_proliferacion;
	
	/**
	 * @brief Reference oxygen concentration
	 * @details Baseline oxygen level for cell behavior calculations
	 * @cancer_research Defines normoxic conditions for comparing hypoxic effects,
	 *                 important for modeling tumor regions at different depths
	 */
	double o2_referencia;
	
	/**
	 * @brief Mechanical boundary interaction flag
	 * @details Controls whether cells interact with domain boundaries
	 * @cancer_research Models how tumors interact with physical barriers in tissues,
	 *                 such as basement membranes or organ boundaries
	 */
	bool interactuar_con_mb;
	
	/**
	 * @brief Side growth flag
	 * @details Controls whether cells can grow adjacent to barriers
	 * @cancer_research Models how cancer cells adapt to physical constraints,
	 *                 relevant for invasion through tissue boundaries
	 */
	bool crecer_al_costado;
	
	/**
	 * @brief Cell identifier number
	 * @details Unique identifier for the cell
	 * @cancer_research Enables tracking of individual cancer cell lineages
	 *                 through time to study clonal evolution
	 */
	int numero_id;

	/**
	 * @brief Set parameters from configuration files
	 * @param xmlname Name of the XML configuration file
	 * @param cellname Name of cell type to configure
	 * @details Loads simulation parameters from specified configuration files
	 * @cancer_research Enables standardized configuration of different cancer
	 *                 simulation scenarios for comparative studies
	 */
	void set_parametros(std::string, std::string);

	/**
	 * @brief Mean value for immune cell properties
	 * @details Average value for stochastic immune cell parameter generation
	 * @cancer_research Models the average effectiveness of tumor-infiltrating
	 *                 lymphocytes for immunotherapy simulations
	 */
	double imm_mean;
	
	/**
	 * @brief Standard deviation for immune cell properties
	 * @details Variance in stochastic immune cell parameter generation
	 * @cancer_research Models heterogeneity in immune cell populations,
	 *                 critical for realistic immunotherapy response prediction
	 */
	double imm_sd;
	
	/**
	 * @brief Immune response activation flag
	 * @details Controls whether immune response is simulated
	 * @cancer_research Enables simulation of immunotherapy or natural
	 *                 immune responses against tumors
	 */
	bool activar_respuesta_inmune;
    
    /**
     * @brief Number of lymphocytes to simulate
     * @details Controls immune cell population size
     * @cancer_research Models immune infiltration density, a key
     *                 prognostic factor in many cancer types
     */
    int cantidad_de_linfocitos;
    
    /**
     * @brief Primary immune response timing
     * @details Time point for initial immune response
     * @cancer_research Models the delay before immune recognition
     *                 of developing tumors
     */
    double tiempo_de_imm;
    
    /**
     * @brief Secondary immune response timing
     * @details Time point for secondary immune response
     * @cancer_research Models sequential immunotherapy treatments
     *                 or adaptive immune response stages
     */
    double tiempo_de_imm_2;

    /**
     * @brief Total simulation time counter
     * @details Tracks the elapsed simulation time
     * @cancer_research Monitors progression through the simulated
     *                 timeframe of cancer development
     */
    double tiempo_total;
    
    /**
     * @brief End time for simulation
     * @details Maximum time point for the simulation
     * @cancer_research Defines the temporal extent of the cancer
     *                 progression being modeled
     */
	double tiempo_final;

};

/**
 * @brief Global parameter instance
 * @details Externally accessible parameter object for simulation configuration
 * @cancer_research Provides centralized access to all cancer model parameters
 *                 throughout the simulation
 */
extern Parametros_globales *pg;

#endif
