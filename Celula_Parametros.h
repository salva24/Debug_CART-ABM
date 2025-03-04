/**
 * @file Celula_Parametros.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Cell parameter definitions for cancer cell simulation
 *
 * @details Defines the parameters that govern individual cancer cell behaviors including
 * cell cycle progression, metabolic activities (secretion and consumption), and responses
 * to microenvironmental conditions. This class enables the customization of cell properties
 * to model different cancer cell types and behaviors.
 * 
 * Inbound Dependencies: Ciclo_Modelo.h, <vector>, <string>
 *
 * Outbound Dependencies: Used by Celula.h for cell parameterization
 * 
 * Usage: Instantiated to define cell-type specific behaviors in cancer simulations
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef _CELULA_PARAMETROS_H
#define _CELULA_PARAMETROS_H

#include "Ciclo_Modelo.h"
#include <vector>
#include <string>

/**
 * @brief Parameter class defining cell-type specific behaviors for cancer simulation
 * @details Encapsulates the key parameters that define how different types of cancer cells
 * behave in the simulation. This includes:
 * 1. Cell cycle behavior (proliferation rates, checkpoints, phases)
 * 2. Metabolic activities (secretion and consumption of substrates)
 * 3. Microenvironmental responses (responses to oxygen levels)
 * 4. Type identification (for heterogeneous tumor modeling)
 * 
 * Cancer research applications include:
 * - Modeling heterogeneous tumor populations with different cell types
 * - Simulating varying metabolic profiles between cancer cell subpopulations
 * - Representing different sensitivities to hypoxia across cell types
 * - Simulating cancer stem cell vs. differentiated cell behaviors
 */
class Celula_Parametros{
	
	public:
		/**
		 * @brief Cell cycle model for this cell type
		 * @details Defines the phases, transition rates, and division behavior
		 * for this specific cell type. Different cancer cell types can have
		 * significantly different cell cycle kinetics, from rapid cycling stem-like
		 * cells to slower cycling differentiated cells.
		 */
		Ciclo_Modelo ciclo;
		
		/**
		 * @brief Secretion rates of biochemical factors
		 * @details Defines how quickly this cell type secretes various factors
		 * into the microenvironment. In cancer, secreted factors can include
		 * growth factors, cytokines, and metabolic byproducts that affect 
		 * neighboring cells and promote tumor progression.
		 */
		double tasas_de_secrecion; 
		
		/**
		 * @brief Consumption rates of microenvironmental substrates
		 * @details Defines how quickly this cell type consumes various substrates
		 * like oxygen, glucose, etc. Cancer cells often exhibit the Warburg effect,
		 * with altered metabolic consumption patterns compared to normal cells.
		 */
		double tasas_de_consumo; 
		
		/**
		 * @brief Saturation densities for uptake/secretion
		 * @details Defines the substrate concentrations at which uptake or secretion
		 * reaches maximum rates. Different cancer cell metabolic profiles can be
		 * modeled by varying these saturation thresholds.
		 */
		double densidades_de_saturacion;
		
		/**
		 * @brief Name identifier for this cell type
		 * @details String identifier for the cell type (e.g., "cancer stem cell", 
		 * "hypoxic tumor cell", "invasive cell"). Useful for tracking and analyzing
		 * heterogeneous tumor populations.
		 */
		std::string nombre;
		
		/**
		 * @brief Numeric type identifier
		 * @details Integer ID for the cell type, used for efficient classification
		 * and processing. Enables quick identification of cell types in large
		 * simulations with heterogeneous tumor populations.
		 */
		int tipo;
		
		/**
		 * @brief Oxygen level for maximum proliferation (mmHg)
		 * @details Cell-type specific oxygen threshold above which this cell type
		 * proliferates at its maximum rate. Different cancer cell types can have
		 * varying sensitivities to hypoxia, with some cells better adapted to
		 * low-oxygen environments.
		 */
		double o2_saturacion_para_la_proliferacion;
		
		/**
		 * @brief Reference oxygen level (mmHg)
		 * @details Baseline oxygen level against which responses are calibrated
		 * for this specific cell type. Different cancer cell types may evolve in
		 * different oxygen environments and thus have different reference points.
		 */
		double o2_referencia;
	
};

#endif
