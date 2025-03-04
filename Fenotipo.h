/**
 * @file Fenotipo.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Cell phenotype definition integrating all behavioral components for cancer research
 *
 * @details This file defines the Fenotipo (Phenotype) class, which coordinates all aspects of
 * cell behavior in the cancer simulation. The phenotype represents the observable characteristics
 * of a cell, integrating cycle progression, physical properties, secretion patterns, and death
 * processes into a cohesive functional unit.
 * 
 * In cancer research, the phenotype concept is critical as it represents the manifestation of
 * genetic alterations and microenvironmental influences that drive cancer progression. The
 * modular design allows representation of the diverse cellular phenotypes found in heterogeneous
 * tumors and modeling of phenotypic plasticity in response to environmental stressors.
 * 
 * Inbound dependencies Ciclo.h, Mecanica.h, Secrecion.h, Muerte.h
 *
 * Outbound dependencies Celula.h
 *
 * @usage Instantiated as a member variable within each Celula object
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __FENOTIPO_H__
#define __FENOTIPO_H__

#include "Ciclo.h"
#include "Mecanica.h"
#include "Secrecion.h"
#include "Muerte.h"

/**
 * @class Fenotipo
 * @brief Integrated cell phenotype model that combines all functional aspects of a cancer cell
 * 
 * @details The Fenotipo class serves as the central coordinator for all aspects of cell behavior,
 * integrating cell cycle progression, volume regulation, geometric properties, mechanical interactions,
 * secretion/consumption patterns, and death pathways into a unified representation of cellular phenotype.
 * 
 * In cancer research, this integration is essential as it allows modeling of:
 * - Phenotypic heterogeneity within tumors (different cancer cell subpopulations)
 * - Phenotypic plasticity in response to microenvironmental changes
 * - Emergent behaviors arising from interactions between phenotypic processes
 * - Treatment effects that may target multiple aspects of the cancer cell phenotype
 * - Transitions between phenotypic states during cancer progression
 * 
 * The modular design enables representation of diverse cancer phenotypes from highly
 * proliferative to invasive to drug-resistant variants, supporting investigation of
 * the complex phenotypic landscape of cancer.
 */
class Fenotipo
{
 private:
 public:
	 
	/**
	 * @brief Cell cycle model controlling proliferation behavior
	 * @details In cancer cells, the cycle is often dysregulated with alterations in
	 * checkpoint controls and progression rates, leading to uncontrolled proliferation
	 */
	Ciclo ciclo; 
	
	/**
	 * @brief Cell volume properties and regulation
	 * @details Cancer cells often exhibit altered volume regulation with changes in
	 * nuclear-to-cytoplasmic ratio, fluid content, and size variability, all of which
	 * affect drug distribution and cellular stress responses
	 */
	Volumen volumen; 
	
	/**
	 * @brief Geometric properties of the cell
	 * @details Cell geometry is altered in cancer, with changes in nuclear size, surface area,
	 * and polarity that affect diagnostic appearance, mechanical properties, and invasive capacity
	 */
	Geometria geometria;
	
	/**
	 * @brief Mechanical properties governing cell-cell and cell-ECM interactions
	 * @details Cancer cells exhibit altered adhesion, stiffness, and motility properties
	 * that enable invasion, migration through tissue, and metastatic spread
	 */
	Mecanica mecanica;
	
	/**
	 * @brief Secretion and consumption patterns for interaction with microenvironment
	 * @details Cancer cells show modified metabolic profiles (Warburg effect) and secrete
	 * factors that remodel the microenvironment to promote growth, angiogenesis, and immune evasion
	 */
	Secrecion secrecion;
	
	/**
	 * @brief Cell death pathways and properties
	 * @details Cancer cells frequently develop resistance to apoptosis and may undergo
	 * alternative death pathways like necrosis, affecting tumor architecture and immune response
	 */
	Muerte muerte;

	/**
	 * @brief Default constructor for the Fenotipo class
	 * @details Creates a phenotype with default values for all components,
	 * initialized to represent a standard cancer cell behavior profile.
	 * In cancer research, this initialization provides a baseline phenotype
	 * that can be modified to represent different cancer cell types or states.
	 */
	Fenotipo(); // done 
	
};

#endif
