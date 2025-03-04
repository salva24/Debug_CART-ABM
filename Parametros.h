/**
 * @file Parametros.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Oxygen parameter definitions for cancer cell simulation
 *
 * @details Defines oxygen thresholds that govern cancer cell behaviors, particularly
 * proliferation and necrosis. These parameters model how cancer cells respond to
 * varying oxygen levels in the tumor microenvironment.
 * 
 * Inbound Dependencies: <string>
 * Outbound Dependencies: Various cell behavior models
 * 
 * Usage: Instantiated by Celula and Microambiente classes to define oxygen-based responses
 *
 * @license: Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __PARAMETROS_H__
#define __PARAMETROS_H__

#include <string>

/**
 * @brief Parameter class defining oxygen thresholds and cell responses
 * @details Models how cancer cells respond to different oxygen concentrations in the tumor
 * microenvironment. Cancer cells exhibit altered metabolic responses to hypoxia compared
 * to normal cells, including the activation of HIF-1α pathways, metabolic reprogramming,
 * and resistance to hypoxia-induced death. These parameters capture key oxygen thresholds
 * that determine cell fate decisions.
 */
class Parametros{

	public:
		// oxygen values (in mmHg)
	//double o2_hypoxia_limite;
	//double o2_hypoxia_respuesta;
	//double o2_hypoxia_saturacion;

	/** 
	 * @brief Oxygen level above which proliferation rate is maximal (mmHg)
	 * @details Cancer cells can continue proliferating at lower oxygen levels than normal cells,
	 * but still require sufficient oxygen for maximal division rates. This parameter defines
	 * the oxygen threshold above which proliferation is not limited by oxygen availability.
	 */
	double o2_saturacion_para_la_proliferacion;
	
	/**
	 * @brief Minimum oxygen level required for proliferation (mmHg)
	 * @details Below this threshold, cancer cells cease proliferation but may remain viable.
	 * This models the phenomenon where tumor cells in hypoxic regions arrest in G0/G1 phase
	 * but can resume division when reoxygenated.
	 */
	double o2_limite_de_proliferacion;

	/**
	 * @brief Reference oxygen level for normal tissue (mmHg)
	 * @details Used as a baseline comparison for oxygen-dependent behaviors. Typically set to
	 * 38 mmHg, corresponding to 5% oxygen, which approximates normal tissue oxygenation.
	 * Cancer cells often evolve in environments with fluctuating oxygen levels, so this
	 * reference point is critical for modeling relative responses.
	 */
	double o2_referencia;

	/**
	 * @brief Oxygen level below which cells begin to die from necrosis (mmHg)
	 * @details When oxygen drops below this threshold, cancer cells begin to undergo necrotic death.
	 * This models tumor necrotic cores, which form in regions far from blood vessels where oxygen
	 * levels are severely depleted. Cancer cells are often more resistant to hypoxia than normal cells.
	 */
	double o2_necrosis_limite;
	
	/**
	 * @brief Oxygen level at which necrosis rate reaches maximum (mmHg)
	 * @details Below this extremely low oxygen threshold, cells undergo necrosis at the maximum rate.
	 * This parameter helps model the gradient of necrotic death observed in tumors, with increasing
	 * rates of cell death as distance from blood vessels increases.
	 */
	double o2_necrosis_max;

	/**
	 * @brief Maximum rate of necrotic cell death under severe hypoxia (1/min)
	 * @details Defines how quickly cells die when oxygen is depleted below o2_necrosis_max.
	 * Cancer cells often develop mechanisms to survive in hypoxic environments, but extreme
	 * hypoxia still leads to necrotic death. This parameter is typically set to correspond to
	 * approximately 6 hours survival time in severe hypoxia.
	 */
	double tasa_necrosis_max;

	/**
	 * @brief Default constructor initializing oxygen parameters
	 * @details Sets default values for oxygen thresholds based on experimental measurements
	 * of cancer cell responses to varying oxygen levels. These defaults model general cancer
	 * cell behavior, but can be adjusted to represent specific cancer types with different
	 * hypoxia tolerance profiles.
	 */
	Parametros();

};

#endif
