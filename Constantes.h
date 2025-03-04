/**
 * @file Constantes.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Global constants for the cell-based cancer simulation model
 *
 * @details This file defines essential constants used throughout the simulation,
 *          including time steps, cycle models, and cell state identifiers.
 * 
 * Inbound dependencies <string>
 *
 * Outbound dependencies Ciclo_Modelo.h, Muerte.h
 * 
 * @usage Include this file when accessing global simulation constants
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef _CONSTANTES_H
#define _CONSTANTES_H

#include <string>

/**
 * @brief Global constants for cancer cell simulation
 * 
 * @details This class defines constants used throughout the simulation to standardize
 *          time steps, cycle models, phase identifiers, and other critical parameters.
 *          These constants are essential for maintaining consistency across different
 *          components of the cancer model.
 * 
 * @cancer_research Consistent parameterization is essential in cancer modeling to:
 *                 - Enable reproducible simulations of tumor growth
 *                 - Standardize cell cycle and death phase identification
 *                 - Provide appropriate time scales for different biological processes
 *                 - Support comparison between different simulation scenarios
 */
class Constantes
{
 public:
 	/**
 	 * @brief Time step for diffusion processes (hours)
 	 * @details Controls the temporal resolution of biochemical diffusion in the tumor microenvironment
 	 * @cancer_research Small time step required for numerical stability in modeling rapid
 	 *                 diffusion of critical factors like oxygen and growth signals in tumors
 	 */
 	static constexpr double dt_difusion = 0.01;
 	
 	/**
 	 * @brief Time step for mechanical processes (hours)
 	 * @details Controls the temporal resolution of cell movement and physical interactions
 	 * @cancer_research Intermediate time step balances computational efficiency with accuracy
 	 *                 in modeling cell-cell forces in densely packed tumor regions
 	 */
 	static constexpr double dt_mecanica = 0.1;
 	
 	/**
 	 * @brief Time step for cell cycle processes (hours)
 	 * @details Controls the temporal resolution of biological processes like cell division
 	 * @cancer_research Larger time step appropriate for slower biological processes like
 	 *                 cell cycle progression, which typically operates on the order of hours to days
 	 */
 	static constexpr double dt_ciclo = 6.0;

	/**
	 * @brief Ki67-based cell cycle model identifier
	 * @details Identifies the cell cycle model based on Ki67 expression, a proliferation marker
	 * @cancer_research Ki67 is a widely used proliferation marker in cancer diagnostics
	 *                 and provides a clinically relevant model for tumor growth
	 */
	static const int ciclo_Ki67= 0;
	
	/**
	 * @brief Basic live/dead cell cycle model identifier
	 * @details Identifies a simplified cell cycle model with basic living state
	 * @cancer_research Simple model useful for cancer simulations focusing on
	 *                 viability rather than detailed proliferation dynamics
	 */
	static const int ciclo_vida= 1;

	/**
	 * @brief Necrosis death cycle model identifier
	 * @details Identifies the necrotic death pathway in the simulation
	 * @cancer_research Necrosis is common in poorly vascularized tumor regions
	 *                 and contributes to inflammation and tumor progression
	 */
	static const int ciclo_de_muerte_necrosis = 100;
	
	/**
	 * @brief Apoptosis death cycle model identifier
	 * @details Identifies the apoptotic death pathway in the simulation
	 * @cancer_research Dysregulated apoptosis is a hallmark of cancer,
	 *                 with cancer cells often resistant to programmed cell death
	 */
	static const int ciclo_de_muerte_apoptosis = 101;

	/**
	 * @brief Identifier for Ki67-positive pre-mitotic phase
	 * @details Identifies cells in the proliferative phase before mitosis
	 * @cancer_research High proportion of cells in this phase indicates
	 *                 aggressive tumor growth and poor prognosis
	 */
	static const int Ki67_positiva_premitotica=0;
	
	/**
	 * @brief Identifier for Ki67-positive post-mitotic phase
	 * @details Identifies newly divided cells still expressing Ki67
	 * @cancer_research Distinguishing pre- and post-mitotic cells allows
	 *                 more accurate modeling of tumor growth dynamics
	 */
	static const int Ki67_positiva_postmitotica=1;
	
	/**
	 * @brief Identifier for Ki67-negative phase
	 * @details Identifies quiescent cells not expressing Ki67
	 * @cancer_research Quiescent cancer cells may be more resistant to
	 *                 treatments targeting actively dividing cells
	 */
	static const int Ki67_negativa=2;
	
	/**
	 * @brief Identifier for basic living cell phase
	 * @details Identifies cells in the basic live/dead model
	 * @cancer_research Simple viability marker for basic cancer models
	 */
	static const int viva=3;

	/**
	 * @brief Identifier for swollen necrotic phase
	 * @details Identifies necrotic cells in the early swelling stage
	 * @cancer_research Early necrotic cells swell due to membrane damage,
	 *                 contributing to mechanical stress in tumors
	 */
	static const int necrotica_hinchada = 100;
	
	/**
	 * @brief Identifier for lysed necrotic phase
	 * @details Identifies necrotic cells that have ruptured
	 * @cancer_research Lysed cells release damage-associated molecular patterns
	 *                 that promote inflammation in the tumor microenvironment
	 */
	static const int necrotica_lisada = 101;
	
	/**
	 * @brief Identifier for generic necrotic phase
	 * @details Identifies cells in necrotic death when specific stage isn't relevant
	 * @cancer_research General necrotic marker for models not requiring
	 *                 fine temporal resolution of necrotic progression
	 */
	static const int necrotica = 102;
	
	/**
	 * @brief Identifier for cellular debris phase
	 * @details Identifies cellular remains after complete cell death
	 * @cancer_research Debris affects nutrient availability and immune cell infiltration
	 *                 in solid tumors, particularly in necrotic cores
	 */
	static const int debris = 104;
	
	/**
	 * @brief Identifier for apoptotic phase
	 * @details Identifies cells undergoing programmed cell death
	 * @cancer_research Apoptotic cells are removed cleanly without inflammation,
	 *                 in contrast to necrotic cells
	 */
	static const int apoptotica = 105;


};


#endif
