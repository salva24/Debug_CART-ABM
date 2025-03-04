/*! @file Fase_Link.h
 *
 * @brief Phase transition definitions for cancer cell cycle modeling.
 *
 * @author Luciana Melina Luque
 *
 * @details
 * This file defines the Fase_Link class which represents transitions between 
 * cell cycle phases. In cancer research, the dysregulation of phase transitions 
 * is a critical aspect of tumor growth dynamics. These transitions determine when 
 * cells progress through their life cycle, including division and death processes,
 * and are often targeted by cancer therapies.
 *
 *
 * Inbound Dependencies:
 *   - Volumen.h - For volume-based transition conditions
 *   - Muerte_Parametros.h - Parameters for death-related transitions
 *
 * Outbound Dependencies:
 *   - Ciclo_Modelo.h - Uses this class to build cycle models
 *   - Fase.h - Links connect phases in the cycle model
 *
 * @usage Used to define transitions between cell cycle phases in Ciclo_Modelo class
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __FASE_LINK_H__
#define __FASE_LINK_H__

#include "Volumen.h"
#include "Muerte_Parametros.h"


/*! \brief Class representing transitions between cell cycle phases
 *
 * The Fase_Link class models connections between distinct phases in a cell's life cycle.
 * In cancer research, these transitions are crucial for understanding tumor growth,
 * as they determine when cells divide, progress through the cell cycle, or undergo death.
 * Many cancer therapies specifically target these phase transitions (e.g., taxanes targeting
 * mitotic transitions, CDK4/6 inhibitors targeting G1/S transitions).
 */
class Fase_Link{
	public:
		/*! \brief Index of the starting phase for this transition
         *
         * Identifies which phase the cell is transitioning from.
         * In cancer research, identifying the source phase helps
         * understand which part of the cell cycle may be dysregulated.
         */
		int indice_fase_inicial;
		
		/*! \brief Index of the ending phase for this transition
         *
         * Identifies which phase the cell is transitioning to.
         * Different cancer types may show characteristic patterns
         * of phase transition probabilities and timing.
         */
		int indice_fase_final;
		
		/*! \brief Flag indicating if this transition has a fixed duration
         *
         * When true, the transition occurs after a fixed time period.
         * When false, the transition occurs stochastically based on a rate.
         * This is important for modeling both deterministic (e.g., M phase) 
         * and stochastic (e.g., G1/S transition) aspects of cancer cell cycles.
         */
		bool duracion_fija;
		
		/*! \brief Function pointer to arrest function
         *
         * This function determines if a cell should be prevented from transitioning
         * between phases. Arrest functions model cell cycle checkpoints, which are
         * often dysfunctional in cancer cells, allowing them to proceed despite
         * cellular damage or unfavorable conditions.
         *
         * \param volumen Reference to cell volume parameters that may affect the arrest
         * \param mp Reference to death parameters that may affect the arrest
         * \return true if the cell should be arrested (prevented from transitioning),
         *         false if the cell is allowed to transition
         */
		bool (*funcion_arrest)( Volumen& volumen, Muerte_parametros& mp );
		
		/*! \brief Default constructor
         *
         * Initializes a new phase transition with default values.
         * All cell cycle models must define the transitions between phases
         * to create a complete cell cycle model.
         */
		Fase_Link();
	
};


#endif
