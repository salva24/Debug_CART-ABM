/*! @file Fase.h
 *
 * @brief Cell phase definitions for cancer cell cycle modeling.
 *
 * @details
 * This file defines the Fase (Phase) class which represents distinct phases in 
 * a cell's life cycle. In cancer research, understanding cell cycle phases is 
 * critical for modeling tumor growth dynamics, as different cancer types show 
 * characteristic patterns of cell cycle regulation and dysregulation. This class
 * forms a fundamental building block for the cell cycle models used in the simulation.
 *
 * @author Luciana Melina Luque
 *
 * Inbound Dependencies:
 *   - Geometria.h - For volume calculations
 *   - Muerte_Parametros.h - Parameters for cell death processes
 *
 * Outbound Dependencies:
 *   - Ciclo_Modelo.h - Uses this class to build cycle models
 *   - Fase_Link.h - Uses this class for phase transitions
 *
 * @usage Used as a building block for cycle models in Ciclo_Modelo class
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __FASE_H__
#define __FASE_H__

#include <string>
#include "Geometria.h"
#include "Muerte_Parametros.h"

/*! \brief Class representing a cell cycle phase
 *
 * The Fase class models a distinct phase in a cell's life cycle.
 * In cancer research and simulation, cell phases are crucial for 
 * understanding tumor growth, as dysregulation of phase transitions
 * is a hallmark of cancer. This class enables modeling various aspects
 * of cell behavior during specific phases, including division, volume changes,
 * and death processes.
 */
class Fase{
	public:
		/*! \brief Index of the phase in the cell cycle model */
		int indice;
		
		/*! \brief Code identifying the specific phase type
		 *
		 * In cancer research, specific phase codes allow tracking cell populations
		 * in different cycle phases, important for analyzing tumor composition.
		 */
		int codigo;
		
		/*! \brief Name of the phase (e.g., "G0", "G1", "S", "G2", "M", "Apoptotic") */
		std::string nombre;
		
		/*! \brief Flag indicating if cell division occurs at the end of this phase
		 *
		 * This is crucial for modeling tumor growth, as it determines
		 * when new cancer cells are generated in the simulation.
		 */
		bool division_al_final_de_la_fase;
		
		/*! \brief Flag indicating if the cell should be removed at the end of this phase
		 *
		 * Used for modeling cell death processes (e.g., apoptosis, necrosis),
		 * which are important aspects of tumor dynamics and response to therapy.
		 */
		bool remover_al_final_de_la_fase;
		
		/*! \brief Flag indicating if cell volume should be updated during this phase
		 *
		 * Volume changes are important in cancer modeling, as different phases
		 * involve cell growth or shrinkage depending on the phase type.
		 */
		bool actualizar_volumen;
		
		/*! \brief Function pointer to phase entry function
		 *
		 * This function is called when a cell enters this phase.
		 * It can modify cell volume and other parameters based on
		 * the specific phase behavior required.
		 *
		 * \param volumen Reference to cell volume parameters to be modified
		 * \param mp Reference to death parameters that may be used
		 */
		void (*funcion_de_entrada)(Volumen& volumen, Muerte_parametros& mp);
		
		/*! \brief Default constructor
		 *
		 * Initializes a new phase with default values.
		 * All cancer cell cycle models start by constructing phases
		 * and then configuring them accordingly.
		 */
		Fase();

};

#endif
