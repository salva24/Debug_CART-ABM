/*! @file Fase_Link.cpp
 *
 * @brief Implementation of phase transition mechanisms for cancer cell cycle modeling.
 *
 * @author Luciana Melina Luque
 *
 * @details
 * This file implements the Fase_Link class which represents transitions between 
 * cell cycle phases. The implementation provides the foundation for modeling 
 * dysregulated cell cycle progression in cancer, including initialization of transition 
 * parameters that may be altered in cancerous cells compared to normal cells.
 *
 * In cancer research, transitions between cell cycle phases are frequently altered, 
 * leading to uncontrolled proliferation. This implementation allows researchers to 
 * simulate these alterations and study their effects on tumor growth dynamics.
 *
 *
 * Inbound Dependencies:
 *   - Fase_Link.h - Class declaration
 *
 * Outbound Dependencies:
 *   - None
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */


#include "Fase_Link.h"

/*! \brief Constructor for Fase_Link
 *
 * Initializes a new phase transition with default settings.
 * In cancer modeling, these default settings provide a baseline
 * that can be modified to represent specific cancer types with
 * different cell cycle characteristics.
 *
 * The constructor initializes:
 * - Phase indices to 0 (unassigned)
 * - Fixed duration to false (stochastic transitions)
 * - Arrest function to NULL (no checkpoint)
 *
 * These parameters allow for flexible configuration of various cancer cell
 * cycle models, such as those with checkpoint deficiencies or altered
 * transition timing that contribute to carcinogenesis.
 */
Fase_Link::Fase_Link(){
	
	indice_fase_inicial = 0;
	indice_fase_final = 0; 
	duracion_fija = false;
	
	funcion_arrest = NULL;
	
	return;
}
