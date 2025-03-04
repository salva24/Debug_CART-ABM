/*! @file Fase.cpp
 *
 * @brief Implementation of the Fase (Phase) class for cancer cell cycle modeling.
 *
 * @details
 * This file implements the Fase class which represents cell cycle phases in
 * cancer simulation. The implementation initializes phase properties with default 
 * values that can be configured for specific phase types. Cell cycle phases
 * are critical for cancer research as they allow modeling of proliferation dynamics,
 * quiescence, and transitions to death phases like apoptosis and necrosis.
 *
 * @author Luciana Melina Luque
 *
 * Inbound Dependencies:
 *   - Fase.h - Class declaration
 *
 * Outbound Dependencies:
 *   - None
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#include "Fase.h"

/*!
 * \brief Constructor for the Fase class
 *
 * Initializes a new phase with default values. In cancer modeling,
 * phases are typically constructed and then configured with specific
 * properties based on their role in the cell cycle or death process.
 *
 * The default initialization creates a "neutral" phase that does not
 * cause division, removal, or volume updates. These properties must be
 * set specifically for each phase type (e.g., mitotic phases for division,
 * apoptotic phases for removal).
 *
 * Default values:
 * - Index and code set to 0
 * - Name set to "Sin nombre" (No name)
 * - No entry function
 * - No division, removal, or volume updates
 */
Fase::Fase(){
	
	indice = 0;
	codigo = 0;
	nombre = "Sin nombre";
	funcion_de_entrada = NULL;
	
	division_al_final_de_la_fase = false;
	remover_al_final_de_la_fase = false;
	actualizar_volumen = false;
	
	return;	
	
}
