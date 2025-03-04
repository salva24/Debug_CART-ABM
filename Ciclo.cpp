/**
 * @file Ciclo.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of the Cycle class
 *
 * @details This file implements the Ciclo (Cycle) class methods.
 *
 * Inbound dependencies Ciclo.h
 *
 * Outbound_dependencies
 *
 * @usage Instantiated as a member variable within each Fenotipo object
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

/**
 * @file Ciclo.cpp
 * @author Luciana Melina Luque
 * @brief Implementation of individual cell cycle progression for cancer simulation
 * @details This file implements the Ciclo (Cycle) class, which tracks the state and
 * progression of an individual cell's cycle. While the Ciclo_Modelo class defines the
 * overall structure of possible cell cycles, this class manages a specific cell's
 * journey through that cycle model.
 * 
 * In cancer research, this implementation enables:
 * - Tracking of cell-specific cycle progression
 * - Modeling heterogeneous proliferation rates within tumors
 * - Simulating cell responses to microenvironmental conditions
 * - Representing cancer-specific cell cycle dysregulation
 * - Integration between cycle progression and death pathways
 * 
 * @inbound_dependencies Ciclo.h
 * @outbound_dependencies None
 * @license Private
 */

#include "Ciclo.h"
#include <cstddef>

/**
 * @brief Default constructor for the Ciclo class
 * @details Initializes a new cycle instance with NULL model pointer and default values
 * for all tracking variables. The cycle must be connected to a specific model using
 * sync_con_ciclo_modelo() before use.
 * 
 * In cancer modeling, this initialization creates the foundation for tracking an
 * individual cancer cell's progression through its life cycle, with appropriate
 * flags for monitoring division and death events.
 * 
 * @return void
 */
Ciclo::Ciclo(){

	pCiclo_Modelo = NULL;
	indice_de_la_fase_actual = 0;
	tiempo_acumulado_en_la_fase = 0;
	flagged_para_dividirse= false;
	flagged_para_remover= false;
    tasa_aleatoria = 0;

	return;
}

/**
 * @brief Synchronizes this cycle instance with a cycle model
 * @param cm Reference to the cycle model to be used
 * @details Links this cycle instance to a specific model and copies the transition rates
 * from the model to this instance. This connection allows the cycle to follow the rules
 * defined by the model while maintaining its own state.
 * 
 * In cancer research, this enables attaching different cycle models (like Ki67 cycle vs
 * basic live/dead cycle) to different cell populations, supporting heterogeneous
 * modeling of cancer cells with diverse proliferation characteristics.
 * 
 * @return void
 */
void Ciclo::sync_con_ciclo_modelo( Ciclo_Modelo& cm ){

	pCiclo_Modelo = &cm;
	tasas_de_transicion = cm.tasas_de_transicion;
	return;

}

//void Ciclo::avanzar_en_el_ciclo( int& indice_de_la_fase_actual, double& tiempo_acumulado_en_la_fase, Volumen& volumen, double dt ){
//
//	pCiclo_Modelo->avanzar_en_el_modelo(indice_de_la_fase_actual, tiempo_acumulado_en_la_fase, volumen, dt);
//	return;
//}

/**
 * @brief Advances the cell through its cycle for a time step
 * @param volumen Reference to the cell's volume
 * @param dt Time step size
 * @param c_tasas_de_transicion Current transition rates, potentially modified by conditions
 * @param mp Death parameters affecting cycle progression
 * @details Progresses the cell through its cycle for the given time step by delegating to
 * the linked cycle model. Updates the phase, accumulated time, and potentially sets
 * flags for division or removal.
 * 
 * In cancer research, this method simulates how cancer cells progress through altered
 * cell cycles, responding to various microenvironmental conditions that may affect
 * transition rates, growth patterns, and susceptibility to death pathways.
 * 
 * @return void
 */
void Ciclo::avanzar_en_el_ciclo(Volumen& volumen, double dt, std::vector< std::vector<double> >& c_tasas_de_transicion, Muerte_parametros& mp ){

	pCiclo_Modelo->avanzar_en_el_modelo(flagged_para_remover, flagged_para_dividirse, indice_de_la_fase_actual, tiempo_acumulado_en_la_fase, volumen, dt, c_tasas_de_transicion, mp);
	return;
}

/**
 * @brief Checks if volume should be updated in current phase
 * @return Boolean indicating if volume update is needed
 * @details Queries the cycle model to determine if the current phase requires
 * volume updates. Different phases may have different growth requirements.
 * 
 * In cancer modeling, reflects how cell growth varies across cycle phases, with
 * cancer cells often showing altered growth patterns in specific phases due to
 * oncogenic signaling or metabolic reprogramming.
 */
bool Ciclo::actualizar_volumen(){

	return pCiclo_Modelo->get_actualizar_volumen(indice_de_la_fase_actual);

}

/**
 * @brief Gets the current phase's transition rate
 * @return Transition rate from current phase
 * @details Retrieves the base transition rate for the current phase from the cycle model.
 * This rate determines how quickly cells move through this phase.
 * 
 * In cancer models, transition rates often reflect oncogenic alterations that
 * accelerate specific phase transitions, such as dysregulated G1/S transitions
 * due to Rb/E2F pathway mutations common in many cancers.
 */
double Ciclo::tasa_de_transicion(){

	double tasa;
	tasa = pCiclo_Modelo->get_tasa_de_transicion(indice_de_la_fase_actual);
	return tasa;
}

/**
 * @brief Updates transition rate between specified phases
 * @param fase_actual Current phase index
 * @param fase_siguiente Target phase index
 * @return Reference to the specific transition rate
 * @details Provides direct access to modify a specific transition rate between
 * the given phases. Uses the model's inverse mapping to find the correct index.
 * 
 * In cancer research, enables representing environmental influences on specific
 * transition rates, such as hypoxia-induced slowing of G1/S transition or
 * acidosis effects on mitotic entry, critical for modeling tumor microenvironment
 * effects on cancer growth dynamics.
 */
double& Ciclo::actualizar_mis_tasas_de_transicion(int fase_actual, int fase_siguiente){

    return tasas_de_transicion[fase_actual][pCiclo_Modelo->get_indice_de_mapa_inverso(fase_actual,fase_siguiente)];


}

/**
 * @brief Gets reference to the current phase object
 * @return Reference to current Fase object
 * @details Provides access to the detailed properties of the current cell cycle phase
 * from the linked cycle model. Allows examination and modification of phase-specific
 * properties.
 * 
 * In cancer research, allows access to phase-specific behaviors that may be
 * dysregulated in cancer cells or targeted by phase-specific therapeutics,
 * such as S-phase specific chemotherapies or mitotic inhibitors.
 */
Fase& Ciclo::fase_actual(){

	return pCiclo_Modelo->fases[indice_de_la_fase_actual];

}

