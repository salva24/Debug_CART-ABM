/**
 * @file Ciclo.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Cell cycle instance definition for individual cancer cells
 *
 * @details This file defines the Ciclo (Cycle) class, which represents the current state
 * of an individual cell's progression through its life cycle. While Ciclo_Modelo defines
 * the structure and rules of cell cycles, this class tracks the actual state of a specific
 * cell as it progresses through that model.
 * 
 * In cancer research, tracking individual cell cycle progression is critical for:
 * - Monitoring proliferation rates in different tumor regions
 * - Analyzing cycle checkpoint dysregulation in cancer cells
 * - Studying effects of therapeutic agents on cell cycle arrest
 * - Modeling cancer cell responses to microenvironmental conditions
 * - Simulating heterogeneous division rates within tumors
 * 
 * Inbound dependencies Ciclos_estandares.h
 *
 * Outbound dependencies Fenotipo.h
 *
 * @usage Instantiated as a member variable within each Fenotipo object
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __CICLO_H__
#define __CICLO_H__


#include "Ciclos_estandares.h"

/**
 * @class Ciclo
 * @brief Tracks and manages an individual cell's progression through its cell cycle
 * 
 * @details The Ciclo class represents the current state of an individual cell's
 * progression through the cell cycle. It maintains a pointer to a Ciclo_Modelo that
 * defines the structure of the cycle (phases, transitions), while tracking the cell's
 * current phase, time spent in that phase, and flags for important state changes like
 * division or death.
 * 
 * In cancer research, this class enables modeling of:
 * - Individual cell cycle dysregulation in cancer
 * - Heterogeneous cell cycle rates in tumors
 * - Cell cycle arrests in response to therapy
 * - Proliferation-quiescence decisions at the single-cell level
 * - Death pathway integration with cell cycle state
 */
class Ciclo{

	public:
	/**
	 * @brief Pointer to the cycle model defining structure and rules
	 * @details Links to the specific cycle model (e.g., Ki67 cycle, live/dead cycle)
	 * that defines the phases and transition rules. In cancer, different models 
	 * can represent various cancer types or states (e.g., highly proliferative vs.
	 * dormant cancer cells).
	 */
	Ciclo_Modelo* pCiclo_Modelo;
	
	/**
	 * @brief Current phase index within the cell cycle
	 * @details Tracks which phase (e.g., G0/G1, S, G2, M) the cell is currently in.
	 * In cancer research, phase distribution analysis reveals proliferation status
	 * of tumor populations and potential therapeutic vulnerabilities.
	 */
	int indice_de_la_fase_actual;
	
	/**
	 * @brief Time accumulated in the current phase
	 * @details Tracks how long the cell has been in its current phase, critical for
	 * determining phase transitions. Cancer cells often show altered phase durations,
	 * particularly shortened G1 phases or prolonged S phases.
	 */
	double tiempo_acumulado_en_la_fase;
	
	/**
	 * @brief Flag indicating if cell is ready to divide
	 * @details Set to true when a cell completes all cycle phases and is ready to divide.
	 * In cancer modeling, tracks proliferation events and controls tumor growth rates.
	 */
	bool flagged_para_dividirse;
	
	/**
	 * @brief Flag indicating if cell should be removed
	 * @details Set to true when a cell is committed to a death pathway and should be removed.
	 * Critical for modeling cancer cell death due to therapy, stress, or immune killing.
	 */
	bool flagged_para_remover;
	
	/**
	 * @brief Transition rates between cell cycle phases
	 * @details Matrix of transition rates between phases, potentially modified by
	 * microenvironmental factors. In cancer, these rates are often dysregulated
	 * due to oncogenic mutations or microenvironmental stressors.
	 */
	std::vector< std::vector<double> > tasas_de_transicion;
    
    /**
     * @brief Random rate value used in stochastic transitions
     * @details Used for probabilistic phase transitions, adding cell-to-cell
     * variability. Enables modeling of heterogeneous division rates observed
     * in cancer cell populations.
     */
    double tasa_aleatoria;
		
	/**
	 * @brief Default constructor
	 * @details Creates a new cycle instance with NULL model pointer and default values.
	 * Must be followed by sync_con_ciclo_modelo() to connect to a specific model.
	 */
	Ciclo();
	
	/**
	 * @brief Synchronizes this cycle instance with a cycle model
	 * @param cm Reference to the cycle model to be used
	 * @details Links this cycle instance to a specific model and copies transition rates.
	 * In cancer modeling, allows assignment of specific cycle types to different
	 * cancer cell populations.
	 */
	void sync_con_ciclo_modelo( Ciclo_Modelo& cm );
	
	/**
	 * @brief Advances the cell through its cycle for a time step
	 * @param volumen Reference to the cell's volume
	 * @param dt Time step size
	 * @param c_tasas_de_transicion Current transition rates, potentially modified by conditions
	 * @param mp Death parameters affecting cycle progression
	 * @details Progresses the cell through its cycle for the given time step,
	 * updating phase, time accumulated, and potentially setting flags for division or death.
	 * In cancer research, models how cancer cells progress through aberrant cycles
	 * and respond to microenvironmental conditions or therapies.
	 */
	void avanzar_en_el_ciclo( Volumen& volumen, double dt, std::vector< std::vector<double> >& c_tasas_de_transicion, Muerte_parametros& mp );
	
	/**
	 * @brief Checks if volume should be updated in current phase
	 * @return Boolean indicating if volume update is needed
	 * @details Different cycle phases may have different volume update requirements.
	 * In cancer modeling, reflects how growth patterns vary across cell cycle phases,
	 * with cancer cells often showing dysregulated growth in specific phases.
	 */
	bool actualizar_volumen();
	
	/**
	 * @brief Gets the current phase's transition rate
	 * @return Transition rate from current phase
	 * @details Retrieves the base transition rate for the current phase.
	 * In cancer models, these rates often reflect oncogenic alterations
	 * that accelerate specific phase transitions.
	 */
	double tasa_de_transicion();
    
    /**
     * @brief Updates transition rate between specified phases
     * @param fase_actual Current phase index
     * @param fase_siguiente Target phase index
     * @return Reference to the specific transition rate
     * @details Provides access to modify specific transition rates based on
     * environmental conditions. In cancer models, enables representing how
     * hypoxia, acidosis, or therapeutics alter specific cycle transitions.
     */
    double& actualizar_mis_tasas_de_transicion(int fase_actual, int fase_siguiente);
    
    /**
     * @brief Gets reference to the current phase object
     * @return Reference to current Fase object
     * @details Provides access to the detailed properties of the current cell cycle phase.
     * In cancer research, allows examination of phase-specific behaviors that may
     * be dysregulated or targeted by therapies.
     */
	Fase& fase_actual();
};

#endif
