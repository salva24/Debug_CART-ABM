/**
 * @file Muerte.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of cell death mechanisms for cancer simulation
 *
 * @details This file implements the Muerte (Death) class methods, which handle
 * the initialization, management, and execution of cell death processes in the
 * cancer simulation. The implementation supports multiple death pathways with
 * stochastic transitions based on configurable rates.
 * 
 * In cancer research, this implementation is significant because:
 * - It enables modeling of differential death responses in cancer cells
 * - It supports stochastic death decisions influenced by microenvironmental factors
 * - It allows for multiple death pathways (apoptosis, necrosis) with distinct behaviors
 * - It facilitates the simulation of death resistance, a hallmark of cancer
 * - It provides mechanisms to model treatment-induced cell death
 * 
 * The implementation uses random number generation to model the probabilistic
 * nature of cell death decisions, which is critical for realistic cancer modeling.
 * 
 * Inbound Dependencies: Muerte.h, Ciclo_Modelo.h, Muerte_Parametros.h, Random.h
 *
 * Outbound Dependencies: None
 * 
 * Usage: Methods are called by Celula instances to manage their death processes
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#include "Muerte.h"

/**
 * @brief Constructor for the Muerte (Death) class
 * @details Initializes a new death management object with empty collections
 * for death rates, cycle models, and parameters. Sets the cell as initially
 * alive (not dead).
 * 
 * In cancer research, this initialization represents the starting state of
 * a cancer cell that has not yet entered any death pathway, which is the
 * default state for newly created cells in the simulation.
 */
Muerte::Muerte(){

	tasas.resize(0);
	ciclos.resize(0);
	parametros.resize(0);

	muerta = false;
	indice_del_ciclo_de_muerte_actual= 0;

	return;

}

/**
 * @brief Adds a new death cycle with default parameters
 * @details Registers a new death cycle with the specified transition rate and
 * cycle model, using default parameters. This method allows the simulation to
 * incorporate different death mechanisms with their specific transition dynamics.
 * 
 * In cancer research, this enables modeling of various death pathways like
 * apoptosis or necrosis, each with characteristic rates that may be altered
 * in cancer cells. For example, reduced apoptosis rates can model the
 * apoptotic resistance common in many cancer types.
 * 
 * @param tasa The base rate at which cells enter this death cycle
 * @param pModelo Pointer to the cell cycle model defining the death process
 * @return The index of the newly added death cycle
 */
int Muerte::agregar_ciclo_de_muerte(double tasa, Ciclo_Modelo* pModelo){

	tasas.push_back(tasa);
	ciclos.push_back(pModelo);
	parametros.resize(tasas.size());

	return tasas.size() - 1;

}

/**
 * @brief Adds a new death cycle with custom parameters
 * @details Registers a new death cycle with the specified transition rate,
 * cycle model, and custom parameters. This method provides more detailed
 * control over death mechanisms by allowing customization of parameters
 * specific to each death pathway.
 * 
 * In cancer research, this customization is crucial for modeling the
 * heterogeneity of cancer cell responses to death signals. Different
 * cancer types and even subpopulations within the same tumor can exhibit
 * varied death parameters, which this method allows researchers to simulate.
 * 
 * @param tasa The base rate at which cells enter this death cycle
 * @param pModelo Pointer to the cell cycle model defining the death process
 * @param muerte_parametros Custom parameters for this specific death mechanism
 * @return The index of the newly added death cycle
 */
int Muerte::agregar_ciclo_de_muerte(double tasa, Ciclo_Modelo* pModelo, Muerte_parametros& muerte_parametros){

	tasas.push_back(tasa);
	ciclos.push_back(pModelo);
	parametros.push_back(muerte_parametros);

	return tasas.size() - 1;

}

/**
 * @brief Finds death cycle index by code identifier
 * @details Searches through the registered death cycles to find one with
 * a matching code identifier. This method enables referencing specific
 * death mechanisms by their unique codes during simulation.
 * 
 * In cancer research, this functionality allows the simulation to
 * target specific death pathways (e.g., apoptosis vs. necrosis)
 * when modeling responses to treatments or environmental stressors
 * that may trigger particular death mechanisms.
 * 
 * @param codigo The code identifier for the death cycle to find
 * @return The index of the matching death cycle, or 0 if not found
 */
int Muerte::encontrar_indice_del_ciclo_de_muerte(int codigo){

	for( unsigned int i=0 ; i < ciclos.size() ; i++ ){

		if( ciclos[i]->codigo == codigo ){
			return i;
		}
	}
	return 0;

}

/**
 * @brief Finds death cycle index by name
 * @details Searches through the registered death cycles to find one with
 * a matching name. This method provides a more human-readable way to
 * reference specific death mechanisms during simulation.
 * 
 * In cancer research, this allows for intuitive referencing of death
 * pathways by descriptive names like "apoptosis" or "necrosis" when
 * analyzing simulation results or configuring cell responses to
 * different stimuli, enhancing the interpretability of the model.
 * 
 * @param nombre The name of the death cycle to find
 * @return The index of the matching death cycle, or 0 if not found
 */
int Muerte::encontrar_indice_del_ciclo_de_muerte(std::string nombre){

	for( unsigned int i=0 ; i < ciclos.size() ; i++ ){

		if( ciclos[i]->nombre == nombre ){
			return i;
		}
	}
	return 0;

}

/**
 * @brief Determines if the cell should enter a death process
 * @details Stochastically evaluates whether the cell should begin dying
 * based on current death rates and elapsed time. The method checks each
 * registered death pathway in sequence, using random number generation
 * to model the probabilistic nature of cell death decisions.
 * 
 * In cancer research, this stochastic approach is essential for realistic modeling:
 * - It captures the probabilistic nature of cell death in biological systems
 * - It allows for modeling reduced death probabilities in cancer cells
 * - It enables simulation of heterogeneous responses within cell populations
 * - It can reflect how microenvironmental factors influence death likelihood
 * - It supports modeling of treatment efficacy through modified death rates
 * 
 * The method returns immediately if the cell is already dead, preventing
 * multiple death pathways from being activated simultaneously.
 * 
 * @param dt The time step for the current simulation iteration
 * @return True if the cell should begin a death process, false otherwise
 */
bool Muerte::chequear_muerte(double dt){
	if( muerta == true ){
		return false;
	}

	// If the cell is alive, evaluate all the death rates for each registered death type.
	unsigned int i = 0;
	while( !muerta && i < tasas.size()){
        double numaleat = rng->RandomNumber();
		
        if( numaleat < tasas[i]*dt ){
            
			// update the Death data structure
			muerta = true;
			indice_del_ciclo_de_muerte_actual = i;
			return muerta;
		}
		i++;
	}

	return muerta;

}

/**
 * @brief Initiates a specific death process
 * @details Explicitly starts a particular death cycle for the cell and
 * marks it as dying. Unlike chequear_muerte which stochastically determines
 * death, this method directly triggers a specific death pathway.
 * 
 * In cancer research, this method is critical for modeling:
 * - Targeted cell killing by immune cells or therapeutic agents
 * - Forced apoptosis in response to specific signals
 * - Necrotic death due to severe environmental stress
 * - Experimental interventions that directly induce cell death
 * - Programmed cell death during development or tissue homeostasis
 * 
 * This deterministic approach complements the stochastic death checks,
 * allowing for both random and directed cell death in the simulation.
 * 
 * @param indice_ciclo_de_muerte The index of the death cycle to initiate
 */
void Muerte::comenzar_muerte(int indice_ciclo_de_muerte){

	muerta = true;
	indice_del_ciclo_de_muerte_actual = indice_ciclo_de_muerte;
	return;

}

/**
 * @brief Gets the currently active death cycle model
 * @details Retrieves the cycle model for the current death process.
 * This method provides access to the specific death pathway the cell
 * is currently undergoing, allowing the simulation to determine how
 * the cell progresses through its death.
 * 
 * In cancer research, this is important for modeling:
 * - Different temporal patterns of cell death (rapid vs. slow)
 * - Morphological changes during death processes
 * - Phase-specific behaviors during cell death
 * - Interactions with neighboring cells during death
 * - Biochemical changes associated with specific death pathways
 * 
 * @return Reference to the current death cycle model
 */
Ciclo_Modelo& Muerte::ciclo_actual(void){

	return *ciclos[indice_del_ciclo_de_muerte_actual];

}

/**
 * @brief Gets the parameters for the current death cycle
 * @details Retrieves the parameter set controlling the current death process.
 * This method provides access to the specific configuration parameters
 * for the active death mechanism, allowing the simulation to apply the
 * appropriate behaviors during cell death.
 * 
 * In cancer research, these parameters are crucial for modeling:
 * - Cell-type specific death characteristics
 * - Microenvironmental influences on death processes
 * - Treatment-modified death behaviors
 * - Tumor-specific alterations to death pathways
 * - Experimental interventions that modify death parameters
 * 
 * @return Reference to the current death parameters
 */
Muerte_parametros& Muerte::parametros_actuales(void){

	return parametros[indice_del_ciclo_de_muerte_actual];
}
