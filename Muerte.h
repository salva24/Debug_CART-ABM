/**
 * @file Muerte.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Cell death modeling system for the cancer simulation
 *
 * @details This file defines the Muerte (Death) class, which manages different
 * cell death mechanisms within the cancer simulation. The class handles multiple
 * death cycles (e.g., apoptosis, necrosis) with associated parameters and transition
 * rates.
 * 
 * In cancer research, accurately modeling cell death is critical because:
 * - Dysregulation of cell death is a hallmark of cancer
 * - Cancer cells often develop resistance to apoptosis
 * - Necrotic regions in tumors influence tumor microenvironment
 * - Different cancer types show characteristic patterns of cell death
 * - Therapeutic approaches often target death pathways
 * 
 * This module enables simulation of various death mechanisms and their
 * impact on tumor growth, progression, and response to environmental stressors.
 * 
 * Inbound Dependencies: Ciclo_Modelo.h, Muerte_Parametros.h, Random.h
 *
 * Outbound Dependencies: Celula.h, Ciclo_Modelo.cpp
 * 
 * Usage: Instantiated within each cell to manage its death processes
 *
 * @license: Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#ifndef __MUERTE_H__
#define __MUERTE_H__

#include<vector>
#include "Ciclo_Modelo.h"
#include "Muerte_Parametros.h"
#include "Random.h"

extern RNG *rng;

class Muerte{
	
	public:
		/**
		 * @brief Transition rates for each death cycle
		 * @details Vector of rates at which cells may enter each death cycle.
		 * These rates can be adjusted dynamically based on environmental conditions,
		 * allowing for modeling of environmentally-induced cell death in cancer.
		 */
		std::vector<double> tasas;
		
		/**
		 * @brief Collection of cell death cycle models
		 * @details Stores pointers to the cycle models (finite state machines) that
		 * define how cells progress through different death mechanisms like apoptosis
		 * and necrosis, which are critical for modeling cancer cell populations.
		 */
		std::vector<Ciclo_Modelo*> ciclos;
		
		/**
		 * @brief Parameters for each death cycle
		 * @details Contains the specific parameters that control each death mechanism's
		 * behavior, allowing for customization of death processes to model different
		 * cancer types and their characteristic patterns of cell death.
		 */
		std::vector<Muerte_parametros> parametros;
		
		/**
		 * @brief Flag indicating if the cell is dead
		 * @details When true, indicates the cell has entered an irreversible death process.
		 * In cancer research, tracking dead cells is important for modeling tumor necrotic
		 * cores and response to treatments.
		 */
		bool muerta;
		
		/**
		 * @brief Index of the currently active death cycle
		 * @details Tracks which death mechanism the cell is currently undergoing.
		 * Allows the simulation to apply the appropriate death behavior for modeling
		 * various causes of cancer cell death.
		 */
		int indice_del_ciclo_de_muerte_actual;
		
		/**
		 * @brief Default constructor
		 * @details Initializes a new Muerte object with default settings.
		 * Sets up the death tracking system for a new cell in the cancer simulation.
		 */
		Muerte();
		
		/**
		 * @brief Adds a new death cycle to the cell
		 * @details Registers a new death cycle with a specified transition rate.
		 * Enables modeling of multiple death pathways relevant to cancer research,
		 * such as apoptosis or necrosis, each with their own transition dynamics.
		 * 
		 * @param tasa The base rate at which cells enter this death cycle
		 * @param pModelo Pointer to the cell cycle model defining the death process
		 * @return The index of the newly added death cycle
		 */
		int agregar_ciclo_de_muerte(double tasa, Ciclo_Modelo* pModelo);
		
		/**
		 * @brief Adds a new death cycle with custom parameters
		 * @details Registers a new death cycle with both a transition rate and custom parameters.
		 * Allows for detailed customization of death processes to model cancer-specific
		 * cell death mechanisms and their unique behaviors.
		 * 
		 * @param tasa The base rate at which cells enter this death cycle
		 * @param pModelo Pointer to the cell cycle model defining the death process
		 * @param muerte_parametros Custom parameters for this specific death mechanism
		 * @return The index of the newly added death cycle
		 */
		int agregar_ciclo_de_muerte(double tasa, Ciclo_Modelo* pModelo, Muerte_parametros& muerte_parametros);
		
		/**
		 * @brief Finds death cycle index by code
		 * @details Retrieves the index of a death cycle based on its code identifier.
		 * Used in cancer models to reference specific death mechanisms during simulation.
		 * 
		 * @param codigo The code identifier for the death cycle
		 * @return The index of the death cycle, or -1 if not found
		 */
		int encontrar_indice_del_ciclo_de_muerte(int codigo);
		
		/**
		 * @brief Finds death cycle index by name
		 * @details Retrieves the index of a death cycle based on its name.
		 * Allows cancer researchers to reference death mechanisms by descriptive names
		 * like "apoptosis" or "necrosis" in the simulation.
		 * 
		 * @param nombre The name of the death cycle to find
		 * @return The index of the death cycle, or -1 if not found
		 */
		int encontrar_indice_del_ciclo_de_muerte(std::string nombre);
		
		/**
		 * @brief Determines if the cell should enter a death process
		 * @details Stochastically evaluates whether the cell should begin dying based on
		 * current death rates and elapsed time. This method models the probabilistic nature
		 * of cell death decisions in cancer, which are influenced by various factors.
		 * 
		 * In cancer research, this represents how cells respond to death signals,
		 * with altered probabilities reflecting cancer cells' characteristic
		 * resistance to death triggers.
		 * 
		 * @param dt The time step for the current simulation iteration
		 * @return True if the cell should begin a death process, false otherwise
		 */
		bool chequear_muerte(double dt);
		
		/**
		 * @brief Initiates a specific death process
		 * @details Starts a particular death cycle for the cell and marks it as dying.
		 * This method is critical for modeling how cancer cells enter different death pathways
		 * in response to environmental stressors, therapeutic agents, or internal signals.
		 * 
		 * @param indice_ciclo_de_muerte The index of the death cycle to initiate
		 */
		void comenzar_muerte(int indice_ciclo_de_muerte);
		
		/**
		 * @brief Gets the currently active death cycle model
		 * @details Retrieves the cycle model for the current death process.
		 * Used to determine how the cell progresses through its death in cancer simulations,
		 * modeling various temporal patterns of cell death observed in different cancer types.
		 * 
		 * @return Reference to the current death cycle model
		 */
		Ciclo_Modelo& ciclo_actual(void);
		
		/**
		 * @brief Gets the parameters for the current death cycle
		 * @details Retrieves the parameter set controlling the current death process.
		 * Allows access to specific configuration of the active death mechanism,
		 * which may vary across different cancer types or microenvironmental conditions.
		 * 
		 * @return Reference to the current death parameters
		 */
		Muerte_parametros& parametros_actuales(void);
	
	
	
};


#endif
