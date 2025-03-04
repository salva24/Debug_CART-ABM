/*! @file Ciclos_estandares.h
 *
 * @brief Standard cell cycle implementations for cancer simulation modeling.
 * 
 * @author Luciana Melina Luque
 *
 * @details
 * This file defines standard cell cycle models that are used in cancer simulations,
 * including proliferative cycles (Ki67, vida/basic life) and death cycles 
 * (apoptosis, necrosis). These models are critical for accurately simulating 
 * tumor growth dynamics and cell death processes in cancer research.
 *
 * Inbound Dependencies:
 *   - Constantes.h - Defines constants used for cell cycle phases
 *   - Ciclo_Modelo.h - Base class for cell cycle models
 *   - Muerte_Parametros.h - Parameters for cell death models
 * 
 * Outbound Dependencies:
 *   - Main.cpp - Uses these cycles for simulation initialization
 *   - Celula.h - Cells use these cycle models for progression
 * 
 * @usage Initialize with crear_ciclo_celular_estandar() and crear_ciclo_de_muerte_estandar()
 * 
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#include "Constantes.h"
#include "Ciclo_Modelo.h"
#include "Muerte_Parametros.h"


/*! \brief Global cell cycle model instances used throughout the simulation */
extern Ciclo_Modelo Ki67, vida, necrosis, apoptosis;

/*! \brief Global death parameter instances used throughout the simulation */
extern Muerte_parametros necrosis_parametros, apoptosis_parametros;

/*! \brief Flags indicating whether standard cycles have been initialized */
extern bool ciclo_celular_estandar_inicializado;
extern bool ciclo_celular_de_muerte_inicializado;

/*! \brief Entry function for the Ki67 positive phase
 *
 * This function is called when a cell enters the Ki67 positive phase of the cell cycle.
 * It adjusts the target volumes to simulate cell growth during proliferation,
 * which is critical for modeling tumor growth dynamics.
 *
 * \param volumen Cell's volume parameters that will be modified
 * \param mp Death parameters (unused but required by function signature)
 */
void Ki67_fase_positiva_funcion_de_entrada(Volumen& volumen, Muerte_parametros& mp); // done

/*! \brief Entry function for the basic life phase
 *
 * Called when a cell enters the basic life cycle phase. Sets volume
 * targets for growth during the cell's living phase.
 *
 * \param volumen Cell's volume parameters that will be modified
 * \param mp Death parameters (unused but required by function signature)
 */
void fase_viva_funcion_de_entrada(Volumen& volumen, Muerte_parametros& mp);

/*! \brief Entry function for the necrosis phase
 *
 * Called when a cell becomes necrotic. Adjusts volume parameters to simulate
 * cell swelling during necrotic death, which occurs in hypoxic/nutrient-deprived 
 * regions of tumors.
 *
 * \param volumen Cell's volume parameters that will be modified
 * \param mp Death parameters used to control the necrosis process
 */
void standard_necrosis_funcion_de_entrada( Volumen& volumen, Muerte_parametros& mp );  // done

/*! \brief Entry function for the apoptosis phase
 *
 * Called when a cell undergoes apoptosis (programmed cell death).
 * Adjusts volume parameters to simulate cell shrinkage during apoptotic death.
 * Critical for modeling therapeutic responses and natural cell turnover in tumors.
 *
 * \param volumen Cell's volume parameters that will be modified
 * \param mp Death parameters used to control the apoptosis process
 */
void standard_apoptosis_funcion_de_entrada( Volumen& volumen, Muerte_parametros& mp );

/*! \brief Entry function for the lysis phase
 *
 * Called when a cell enters the lysed state. Adjusts volume parameters to
 * simulate cellular breakdown after death.
 *
 * \param volumen Cell's volume parameters that will be modified
 * \param mp Death parameters used to control the lysis process
 */
void standard_lysis_funcion_de_entrada( Volumen& volumen, Muerte_parametros& mp ); // done

/*! \brief Arrest function for the necrosis phase
 *
 * This function determines whether a necrotic cell should proceed to the next phase.
 * The arrest occurs based on volume thresholds, representing cell bursting.
 *
 * \param volumen Cell's volume parameters
 * \param mp Death parameters
 * \return true if the cell should arrest (not transition), false otherwise
 */
bool standard_necrosis_funcion_de_arrest( Volumen& volumen, Muerte_parametros& mp ); // done

/*! \brief Creates and initializes standard cell cycles (Ki67, vida)
 *
 * Initializes proliferative cell cycle models used in the simulation.
 * These cycles model how cancer cells progress through different growth phases,
 * which is essential for tumor growth dynamics research.
 *
 * \return true if initialization was successful, false if already initialized
 */
bool crear_ciclo_celular_estandar( void );

/*! \brief Creates and initializes standard death cycles (necrosis, apoptosis)
 *
 * Initializes cell death models used in the simulation. These cycles model how cancer
 * cells die through different mechanisms, which is crucial for understanding
 * tumor behavior and treatment responses.
 *
 * \return true if initialization was successful, false if already initialized
 */
bool crear_ciclo_de_muerte_estandar( void );
