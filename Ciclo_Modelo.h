/**
 * @file Ciclo_Modelo.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Cell Cycle Model abstract class for cancer simulation
 *
 * @details This file defines the Ciclo_Modelo class, which represents the cell cycle model
 * used in cancer research simulations. It manages phases, transitions between phases,
 * and cell division/death events based on phase transitions.
 * 
 * Inbound Dependencies: Fase.h, Fase_Link.h, Random.h
 *
 * Outbound Dependencies: Celula.h, Volumen.h, Muerte_parametros.h
 * 
 * @usage Instantiated by Celula instances to model their division cycle
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#ifndef __CICLO_MODELO_H__
#define __CICLO_MODELO_H__

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>

#include "Fase.h"
#include "Fase_Link.h"
#include "Random.h"

extern RNG *rng;


class Ciclo;

/**
 * @class Ciclo_Modelo
 * @brief Models the cell cycle with phases and transitions
 * @details
 * The Ciclo_Modelo class implements a finite state machine representing
 * the cell cycle with multiple phases and transitions between phases.
 * 
 * In cancer research, this model is critical for:
 * - Representing heterogeneous cell cycle dynamics in tumor populations
 * - Modeling different proliferation rates across tumor regions
 * - Simulating the effects of drugs targeting specific cell cycle phases
 * - Capturing growth arrest in response to environmental factors like hypoxia
 * 
 * Each phase is represented by a Fase object, and transitions between phases
 * are governed by transition rates and optional arrest functions that
 * can halt transitions based on cell volume or death parameters.
 */
class Ciclo_Modelo{
	private:
		/** @brief Maps to find the correct index for phase transitions */
		std::vector< std::unordered_map<int,int> > mapas_de_indice_inverso; 
 
 	public:
 		/** @brief Name of the cell cycle model */
 		std::string nombre;
		/** @brief Time units used for transition rates (e.g., "min") */
		std::string unidades_tiempo;
 		/** @brief Unique identifier code for the model */
 		int codigo; 
 		/** @brief Collection of all phases in the cell cycle */
 		std::vector<Fase> fases; 
 		/** @brief Collection of links between phases */
 		std::vector< std::vector<Fase_Link> > fase_links; 
		/** @brief Transition rates between phases */
		std::vector< std::vector<double> > tasas_de_transicion;

		/**
		 * @brief Default constructor
		 * @details Initializes an empty cell cycle model with default name and units
		 */
		Ciclo_Modelo(); 
		
		/**
		 * @brief Adds a new phase to the cell cycle model
		 * @param codigo Unique identifier for the phase
		 * @param nombre Descriptive name of the phase
		 * @return Index of the newly added phase
		 */
		int agregar_fase( int codigo, std::string nombre);
		
		/**
		 * @brief Creates a link between two phases
		 * @param indice_fase_inicial Index of the source phase
		 * @param indice_fase_final Index of the destination phase
		 * @param funcion_arrest Optional function that can arrest the transition
		 * @return Index of the newly created link
		 * @details In cancer modeling, arrest functions can represent growth controls
		 * like contact inhibition or environmental stress responses
		 */
		int agregar_link(int indice_fase_inicial, int indice_fase_final, bool (*funcion_arrest)( Volumen& volumen, Muerte_parametros& mp ) );
		
		/**
		 * @brief Gets reference to the transition rate between phases
		 * @param indice_fase_inicial Index of the source phase
		 * @param indice_fase_final Index of the destination phase
		 * @return Reference to the transition rate value
		 */
		double& tasa_de_transicion( int indice_fase_inicial, int indice_fase_final);
		
		/**
		 * @brief Gets the phase link between two phases
		 * @param indice_fase_inicial Index of the source phase
		 * @param indice_fase_final Index of the destination phase
		 * @return Reference to the phase link
		 */
		Fase_Link& fase_link(int indice_fase_inicial, int indice_fase_final);
		
		/**
		 * @brief Displays the cell cycle model structure
		 * @param os Output stream to write to
		 * @return Reference to the output stream
		 */
	    std::ostream& mostrar_ciclo( std::ostream& os );
	    
	    /**
	     * @brief Advances the cell in its current cell cycle phase
	     * @param flagged_para_remover Set to true if cell should be removed
	     * @param flagged_para_dividirse Set to true if cell should divide
	     * @param indice_de_la_fase_actual Current phase index (updated if phase changes)
	     * @param tiempo_acumulado_en_la_fase Time spent in current phase (updated)
	     * @param volumen Cell volume parameters
	     * @param dt Time step for simulation
	     * @param c_tasas_de_transicion Current transition rates (possibly modified by microenvironment)
	     * @param mp Death parameters for the cell
	     * @details 
	     * This is a key method for cancer simulation as it determines:
	     * - When and if cells divide (tumor growth)
	     * - When cells transition between phases (affecting drug sensitivity)
	     * - How microenvironmental factors affect cycle progression
	     */
	    void avanzar_en_el_modelo(bool& flagged_para_remover, bool& flagged_para_dividirse, int& indice_de_la_fase_actual, double& tiempo_acumulado_en_la_fase, Volumen& volumen, double dt, std::vector< std::vector<double> >& c_tasas_de_transicion, Muerte_parametros& mp);
	    
	    /**
	     * @brief Gets the transition rate for a phase
	     * @param indice_fase_inicial Index of the phase
	     * @return Transition rate from the specified phase
	     */
	    double get_tasa_de_transicion( int indice_fase_inicial);
	    
	    /**
	     * @brief Checks if volume should be updated in the given phase
	     * @param indice_fase_inicial Index of the phase to check
	     * @return True if volume should be updated, false otherwise
	     */
	    bool get_actualizar_volumen(int indice_fase_inicial);
	    
	    /**
	     * @brief Finds the index of a phase by its code
	     * @param codigo The code to search for
	     * @return Index of the phase with the given code
	     */
	    int encontrar_indice_de_la_fase(int codigo);
	    
	    /**
	     * @brief Gets the inverse map index for a phase transition
	     * @param fase_uno Source phase index
	     * @param fase_dos Destination phase index
	     * @return Index in the internal map
	     */
        int get_indice_de_mapa_inverso(int fase_uno, int fase_dos);
};

#endif
