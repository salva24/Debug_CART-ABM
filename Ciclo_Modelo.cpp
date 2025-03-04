/**
 * @file Ciclo_Modelo.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of the Cell Cycle Model class for cancer simulation
 *
 * @details
 * This file implements the Ciclo_Modelo class methods, including phase creation,
 * phase transition management, and cell cycle progression logic. The cell cycle
 * model is represented as a finite state machine with phases and transitions.
 * 
 * In cancer research, cell cycle models are essential for studying:
 * - Differential proliferation rates in tumor regions
 * - Effects of targeting specific cycle phases with therapeutics
 * - Responses to growth factors and inhibitors
 * - Cell cycle arrest due to stress conditions
 * 
 * Inbound Dependencies: Ciclo_Modelo.h
 *
 * Outbound Dependencies: Volumen.h, Muerte_parametros.h
 * 
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#include "Ciclo_Modelo.h"

/**
 * @brief Default constructor for the Ciclo_Modelo class
 * @details
 * Initializes a new cell cycle model with:
 * - Empty phase collections
 * - Default name "Sin nombre"
 * - Default time units "min"
 * - Default code 9999
 * 
 * The model is not functional until phases and transitions are added.
 */
Ciclo_Modelo::Ciclo_Modelo(){

	mapas_de_indice_inverso.resize(0);

	nombre = "Sin nombre";
	unidades_tiempo = "min";
	codigo = 9999;

	fases.resize(0);
	fase_links.resize(0);
	tasas_de_transicion.resize(0);


	return;

}

/**
 * @brief Adds a new phase to the cell cycle model
 * @details
 * Creates a new phase in the model and initializes all related data structures.
 * In cancer modeling, phases typically represent key cell cycle stages like:
 * - G0 (quiescence)
 * - G1 (growth and preparation for DNA synthesis)
 * - S (DNA synthesis)
 * - G2 (preparation for mitosis)
 * - M (mitosis/cell division)
 * 
 * @param codigo Unique identifier for the phase
 * @param nombre Descriptive name of the phase
 * @return Index of the newly added phase
 */
int Ciclo_Modelo::agregar_fase( int codigo, std::string nombre  ){

	int n = fases.size();

	//redimensionar las estructuras de datos
	fases.resize(n+1);
	fase_links.resize(n+1);
	tasas_de_transicion.resize(n+1);
	fase_links[n].resize(0);
	tasas_de_transicion[n].resize(0);

	mapas_de_indice_inverso.resize(n+1);
	mapas_de_indice_inverso[n].clear();

	//actualizar la fase n
	fases[n].codigo = codigo;
	fases[n].indice = n;
	fases[n].nombre.assign( nombre );

	return n;

}

/**
 * @brief Creates a link between two phases in the cell cycle
 * @details
 * Establishes a transition pathway between phases and assigns an optional
 * arrest function that can block the transition under certain conditions.
 * 
 * In cancer research, phase transitions can be modulated by:
 * - Growth factors in the microenvironment
 * - Stress conditions like hypoxia
 * - Checkpoint proteins like p53
 * - Drug effects on specific cycle phases
 * 
 * @param indice_fase_inicial Index of the source phase
 * @param indice_fase_final Index of the destination phase
 * @param funcion_arrest Function pointer to a condition that can arrest the transition
 * @return Index of the newly created link
 */
int Ciclo_Modelo::agregar_link(int indice_fase_inicial, int indice_fase_final, bool (*funcion_arrest)( Volumen& volumen, Muerte_parametros& mp ) ){


	//Primero redimensionar los links
	int n =	fase_links[indice_fase_inicial].size();
	fase_links[indice_fase_inicial].resize(n+1);
	tasas_de_transicion[indice_fase_inicial].resize(n+1);


	//ahora, actualizo el nuevo link de las fases

	fase_links[indice_fase_inicial][n].indice_fase_inicial = indice_fase_inicial;
	fase_links[indice_fase_inicial][n].indice_fase_final = indice_fase_final;
	fase_links[indice_fase_inicial][n].funcion_arrest = funcion_arrest;

	//ahora, actualizo el mapa inverso
	mapas_de_indice_inverso[indice_fase_inicial][indice_fase_final] = n;


	return n;

}

/**
 * @brief Gets a reference to the transition rate between two phases
 * @details
 * Retrieves the rate at which cells transition from one phase to another.
 * These rates determine the tumor growth kinetics and can be adjusted to
 * model different cancer types or responses to treatment.
 * 
 * @param indice_fase_inicial Index of the source phase
 * @param indice_fase_final Index of the destination phase
 * @return Reference to the transition rate value
 */
double& Ciclo_Modelo::tasa_de_transicion( int indice_fase_inicial, int indice_fase_final){

	return tasas_de_transicion[indice_fase_inicial][mapas_de_indice_inverso[indice_fase_inicial][indice_fase_final]];

}

/**
 * @brief Gets the transition rate for a specific phase
 * @details
 * Retrieves the total transition rate from a specific phase, which determines
 * how quickly cells exit this phase. In cancer modeling, slower phases like
 * G1 or G0 often have lower rates, while faster phases have higher rates.
 * 
 * @param indice_fase_inicial Index of the phase
 * @return Transition rate from the specified phase
 */
double Ciclo_Modelo::get_tasa_de_transicion( int indice_fase_inicial){

	double tasa;
	tasa = 0.0;
	for (unsigned int k=0; k < fase_links[indice_fase_inicial].size(); k++){
		tasa = tasas_de_transicion[indice_fase_inicial][k];
	}
	return tasa;

}

/**
 * @brief Checks if volume should be updated in a specific phase
 * @details
 * Some cell cycle phases involve volume changes (like growth in G1),
 * while others don't. This method checks if the specified phase
 * should update cell volume parameters.
 * 
 * In cancer modeling, modified growth rates in specific phases
 * can represent abnormal proliferation characteristics.
 * 
 * @param indice_fase_inicial Index of the phase to check
 * @return True if volume should be updated, false otherwise
 */
bool Ciclo_Modelo::get_actualizar_volumen(int indice_fase_inicial){

	return fases[indice_fase_inicial].actualizar_volumen;
}

/**
 * @brief Gets the phase link object between two phases
 * @details
 * Retrieves the link object that contains transition properties between phases.
 * This allows access to arrest functions and other properties of the transition.
 * 
 * @param indice_fase_inicial Index of the source phase
 * @param indice_fase_final Index of the destination phase
 * @return Reference to the phase link object
 */
Fase_Link& Ciclo_Modelo::fase_link(int indice_fase_inicial, int indice_fase_final){

	return fase_links[indice_fase_inicial][mapas_de_indice_inverso[indice_fase_inicial][indice_fase_final]];

}

/**
 * @brief Displays the cell cycle model structure
 * @details
 * Creates a text representation of the cell cycle model showing all phases,
 * their connections, and transition rates. Phases that lead to cell division
 * are marked with an asterisk (*).
 * 
 * This is useful for debugging and verifying the correct structure of
 * cancer-specific cell cycle implementations.
 * 
 * @param os Output stream to write to
 * @return Reference to the output stream
 */
std::ostream& Ciclo_Modelo::mostrar_ciclo( std::ostream& os ){

	os << "Ciclo Celular: " << nombre << "(Codigo: " << codigo << ")" << std::endl;
	os << "Fases y Links: (el * denota division celular en esta fase)" << std::endl;
	for(unsigned int i=0; i< fases.size(); i++){

		os << "La fase " << i << " (" << fases[i].nombre << ") ";

		if(fases[i].division_al_final_de_la_fase){
			os << " * ";
		}
		os << "se conecta con: " << std::endl;
		for( unsigned int k=0 ; k < fase_links[i].size() ; k++ )
		{
			int j = fase_links[i][k].indice_fase_final;

			os << "\tla fase " << j << " (" << fases[j].nombre << ") con tasa de transicion " << tasa_de_transicion(i,j) << " " << unidades_tiempo << "^-1;" << std::endl;
		}
		os << std::endl;
	}
	return os;
}

/**
 * @brief Advances a cell through its current cycle phase
 * @details
 * This is the core method that handles cell cycle progression. It:
 * 1. Updates time spent in the current phase
 * 2. Checks for possible transitions to other phases
 * 3. Applies arrest functions that might block transitions
 * 4. Handles cell division/removal events when phases end
 * 5. Executes phase entry functions when transitioning
 * 
 * In cancer research, this method is critical as it models:
 * - Differential proliferation rates across tumor regions
 * - Cell cycle responses to microenvironmental conditions
 * - Cell division events that drive tumor growth
 * - Effects of cycle-specific drugs and treatments
 * 
 * @param flagged_para_remover Flag set if cell should be removed
 * @param flagged_para_dividirse Flag set if cell should divide
 * @param indice_de_la_fase_actual Current phase index (updated if changed)
 * @param tiempo_acumulado_en_la_fase Time in current phase (updated)
 * @param volumen Cell volume parameters
 * @param dt Time step for simulation
 * @param c_tasas_de_transicion Current transition rates
 * @param mp Cell death parameters
 */
void Ciclo_Modelo::avanzar_en_el_modelo(bool& flagged_para_remover, bool& flagged_para_dividirse, int& indice_de_la_fase_actual, double& tiempo_acumulado_en_la_fase, Volumen& volumen, double dt,std::vector< std::vector<double> >& c_tasas_de_transicion, Muerte_parametros& mp){

	int i = indice_de_la_fase_actual;
	tiempo_acumulado_en_la_fase += dt;

	/* Evalúo cada fase linkeada:
	   Avanzo a la fase SI la probabilidad está en el rango.*/
	int j;
	for ( unsigned int k=0; k < fase_links[i].size(); k++){
		j = fase_links[i][k].indice_fase_final;

		bool transicion_parada = false;
		if(fase_links[i][k].funcion_arrest){

			transicion_parada = fase_links[i][k].funcion_arrest(volumen, mp);
		}

		if(!transicion_parada){


			bool continuar_transicion = false;
			if(fase_links[i][k].duracion_fija){

				if(tiempo_acumulado_en_la_fase > 1.0 / c_tasas_de_transicion[i][k]){
					continuar_transicion = true;
				}
			}
			else{

				double prob = c_tasas_de_transicion[i][k]*dt;

				if(rng->RandomNumber() <= prob){
					continuar_transicion = true;
				}
			}

			if (continuar_transicion){

				//chequeo si me tengo que dividir
				if(fases[i].division_al_final_de_la_fase){

					//volumen.dividir();
					flagged_para_dividirse = true;
				}

				if(fases[i].remover_al_final_de_la_fase){

					flagged_para_remover = true;
					return;
				}

				//me muevo a la próxima fase y reseteo el tiempo acumulado
				indice_de_la_fase_actual = j;
				tiempo_acumulado_en_la_fase = 0.0;

				//si la fase tiene alguna función de entrada, la ejecuto
				if(fases[j].funcion_de_entrada){
					fases[j].funcion_de_entrada(volumen, mp);
				}

				return;
			}
		}
	}



	return;

}

/**
 * @brief Finds the index of a phase by its code
 * @details
 * Searches through all phases to find one with the specified code.
 * This is useful when initializing cells with specific cycle phases
 * or when transitioning between different cycle models.
 * 
 * In cancer research, this helps identify particular phases that
 * might be targeted by treatments or affected by mutations.
 * 
 * @param codigo The code to search for
 * @return Index of the phase with the given code, or 0 if not found
 */
int Ciclo_Modelo::encontrar_indice_de_la_fase(int codigo){

	for( unsigned int i=0 ; i < fases.size() ; i++ )
	{
		if( fases[i].codigo == codigo )
		{ return i; }
	}
	return 0;

}

/**
 * @brief Gets the inverse map index for a phase transition
 * @details
 * Retrieves the internal index mapping for a transition between phases.
 * This is used internally to access the correct transition rate and link objects.
 * 
 * @param fase_uno Source phase index
 * @param fase_dos Destination phase index
 * @return Internal index for the phase transition
 */
int Ciclo_Modelo::get_indice_de_mapa_inverso(int fase_uno, int fase_dos){

    return mapas_de_indice_inverso[fase_uno][fase_dos];

}
