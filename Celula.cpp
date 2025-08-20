/**
 * @file Celula.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of the Cell (Celula) class and Lymphocyte (Linfocito) class for cancer simulation
 *
 * @details This file implements the core cell agent behaviors for the cancer simulation model, including:
 *          - Cell cycle progression and volume dynamics
 *          - Cell death processes (apoptosis and necrosis)
 *          - Cell division and proliferation
 *          - Cell-cell interactions and mechanics
 *          - Microenvironment sensing and response to oxygen
 *          - Cell secretion and consumption of biochemical factors
 *          - Lymphocyte-specific behaviors (immune cell functions)
 * 
 * The implementation focuses on modeling cancer-specific behaviors such as:
 *          - Oxygen-dependent proliferation (modeling hypoxia effects in tumors)
 *          - Immune cell recognition of cancer cells via oncoprotein expression
 *          - Cancer cell-immune cell interactions (adhesion and killing)
 *          - Mechanical interactions between cells (critical for tumor growth)
 * 
 * Inbound dependencies Celula.h, Vector.h, Microambiente.h, Constantes.h, Fase.h, Fenotipo.h
 *
 * Outbound dependencies Microambiente.cpp, Muerte.cpp, Ciclo_Modelo.cpp
 *
 * @usage Cell agents are created and managed by the simulation system with their behaviors
 *        updated at each time step through the cell container.
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#include "Celula.h"

/**
 * @brief Default constructor for the Cell class
 * @details Initializes a cell with default parameters and properties. This creates an
 * inactive cell that needs further initialization before being used in simulation.
 * In cancer research, this provides the basic foundation for modeling both normal and
 * cancer cells with their respective behaviors.
 */
Celula::Celula(){

	nombre = "Sin Nombre";
    tipo = 0;
	id = -1;
    indice = -1;
    voxel = -1;
    madre = 0;
    tiempo_desde_el_ultimo_ciclo = 0.0;
    tiempo_desde_la_ultima_mecanica = 0.0;
    hora_del_ultimo_ciclo = 0.0;
    hora_de_la_ultima_mecanica = 0.0;
    adherida = false;
	posicion.x=-0.0;
	posicion.y=0.0;
	posicion.z=0.0;
    celula_adherida = NULL;

    /*Microambiente y oxígeno*/
    es_activa=true;
    microambiente = NULL;
	if( get_microambiente_default() != NULL )
	{ microambiente = get_microambiente_default(); }
	registrar_microambiente(microambiente);


	return;
}

/**
 * @brief Updates the cell volume according to phenotype parameters
 * @details Manages fluid and solid volume components for the cell, including nuclear and 
 * cytoplasmic components. This method is critical for cancer modeling as aberrant cell growth
 * is a hallmark of cancer.
 * 
 * The method handles:
 * - Changes in fluid fraction
 * - Nuclear and cytoplasmic volume changes
 * - Volume-based geometry updates
 * - Calcification fraction changes
 * 
 * @param fenotipo The cell phenotype containing volume parameters
 * @param dt Time step size for the update (in minutes)
 */
void Celula::actualizar_volumen( Fenotipo& fenotipo, double dt)
{
		// if(fenotipo.ciclo.pCiclo_Modelo->nombre != "Vida"){
// 	bool debug= fenotipo.ciclo.actualizar_volumen() ;
// 	auto st="No";
// 	if(debug){
// 		st="Si";
// 	}
// std::cout << "#" << fenotipo.ciclo.pCiclo_Modelo->nombre << " " << fenotipo.volumen.citoplasma_tasa_de_cambio<< st<< tipo<< "#" << std::endl;
		
	// }
	
	fenotipo.volumen.fluido += dt* fenotipo.volumen.fluido_tasa_de_cambio *
	(fenotipo.volumen.target_fraccion_fluido * fenotipo.volumen.total - fenotipo.volumen.fluido);


	//Si el volumen del fluido es negativo, set a 0
	if (fenotipo.volumen.fluido < 0.0){
		fenotipo.volumen.fluido = 0.0;
	}

	fenotipo.volumen.nuclear_fluido = (fenotipo.volumen.nuclear / fenotipo.volumen.total) *
	(fenotipo.volumen.fluido);
	fenotipo.volumen.citoplasmatico_fluido = fenotipo.volumen.fluido - fenotipo.volumen.nuclear_fluido;


	fenotipo.volumen.nuclear_solido += dt* fenotipo.volumen.nucleo_tasa_de_cambio *
	(fenotipo.volumen.target_nucleo_solido - fenotipo.volumen.nuclear_solido);

	if(fenotipo.volumen.nuclear_solido <0.0){
		fenotipo.volumen.nuclear_solido= 0.0;
	}


	fenotipo.volumen.target_citoplasma_solido = fenotipo.volumen.target_relacion_citoplasma_nucleo *
	fenotipo.volumen.target_nucleo_solido;

	fenotipo.volumen.citoplasmatico_solido += dt* fenotipo.volumen.citoplasma_tasa_de_cambio *
	(fenotipo.volumen.target_citoplasma_solido - fenotipo.volumen.citoplasmatico_solido);

	if(fenotipo.volumen.citoplasmatico_solido <0.0){
		fenotipo.volumen.citoplasmatico_solido= 0.0;
	}

	fenotipo.volumen.solido = fenotipo.volumen.nuclear_solido + fenotipo.volumen.citoplasmatico_solido;

	fenotipo.volumen.nuclear = fenotipo.volumen.nuclear_solido + fenotipo.volumen.nuclear_fluido;
	fenotipo.volumen.citoplasmatico = fenotipo.volumen.citoplasmatico_solido + fenotipo.volumen.citoplasmatico_fluido;

	fenotipo.volumen.fraccion_calcificada = dt* fenotipo.volumen.tasa_de_calcificacion *
	(1- fenotipo.volumen.fraccion_calcificada);

	fenotipo.volumen.total = fenotipo.volumen.citoplasmatico + fenotipo.volumen.nuclear;

	fenotipo.volumen.fraccion_de_fluido = fenotipo.volumen.fluido / (1e-16 + fenotipo.volumen.total);


   	fenotipo.geometria.actualizar(fenotipo.volumen);

	return;
}

/**
 * @brief Updates cell cycle and death parameters based on local oxygen concentration
 * @details This method implements the cell's response to oxygen levels in the microenvironment,
 * which is critical for modeling cancer cell behavior in hypoxic conditions. The method:
 * - Identifies the appropriate cell cycle model (Ki67 or simple live/dead)
 * - Retrieves local oxygen concentration
 * - Adjusts proliferation rates based on oxygen availability (modeling hypoxia effects)
 * - Adjusts necrosis rates when oxygen is severely limited
 *
 * This is particularly important for cancer research as hypoxia influences tumor growth,
 * migration, and treatment response.
 * 
 * @param fenotipo The cell phenotype containing cycle and death parameters
 * @param dt Time step size for the update (in minutes)
 */
void Celula::actualizar_parametros_de_celula_y_muerte_con_o2(Fenotipo& fenotipo, double dt){

	if(fenotipo.muerte.muerta == true){
		return;
	}
	static bool indices_iniciados = false;
	static int indice_fase_inicial;
	static int indice_fase_final;
	static int indice_necrosis;

	static int indice_del_oxigeno = microambiente->encontrar_indice_de_densidad("oxigeno");
	if(indices_iniciados == false){

		if(fenotipo.ciclo.pCiclo_Modelo->codigo == Constantes::ciclo_Ki67){
			indice_fase_inicial = fenotipo.ciclo.pCiclo_Modelo->encontrar_indice_de_la_fase(Constantes::Ki67_negativa);
			indice_necrosis = fenotipo.muerte.encontrar_indice_del_ciclo_de_muerte(Constantes::ciclo_de_muerte_necrosis);
			indice_fase_final = fenotipo.ciclo.pCiclo_Modelo->encontrar_indice_de_la_fase(Constantes::Ki67_positiva_premitotica);
			indices_iniciados = true;

		}

		if( fenotipo.ciclo.pCiclo_Modelo->codigo == Constantes::ciclo_vida ){

			indice_fase_inicial = fenotipo.ciclo.pCiclo_Modelo->encontrar_indice_de_la_fase(Constantes::viva);
			indice_necrosis = fenotipo.muerte.encontrar_indice_del_ciclo_de_muerte(Constantes::ciclo_de_muerte_necrosis);
			indice_fase_final = fenotipo.ciclo.pCiclo_Modelo->encontrar_indice_de_la_fase(Constantes::viva);
			indices_iniciados = true;
		}
	}
	double pO2 = (vector_de_densidades_mas_cercano())[indice_del_oxigeno];
	// pO2 = 30.;//Debug

	double multiplicador = 1.0;
	if(pO2 < parametros.o2_saturacion_para_la_proliferacion){
		multiplicador = (pO2 - parametros.o2_limite_de_proliferacion) /
		(parametros.o2_saturacion_para_la_proliferacion - parametros.o2_limite_de_proliferacion);
	}

	if(pO2 < parametros.o2_limite_de_proliferacion){
		multiplicador= 0.0;
	}

	fenotipo.ciclo.tasas_de_transicion[indice_fase_inicial][indice_fase_final] = multiplicador *
	fenotipo.ciclo.pCiclo_Modelo->tasa_de_transicion(indice_fase_inicial, indice_fase_final);

	multiplicador = 0.0;
	if(pO2 < parametros.o2_necrosis_limite){
		multiplicador = (parametros.o2_necrosis_limite - pO2) /
		(parametros.o2_necrosis_limite - parametros.o2_necrosis_max);
	}
	if(pO2 < parametros.o2_necrosis_max){
		multiplicador = 1.0;
	}
	fenotipo.muerte.tasas[indice_necrosis] = multiplicador * parametros.tasa_necrosis_max;
}

/**
 * @brief Updates cell cycle and death parameters based on oxygen and oncoprotein levels
 * @details Similar to actualizar_parametros_de_celula_y_muerte_con_o2, but also incorporates
 * the effect of oncoprotein expression on cell cycle progression. This method is specific to
 * modeling cancer cells (type 0) and provides a mechanism to simulate oncogene-driven 
 * proliferation, which can be modulated by the microenvironment.
 *
 * Key cancer research aspects:
 * - Combines hypoxia response with oncogene activity
 * - Models how oncogene expression can override normal oxygen dependency
 * - Simulates tumor growth drivers at the molecular and microenvironmental levels
 * 
 * @param fenotipo The cell phenotype containing cycle, death and secretion parameters
 * @param dt Time step size for the update (in minutes)
 */
void Celula::actualizar_parametros_de_celula_y_muerte_con_o2_y_oncoproteina(Fenotipo& fenotipo, double dt){

	if(fenotipo.muerte.muerta == true || tipo != 0){
		return;
	}
	static bool indices_iniciados = false;
	static int indice_fase_inicial;
	static int indice_fase_final;
	static int indice_necrosis;

	static int indice_del_oxigeno = microambiente->encontrar_indice_de_densidad("oxigeno");
	if(indices_iniciados == false){

		if(fenotipo.ciclo.pCiclo_Modelo->codigo == Constantes::ciclo_Ki67){
			indice_fase_inicial = fenotipo.ciclo.pCiclo_Modelo->encontrar_indice_de_la_fase(Constantes::Ki67_negativa);
			indice_necrosis = fenotipo.muerte.encontrar_indice_del_ciclo_de_muerte(Constantes::ciclo_de_muerte_necrosis);
			indice_fase_final = fenotipo.ciclo.pCiclo_Modelo->encontrar_indice_de_la_fase(Constantes::Ki67_positiva_premitotica);
			indices_iniciados = true;

		}

		if( fenotipo.ciclo.pCiclo_Modelo->codigo == Constantes::ciclo_vida ){

			indice_fase_inicial = fenotipo.ciclo.pCiclo_Modelo->encontrar_indice_de_la_fase(Constantes::viva);
			indice_necrosis = fenotipo.muerte.encontrar_indice_del_ciclo_de_muerte(Constantes::ciclo_de_muerte_necrosis);
			indice_fase_final = fenotipo.ciclo.pCiclo_Modelo->encontrar_indice_de_la_fase(Constantes::viva);
			indices_iniciados = true;
		}
	}
	double pO2 = (vector_de_densidades_mas_cercano())[indice_del_oxigeno];
	// pO2 = 30.;//Debug
	
	double multiplicador = 1.0;
	if(pO2 < parametros.o2_saturacion_para_la_proliferacion){
		multiplicador = (pO2 - parametros.o2_limite_de_proliferacion) /
		(parametros.o2_saturacion_para_la_proliferacion - parametros.o2_limite_de_proliferacion);
	}

	if(pO2 < parametros.o2_limite_de_proliferacion){
		multiplicador= 0.0;
	}

    fenotipo.ciclo.actualizar_mis_tasas_de_transicion(indice_fase_inicial, indice_fase_final) = multiplicador *fenotipo.ciclo.tasa_aleatoria;

	multiplicador = 0.0;
	if(pO2 < parametros.o2_necrosis_limite){
		multiplicador = (parametros.o2_necrosis_limite - pO2) /
		(parametros.o2_necrosis_limite - parametros.o2_necrosis_max);
	}
	if(pO2 < parametros.o2_necrosis_max){
		multiplicador = 1.0;
	}
	fenotipo.muerte.tasas[indice_necrosis] = multiplicador * parametros.tasa_necrosis_max;




      fenotipo.ciclo.actualizar_mis_tasas_de_transicion(indice_fase_inicial, indice_fase_final) *= fenotipo.secrecion.oncoproteina;
	//Debug
	//output info only for cancer cells
	if (tipo==0){
		int half_day = pg->tiempo_total / (60*12);
		std::ofstream file("out/simulation_data"+std::to_string(half_day)+".csv", std::ios::app);
		if (file.is_open()) {
		file  << pO2 << "," 
			<< fenotipo.secrecion.oncoproteina  << "," << fenotipo.ciclo.tasa_aleatoria <<"," << fenotipo.ciclo.actualizar_mis_tasas_de_transicion(indice_fase_inicial, indice_fase_final) << "," << multiplicador * parametros.tasa_necrosis_max << "\n";
		}
	}
	//End Debug
      

}

/**
 * @brief Advances cell phenotype functions based on standard model
 * @details This method updates the cell's phenotype over time, handling:
 * - Cell cycle progression
 * - Volume changes
 * - Death processes
 * - Division flagging
 * - Cell removal flagging
 * 
 * It uses oxygen-dependent parameters but doesn't directly check oxygen levels.
 * In cancer research, this represents the basic mechanisms of tumor growth and cell turnover.
 * 
 * @param hora_global Current simulation time
 * @param dt_ciclo Time step for cycle advancement
 */
// void Celula::avanzar_funciones_del_fenotipo(double hora_global, double dt_ciclo){


// 		actualizar_parametros_de_celula_y_muerte_con_o2(fenotipo, dt_ciclo);

// 		if(fenotipo.muerte.chequear_muerte(dt_ciclo) == true){
// 			fenotipo.ciclo.sync_con_ciclo_modelo(fenotipo.muerte.ciclo_actual());

// 			fenotipo.secrecion.set_todas_las_secreciones_a_cero();
// 			fenotipo.secrecion.multiplicar_los_consumos_por_un_factor(0.10);

// 			if(fenotipo.ciclo.fase_actual().funcion_de_entrada){

// 				fenotipo.ciclo.fase_actual().funcion_de_entrada(fenotipo.volumen, fenotipo.muerte.parametros_actuales());

// 			}

// 		}

// 		fenotipo.ciclo.avanzar_en_el_ciclo( fenotipo.volumen, tiempo_desde_el_ultimo_ciclo, fenotipo.ciclo.tasas_de_transicion, fenotipo.muerte.parametros_actuales() );


// 		if(fenotipo.ciclo.actualizar_volumen()){

// 			actualizar_volumen(fenotipo, tiempo_desde_el_ultimo_ciclo);
// 			fenotipo.volumen.cambio_el_volumen=true;
// 		}

// 		if(fenotipo.ciclo.flagged_para_dividirse){

//             if(fenotipo.muerte.muerta){
//                     std::cout<< "Se divide celula muerta id: " << id << "\n";
//                     std::cin.get();
//             }

// 			celulas_listas_para_dividirse.push_back(this);
// 			fenotipo.ciclo.flagged_para_dividirse=false;
// 		}

// 		if(fenotipo.ciclo.flagged_para_remover){
// 			celulas_listas_para_remover.push_back(this);
// 			fenotipo.ciclo.flagged_para_remover=false;
// 		}


// 	return;
// }

/**
 * @brief Advances cell phenotype functions with explicit oxygen dependency
 * @details Enhanced version of the basic phenotype advancement that explicitly checks
 * oxygen levels to determine if the cell can proliferate. This models the common
 * cancer phenotype where cells stop dividing in hypoxic regions but may still survive.
 * 
 * Key cancer research aspects:
 * - Models oxygen-dependent cell cycle progression
 * - Implements different behavior for dead vs. living cells
 * - Captures cellular adaptation to hypoxic tumor microenvironments
 * 
 * @param hora_global Current simulation time
 * @param dt_ciclo Time step for cycle advancement
 */
// void Celula::avanzar_funciones_del_fenotipo_con_O2(double hora_global, double dt_ciclo){

//         static int indice_del_oxigeno = microambiente->encontrar_indice_de_densidad("oxigeno");
//         double pO2 = (vector_de_densidades_mas_cercano())[indice_del_oxigeno];
// 		actualizar_parametros_de_celula_y_muerte_con_o2(fenotipo, dt_ciclo);
// 		if(fenotipo.muerte.chequear_muerte(dt_ciclo) == true){
// 			fenotipo.ciclo.sync_con_ciclo_modelo(fenotipo.muerte.ciclo_actual());
// 			fenotipo.secrecion.set_todas_las_secreciones_a_cero();
// 			fenotipo.secrecion.multiplicar_los_consumos_por_un_factor(0.10);

// 			if(fenotipo.ciclo.fase_actual().funcion_de_entrada){

// 				fenotipo.ciclo.fase_actual().funcion_de_entrada(fenotipo.volumen, fenotipo.muerte.parametros_actuales());

// 			}

// 		}
//         if(!fenotipo.muerte.muerta && pO2 > parametros.o2_limite_de_proliferacion){

//             fenotipo.ciclo.avanzar_en_el_ciclo( fenotipo.volumen, tiempo_desde_el_ultimo_ciclo, fenotipo.ciclo.tasas_de_transicion, fenotipo.muerte.parametros_actuales() );

//             if(fenotipo.ciclo.actualizar_volumen()){

//                 actualizar_volumen(fenotipo, tiempo_desde_el_ultimo_ciclo);
//                 fenotipo.volumen.cambio_el_volumen=true;
//             }

//             if(fenotipo.ciclo.flagged_para_dividirse){

//                 celulas_listas_para_dividirse.push_back(this);
//                 fenotipo.ciclo.flagged_para_dividirse=false;

//             }

//             if(fenotipo.ciclo.flagged_para_remover){
//                 celulas_listas_para_remover.push_back(this);
//                 fenotipo.ciclo.flagged_para_remover=false;
//             }

//         }else if(fenotipo.muerte.muerta){

//             fenotipo.ciclo.avanzar_en_el_ciclo( fenotipo.volumen, tiempo_desde_el_ultimo_ciclo, fenotipo.ciclo.tasas_de_transicion, fenotipo.muerte.parametros_actuales() );

//             if(fenotipo.ciclo.actualizar_volumen()){

//                 actualizar_volumen(fenotipo, tiempo_desde_el_ultimo_ciclo);
//                 fenotipo.volumen.cambio_el_volumen=true;
//             }

//             if(fenotipo.ciclo.flagged_para_remover){
//                 celulas_listas_para_remover.push_back(this);
//                 fenotipo.ciclo.flagged_para_remover=false;
//             }

//         }
// 	return;
// }

/**
 * @brief Advances cell phenotype functions considering oxygen and oncoprotein
 * @details Most sophisticated version of phenotype advancement that takes into account:
 * - Oxygen availability
 * - Oncoprotein expression
 * - Cell attachment status
 * - Cell type specific behaviors (especially for lymphocytes)
 * 
 * This method is particularly important for cancer research because it integrates:
 * - Molecular factors (oncoprotein)
 * - Microenvironmental factors (oxygen)
 * - Cellular interactions (attachment)
 * - Cell-specific behaviors that may vary between tumor and immune cells
 * 
 * @param hora_global Current simulation time
 * @param dt_ciclo Time step for cycle advancement
 */
void Celula::avanzar_funciones_del_fenotipo_con_O2_y_oncoproteina(double hora_global, double dt_ciclo){


		// if(fenotipo.ciclo.pCiclo_Modelo->nombre != "Vida"){
		// bool debug= fenotipo.ciclo.actualizar_volumen() ;
		// auto st="No";
		// if(debug){
		// 	st="Si";
		// }
		// std::cout << "#" << fenotipo.ciclo.pCiclo_Modelo->nombre << " " << fenotipo.volumen.citoplasma_tasa_de_cambio<< st<< tipo<< "#" << std::endl;
			
		// }

        if(tipo!=2 && adherida==true){
            return;
        }


        static int indice_del_oxigeno = microambiente->encontrar_indice_de_densidad("oxigeno");
        double pO2 = (vector_de_densidades_mas_cercano())[indice_del_oxigeno];
		// pO2 = 30.;//Debug
		
		// Print detallado para encontrar IDs altos
		// static int contador_print = 0;
		// contador_print++;
		// if (contador_print % 1000 == 0) {  // Print cada 1000 llamadas para no saturar
		// 	std::cout << "CELULA ID: " << id << " tipo: " << tipo << " pO2: " << pO2 << " adherida: " << adherida << std::endl;
		// }
		
		// // Print específico para IDs altos
		// if (id > 3963) {
		// 	std::cout << "*** ID ALTO ENCONTRADO: " << id << " tipo: " << tipo << " pO2: " << pO2 << " adherida: " << adherida << std::endl;
		// }
		
		// // Print específico para células de tipo 0 con IDs altos (divisiones)
		// if (tipo == 0 && id > 3963) {
		// 	std::cout << "*** CELULA TIPO 0 CON ID ALTO (DIVISION): " << id << " tipo: " << tipo << " pO2: " << pO2 << " madre: " << madre << std::endl;
		// }
		
		// if (tipo == 0) {
		// 	std::cout << "Nivel de oxígeno (pO2) celula id " << id << ": " << pO2 << std::endl;
		// }
        actualizar_parametros_de_celula_y_muerte_con_o2_y_oncoproteina(fenotipo, dt_ciclo);

		if(fenotipo.muerte.chequear_muerte(dt_ciclo) == true){
			fenotipo.ciclo.sync_con_ciclo_modelo(fenotipo.muerte.ciclo_actual());
            fenotipo.ciclo.indice_de_la_fase_actual = 0;
            fenotipo.ciclo.tiempo_acumulado_en_la_fase = 0.0;
			fenotipo.secrecion.set_todas_las_secreciones_a_cero();
			fenotipo.secrecion.multiplicar_los_consumos_por_un_factor(0.10);
			if(tipo==2){
                es_movil(false);
                if(adherida==true){
                    adherida=false;
                    celula_adherida->adherida=false;
                    celula_adherida->celula_adherida=NULL;
                    celula_adherida=NULL;
                }
            }

			if(fenotipo.ciclo.fase_actual().funcion_de_entrada){
				fenotipo.ciclo.fase_actual().funcion_de_entrada(fenotipo.volumen, fenotipo.muerte.parametros_actuales());
			}
		}

		if(tipo==2 && fenotipo.muerte.muerta == false){
            fenotipo.muerte.tasas[0] = 1/abs(((1/fenotipo.muerte.tasas[0]) - (dt_ciclo*dt_ciclo)));
		}

        if(fenotipo.muerte.muerta == false && pO2 > parametros.o2_limite_de_proliferacion && adherida==false){

            fenotipo.ciclo.avanzar_en_el_ciclo( fenotipo.volumen, tiempo_desde_el_ultimo_ciclo, fenotipo.ciclo.tasas_de_transicion, fenotipo.muerte.parametros_actuales() );


            if(fenotipo.ciclo.actualizar_volumen()){

                actualizar_volumen(fenotipo, tiempo_desde_el_ultimo_ciclo);
                fenotipo.volumen.cambio_el_volumen=true;
            }

            if(fenotipo.ciclo.flagged_para_dividirse){

                celulas_listas_para_dividirse.push_back(this);
                fenotipo.ciclo.flagged_para_dividirse=false;
            }

            if(fenotipo.ciclo.flagged_para_remover){
                celulas_listas_para_remover.push_back(this);
                fenotipo.ciclo.flagged_para_remover=false;
            }

        }else if(fenotipo.muerte.muerta == true){

			// // if(fenotipo.ciclo.pCiclo_Modelo->nombre != "Vida"){
			// bool debug= fenotipo.ciclo.actualizar_volumen() ;
			// auto st="No";
			// if(debug){
			// 	st="Si";
			// }
			// std::cout << "#" << fenotipo.ciclo.pCiclo_Modelo->nombre << " " << fenotipo.volumen.citoplasma_tasa_de_cambio<< st<< tipo<< "#" << std::endl;
				
			// // }

            fenotipo.ciclo.avanzar_en_el_ciclo( fenotipo.volumen, tiempo_desde_el_ultimo_ciclo, fenotipo.ciclo.tasas_de_transicion, fenotipo.muerte.parametros_actuales() );

            if(fenotipo.ciclo.actualizar_volumen()){

                actualizar_volumen(fenotipo, tiempo_desde_el_ultimo_ciclo);
                fenotipo.volumen.cambio_el_volumen=true;
            }

            if(fenotipo.ciclo.flagged_para_remover){
                celulas_listas_para_remover.push_back(this);
                fenotipo.ciclo.flagged_para_remover=false;
            }

        }

	// Estadísticas de IDs cada cierto tiempo
	static int contador_estadisticas = 0;
	static int max_id_visto = 0;
	contador_estadisticas++;
	
	if (id > max_id_visto) {
		max_id_visto = id;
	}
	
	//Debug
	// if (contador_estadisticas % 5000 == 0) {
	// 	std::cout << "=== ESTADISTICAS IDs === Max ID visto: " << max_id_visto << " (llamadas: " << contador_estadisticas << ")" << std::endl;
	// }

	return;
}



/**
 * @brief Creates a new cell and registers it in the simulation
 * @details Factory function that creates a cell, registers it in the global cell collection
 * and assigns a unique ID. This is essential for cancer simulations as it manages
 * the creation of new cells during proliferation and initialization.
 * 
 * @return Pointer to the newly created cell
 */
Celula* crear_celula(void){

    Celula* pNuevo;
	pNuevo = new Celula;
	todas_las_celulas.push_back(pNuevo);
	celulas_para_registrar_en_voxeles.push_back(pNuevo);
	pNuevo->indice=todas_las_celulas.size()-1;
	pg->numero_id += 1;
	pNuevo->id=pg->numero_id;


	return pNuevo;
}

/**
 * @brief Creates a new lymphocyte and registers it in the simulation
 * @details Similar to crear_celula but creates a lymphocyte (immune cell) specifically.
 * This is important for cancer immunology research as it enables modeling of
 * immune cell interactions with tumor cells.
 * 
 * @return Pointer to the newly created lymphocyte
 */
Linfocito* crear_linfocito(void){
    Linfocito* pNuevo;
	pNuevo = new Linfocito;
	todas_las_celulas.push_back(pNuevo);
	celulas_para_registrar_en_voxeles.push_back(pNuevo);
	pNuevo->indice=todas_las_celulas.size()-1;
	pg->numero_id += 1;
	pNuevo->id=pg->numero_id;

	return pNuevo;
}

/**
 * @brief Sets the cell position in 3D space with boundary handling
 * @details Positions the cell at specified coordinates, handling periodic boundary
 * conditions if enabled. This is critical for spatial modeling of tumors and their
 * microenvironment in cancer research.
 * 
 * @param x X-coordinate
 * @param y Y-coordinate
 * @param z Z-coordinate
 * @return True if position was successfully set
 */
bool Celula::set_posicion(double x, double y, double z){

	if (pg->condiciones_de_periodicidad) {
		if (x < pg->rango_en_X[0] && pg->condiciones_de_periodicidad_x )  {
			x = pg->rango_en_X[1] - (pg->rango_en_X[0] - x);
		} else if (x > pg->rango_en_X[1] && pg->condiciones_de_periodicidad_x ) {
			x = pg->rango_en_X[0] + (x - pg->rango_en_X[1]);
		}
		if (y < pg->rango_en_Y[0] && pg->condiciones_de_periodicidad_y ) {
			y = pg->rango_en_Y[1] - (pg->rango_en_Y[0] - y);
		} else if (y > pg->rango_en_Y[1] && pg->condiciones_de_periodicidad_y ) {
			y = pg->rango_en_Y[0] + (y - pg->rango_en_Y[1]);
		}
		if (z < pg->rango_en_Z[0] && pg->condiciones_de_periodicidad_z ) {
			z = pg->rango_en_Z[1] - (pg->rango_en_Z[0] - z);
		} else if (z > pg->rango_en_Z[1] && pg->condiciones_de_periodicidad_z ) {
			z = pg->rango_en_Z[0] + (z - pg->rango_en_Z[1]);
		}
	}
	posicion.x=x;
	posicion.y=y;
	posicion.z=z;

	return true;
}

bool Celula::set_posicion(Vector posicion){

    return set_posicion(posicion.x, posicion.y, posicion.z);
}


Celula* Celula::dividir(){

	Celula* hija = crear_celula();
	hija->tipo = tipo;
	hija->madre = id;
	hija->nombre = nombre;
	hija->hora_del_ultimo_ciclo = hora_del_ultimo_ciclo;
	hija->hora_de_la_ultima_mecanica = hora_de_la_ultima_mecanica;
	hija->parametros = parametros;

	// Debug: Print cuando se crea una célula por división
	// std::cout << "*** DIVISION CELULAR: Madre ID " << id << " (tipo " << tipo << ") -> Hija ID " << hija->id << " (tipo " << hija->tipo << ")" << std::endl;


	double angulo_temporal = 6.28318530717959*rng->RandomNumber();
	double phi_temporal = 3.1415926535897932384626433832795*rng->RandomNumber();

	double radio = fenotipo.geometria.radio;
	Vector vector_random (0.0, 0.0, 0.0);

    if(pg->crecer_al_costado){
    	if(rng->RandomNumber()<0.05){
            vector_random.x= cos( angulo_temporal ) * sin( phi_temporal );
            vector_random.y= sin( angulo_temporal ) * sin( phi_temporal );
            vector_random.z= cos( phi_temporal );
            }else{
                vector_random.x= cos( angulo_temporal );
                vector_random.y= sin( angulo_temporal );
                vector_random.z= 0;
                }
        }else{
            vector_random.x= cos( angulo_temporal ) * sin( phi_temporal );
            vector_random.y= sin( angulo_temporal ) * sin( phi_temporal );
            vector_random.z= cos( phi_temporal );
        }


	//Le asigno una posición a la hija
	hija->set_posicion(posicion.x + 0.206299474 * radio*vector_random.x,
						 posicion.y + 0.206299474 * radio*vector_random.y,
						 posicion.z + 0.206299474 * radio*vector_random.z);
	hija->actualizar_voxel_del_microambiente();

	//Modifico la posición de la madre
	set_posicion(posicion.x - 0.206299474*radio*vector_random.x,
		posicion.y - 0.206299474*radio*vector_random.y,
		posicion.z - 0.206299474*radio*vector_random.z);
	actualizar_voxel_del_microambiente();

	fenotipo.volumen.dividir();
	fenotipo.geometria.actualizar(fenotipo.volumen);
	hija->fenotipo = fenotipo;
    hija->adherida = false;
    
    if(tipo==0){
        hija->fenotipo.ciclo.tasa_aleatoria=1/(rng->NormalRandom_CM(38.6, 3.7)*60);
        hija->fenotipo.ciclo.actualizar_mis_tasas_de_transicion(0,0) = fenotipo.ciclo.tasa_aleatoria;
    }    


	return hija;
}



/**
 * @brief Removes a cell from the simulation
 * @details Deletes a cell and updates the global cell collection. This is used when
 * cells are removed from the simulation due to completion of death processes.
 * In cancer research, this models cell clearance after cell death, which affects 
 * tumor dynamics and treatment responses.
 * 
 * @param indice Index of the cell to remove in the todas_las_celulas array
 */
void Celula::morir(int indice){
	delete todas_las_celulas[indice];

	todas_las_celulas[ todas_las_celulas.size()-1 ]->indice=indice;

	todas_las_celulas[indice] = todas_las_celulas[ todas_las_celulas.size()-1 ];

	todas_las_celulas.pop_back();

	return;


}


/**
 * @brief Calculates mechanical interaction forces between two cells
 * @details Implements the core mechanical interactions between cells, including:
 * - Repulsion when cells overlap (preventing excessive compression)
 * - Adhesion when cells are within adhesion distance
 * 
 * These forces are critical for cancer modeling as they determine:
 * - Tumor structure and packing density
 * - Mechanical stress in the tumor microenvironment
 */
void Celula::agregar_potenciales(Celula* otra_celula){

	if( this->id == otra_celula->id ){
		return;
	}


	desplazamiento.x = posicion.x - otra_celula->posicion.x;
	desplazamiento.y = posicion.y - otra_celula->posicion.y;
	desplazamiento.z = posicion.z - otra_celula->posicion.z;

    //Agrego esto por las Condiciones Periódicas de Contorno (funciona también sin CPC)
    desplazamiento.x = desplazamiento.x - (2*pg->rango_en_X[1])*round(desplazamiento.x/(2*pg->rango_en_X[1]));
	desplazamiento.y = desplazamiento.y - (2*pg->rango_en_Y[1])*round(desplazamiento.y/(2*pg->rango_en_Y[1]));
	desplazamiento.z = desplazamiento.z - (2*pg->rango_en_Z[1])*round(desplazamiento.z/(2*pg->rango_en_Z[1]));




	double distancia = desplazamiento.x * desplazamiento.x +
	desplazamiento.y * desplazamiento.y + desplazamiento.z * desplazamiento.z;


	distancia = std::max(sqrt(distancia), 0.00001);

	double R = fenotipo.geometria.radio + otra_celula->fenotipo.geometria.radio;

	double temp_r;
	if( distancia > R )
	{
		temp_r=0;
	}
	else
	{
		temp_r = -distancia; // -d
		temp_r /= R; // -d/R
		temp_r += 1.0; // 1-d/R
		temp_r *= temp_r; // (1-d/R)^2

	}

	if(this->tipo == otra_celula->tipo){
	double repulsion_efectiva = sqrt( fenotipo.mecanica.fuerza_de_repulsion_cc * otra_celula->fenotipo.mecanica.fuerza_de_repulsion_cc );
	temp_r *= repulsion_efectiva;
    }else{
	double repulsion_efectiva = sqrt( fenotipo.mecanica.fuerza_de_repulsion_co * otra_celula->fenotipo.mecanica.fuerza_de_repulsion_co );
	temp_r *= repulsion_efectiva;
    }


	double distancia_de_interaccion_max = fenotipo.mecanica.distancia_de_adhesion_maxima_relativa * fenotipo.geometria.radio +
		otra_celula->fenotipo.mecanica.distancia_de_adhesion_maxima_relativa * otra_celula->fenotipo.geometria.radio;

	if(distancia < distancia_de_interaccion_max )
	{
		double temp_a = -distancia; // -d
		temp_a /= distancia_de_interaccion_max; // -d/S
		temp_a += 1.0; // 1 - d/S
		temp_a *= temp_a; // (1-d/S)^2

		if(this->tipo == otra_celula->tipo){
		double adhesion_efectiva = sqrt( fenotipo.mecanica.fuerza_de_adhesion_cc * otra_celula->fenotipo.mecanica.fuerza_de_adhesion_cc );
		temp_a *= adhesion_efectiva;
		}else{
		double adhesion_efectiva = sqrt( fenotipo.mecanica.fuerza_de_adhesion_co * otra_celula->fenotipo.mecanica.fuerza_de_adhesion_co );
		temp_a *= adhesion_efectiva;
		}

		temp_r -= temp_a;
	}
	if( fabs(temp_r) < 1e-16 )
	{ return; }
	temp_r /= distancia;
	
    //axpy(&velocidad, temp_r, desplazamiento);
    velocidad.x += temp_r * desplazamiento.x;
	velocidad.y += temp_r * desplazamiento.y;
	velocidad.z += temp_r * desplazamiento.z;
    
    //Le agrego el potencial a la otra célula también
    otra_celula->velocidad.x -= temp_r * desplazamiento.x;
	otra_celula->velocidad.y -= temp_r * desplazamiento.y;
	otra_celula->velocidad.z -= temp_r * desplazamiento.z; 

	return;
}

/**
 * @brief Calculates and applies bottom boundary mechanical interactions
 * @details Computes the mechanical forces between the cell and the bottom boundary
 * of the simulation domain. This includes both adhesive and repulsive forces based 
 * on distance, which are then applied to the cell's velocity.
 * 
 * In cancer research, this method is significant for modeling:
 * - Tumor-basement membrane interactions
 * - Cancer cell adhesion to boundaries (critical in invasion modeling)
 * - Mechanical constraints on tumor growth
 * - Boundary effects on tumor morphology
 * 
 * The method calculates adhesion and repulsion forces separately and combines them,
 * with forces becoming stronger at closer distances following a quadratic relationship.
 */
void Celula::agregar_potenciales_mb(){

    if( pg->interactuar_con_mb==false ){
		return;
	}


	desplazamiento.x = posicion.x - posicion.x;
	desplazamiento.y = posicion.y - posicion.y;
	desplazamiento.z = posicion.z - pg->rango_en_Z[0];

	double distancia = desplazamiento.x * desplazamiento.x +
	desplazamiento.y * desplazamiento.y + desplazamiento.z * desplazamiento.z;


	distancia = std::max(sqrt(distancia), 0.00001);


	double distancia_de_interaccion_max = fenotipo.mecanica.distancia_de_adhesion_maxima_relativa * fenotipo.geometria.radio;

	double temp_a=0;
	if(distancia < distancia_de_interaccion_max)
	{
		temp_a= (1- distancia/distancia_de_interaccion_max);
		temp_a*=temp_a;
		temp_a*=-fenotipo.mecanica.fuerza_de_adhesion_mb;
	}
	double temp_r = 0;
	if(distancia < fenotipo.geometria.radio)
	{
		temp_r = (1- distancia/fenotipo.geometria.radio);
		temp_r *= temp_r;
		temp_r *= fenotipo.mecanica.fuerza_de_repulsion_mb;
	}
	temp_r += temp_a;
	if( fabs( temp_r ) < 1e-16 )
	{ return; }

	//axpy(&velocidad, temp_r, desplazamiento);
    velocidad.x += temp_r * desplazamiento.x;
	velocidad.y += temp_r * desplazamiento.y;
	velocidad.z += temp_r * desplazamiento.z;
	return;


}


/**
 * @brief Calculates and applies mechanical interactions with all domain boundaries
 * @details Computes the mechanical forces between the cell and all six boundaries 
 * of the simulation domain (±X, ±Y, ±Z). For each non-periodic boundary, this method
 * calculates both adhesive and repulsive forces based on distance, then applies
 * these forces to the cell's velocity.
 * 
 * In cancer research, this comprehensive boundary interaction model is crucial for:
 * - Studying tumor growth constraints in confined environments
 * - Modeling tumor-basement membrane interactions in all directions
 * - Simulating cancer cell invasion through tissue boundaries
 * - Investigating the effects of spatial constraints on tumor morphology
 * - Replicating in vitro experimental setups with defined boundaries
 * 
 * The method checks each boundary separately and only applies forces for non-periodic
 * boundaries when the cell is within interaction distance, collecting all relevant
 * displacement vectors and applying forces to each.
 */
void Celula::agregar_potenciales_mb_2(){

    if( pg->condiciones_de_periodicidad_x==true && pg->condiciones_de_periodicidad_y==true && pg->condiciones_de_periodicidad_z==true ){
		return;
	}
    std::cout << "Entre a calcular membrana \n";
	double distancia_de_interaccion_max = fenotipo.mecanica.distancia_de_adhesion_maxima_relativa * fenotipo.geometria.radio;

	double distancia;
	std::vector<double> distancias;
	std::vector<Vector> desplazamientos;

	//Comparación con x_0
	desplazamiento.x = posicion.x - pg->rango_en_X[0];
	desplazamiento.y = 0;
	desplazamiento.z = 0;

	distancia = desplazamiento.x * desplazamiento.x;
	distancia = std::max(sqrt(distancia), 0.00001);

	if(distancia < distancia_de_interaccion_max && !pg->condiciones_de_periodicidad_x){
	distancias.push_back(distancia);
	desplazamientos.push_back(desplazamiento);

	}

    //Comparación con y_0
	desplazamiento.x = 0;
	desplazamiento.y = posicion.y - pg->rango_en_Y[0];
	desplazamiento.z = 0;

	distancia = desplazamiento.y * desplazamiento.y;
	distancia = std::max(sqrt(distancia), 0.00001);

	if(distancia < distancia_de_interaccion_max && !pg->condiciones_de_periodicidad_y){
	distancias.push_back(distancia);
	desplazamientos.push_back(desplazamiento);

	}

    //Comparación con z_0
	desplazamiento.x = 0;
	desplazamiento.y = 0;
	desplazamiento.z = posicion.z - pg->rango_en_Z[0];

	distancia = desplazamiento.z * desplazamiento.z;
	distancia = std::max(sqrt(distancia), 0.00001);

	if(distancia < distancia_de_interaccion_max && !pg->condiciones_de_periodicidad_z){
	distancias.push_back(distancia);
	desplazamientos.push_back(desplazamiento);

	}

    //Comparación con x_Max
	desplazamiento.x = pg->rango_en_X[1]- posicion.x;
	desplazamiento.y = 0;
	desplazamiento.z = 0;

	distancia = desplazamiento.x * desplazamiento.x;
	distancia = std::max(sqrt(distancia), 0.00001);

	if(distancia < distancia_de_interaccion_max && !pg->condiciones_de_periodicidad_x){
	distancias.push_back(distancia);
	desplazamientos.push_back(desplazamiento);

	}

    //Comparación con y_Max
	desplazamiento.x = 0;
	desplazamiento.y = pg->rango_en_Y[1]- posicion.y;
	desplazamiento.z = 0;

	distancia = desplazamiento.y * desplazamiento.y;
	distancia = std::max(sqrt(distancia), 0.00001);

	if(distancia < distancia_de_interaccion_max && !pg->condiciones_de_periodicidad_y){
	distancias.push_back(distancia);
	desplazamientos.push_back(desplazamiento);

	}

    //Comparación con z_Max
	desplazamiento.x = 0;
	desplazamiento.y = 0;
	desplazamiento.z = pg->rango_en_Z[1]- posicion.z;

	distancia = desplazamiento.z * desplazamiento.z;
	distancia = std::max(sqrt(distancia), 0.00001);

	if(distancia < distancia_de_interaccion_max && !pg->condiciones_de_periodicidad_z){
	distancias.push_back(distancia);
	desplazamientos.push_back(desplazamiento);

	}

	if(distancias.size()>0){
        for(unsigned int i=0; i < distancias.size(); i++){

            double temp_a=0;
            // Adhesion a la membrana basal
            temp_a= (1- distancias[i]/distancia_de_interaccion_max);
            temp_a*=temp_a;
            temp_a*=-fenotipo.mecanica.fuerza_de_adhesion_mb;

            // Repulsion de la membrana basal
            double temp_r = 0;
            if(distancias[i] < fenotipo.geometria.radio){
                temp_r = (1- distancias[i]/fenotipo.geometria.radio);
                temp_r *= temp_r;
                temp_r *= fenotipo.mecanica.fuerza_de_repulsion_mb;
            }
            temp_r += temp_a;
            if( fabs( temp_r ) > 1e-16 )
            { axpy(&velocidad, temp_r, desplazamientos[i]); }


        }

	}

	return;


}


/**
 * @brief Updates the cell's position based on accumulated velocity
 * @details Applies the calculated velocity to update the cell's position using a
 * modified Euler integration scheme. This method also handles periodic boundary
 * conditions and updates the cell's voxel assignment.
 * 
 * In cancer research, this method is critical for modeling:
 * - Tumor growth patterns and spatial organization
 * - Cell migration and invasion
 * - Mechanical interactions within the tumor microenvironment
 * - Spatial heterogeneity in tumor tissue
 * 
 * @param dt Time step size for the update
 */
void Celula::actualizar_posicion( double dt ){

	static double d1;
	static double d2;
	static bool constants_defined = false;
	if( constants_defined == false )
	{
		d1 = dt;
		d1 *= 1.5;
		d2 = dt;
		d2 *= -0.5;
		constants_defined = true;
	}


    posicion.x += d1 * velocidad.x;
	posicion.y += d1 * velocidad.y;
	posicion.z += d1 * velocidad.z;
	
    posicion.x += d2 * velocidad_anterior.x;
	posicion.y += d2 * velocidad_anterior.y;
	posicion.z += d2 * velocidad_anterior.z;    

	velocidad_anterior = velocidad;

	velocidad.x = 0;
	velocidad.y = 0;
	velocidad.z = 0;

	set_posicion(posicion.x, posicion.y, posicion.z);
	actualizar_voxel_del_microambiente();

	return;

}


/**
 * @brief Initiates the cell death process
 * @details Transitions the cell to a death pathway by synchronizing with the appropriate
 * death cycle model, resetting the cell's phase, and adjusting its secretion rates.
 * For lymphocytes, also handles detachment from target cells.
 * 
 * In cancer research, this method is critical for modeling:
 * - Cell death in response to therapy
 * - Programmed cell death (apoptosis) mechanisms
 * - Stress-induced cell death (necrosis)
 * - Cell clearance from the tumor microenvironment
 * - Immune cell detachment after target cell engagement
 * 
 * @param indice_del_ciclo_de_muerte Index of the death cycle to initiate
 */
void Celula::comenzar_muerte( int indice_del_ciclo_de_muerte )
{
    fenotipo.muerte.muerta = true;
	fenotipo.muerte.comenzar_muerte(indice_del_ciclo_de_muerte);
	fenotipo.ciclo.sync_con_ciclo_modelo(fenotipo.muerte.ciclo_actual());
	fenotipo.ciclo.indice_de_la_fase_actual = 0;
	fenotipo.ciclo.tiempo_acumulado_en_la_fase = 0.0;
	fenotipo.secrecion.set_todas_las_secreciones_a_cero();
	fenotipo.secrecion.multiplicar_los_consumos_por_un_factor(0.10);
    if(tipo==2){
        es_movil(false);
        if(adherida==true){
            adherida=false;
            celula_adherida->adherida=false;
            celula_adherida->celula_adherida=NULL;
            celula_adherida=NULL;
        }
	}

	if( fenotipo.ciclo.fase_actual().funcion_de_entrada)
	{
		fenotipo.ciclo.fase_actual().funcion_de_entrada(fenotipo.volumen, fenotipo.muerte.parametros_actuales());
	}

	return;
}


/* 
 * ===============================================================================================
 * CELL-MICROENVIRONMENT INTERACTIONS
 * ===============================================================================================
 *
 * The following section contains methods for cell interactions with the biochemical microenvironment.
 * These methods are critical for modeling how cancer cells sense and modify their surroundings,
 * including oxygen consumption, secretion of factors, and responses to environmental conditions.
 *
 * In cancer research, these interactions help model:
 * - Metabolic adaptations (e.g., Warburg effect)
 * - Hypoxic regions within tumors
 * - Angiogenic signaling 
 * - Extracellular matrix modifications
 * - Chemotactic migration
 */

/*MICROAMBIENTE Y OXÍGENO*/

/**
 * @brief Registers the cell with the microenvironment
 * @details Sets up the cell's secretion and consumption vectors to match the
 * microenvironment's substrates. This is essential for the cell to interact
 * with the biochemical environment.
 * 
 * In cancer research, this enables modeling:
 * - Cell metabolism (nutrients uptake)
 * - Signaling factor secretion
 * - Hypoxic responses
 * - Cell-cell communication via diffusible factors
 * 
 * @param microambiente Pointer to the microenvironment (unused parameter - uses class member)
 */
void Celula::registrar_microambiente( Microambiente* ){

	fenotipo.secrecion.tasas_de_secrecion.resize( microambiente->vector_de_densidades(0).size() , 0.0 );
	fenotipo.secrecion.densidades_de_saturacion.resize( microambiente->vector_de_densidades(0).size() , 0.0 );
	fenotipo.secrecion.tasas_de_consumo.resize( microambiente->vector_de_densidades(0).size() , 0.0 );
	fenotipo.secrecion.tasas_de_exportacion_neta.resize( microambiente->vector_de_densidades(0).size() , 0.0 );

	temp_celula_fuente_sumidero_solver1.resize( microambiente->vector_de_densidades(0).size() , 0.0 );
	temp_celula_fuente_sumidero_solver2.resize( microambiente->vector_de_densidades(0).size() , 1.0 );

	temp_celula_fuente_sumidero_exportacion_solver1.resize( microambiente->vector_de_densidades(0).size() , 0.0 );
	temp_celula_fuente_sumidero_exportacion_solver1.resize( microambiente->vector_de_densidades(0).size() , 0.0 );


	return;

}

/**
 * @brief Updates the cell's voxel index in the microenvironment
 * @details Maps the cell to the appropriate voxel in the microenvironment grid based
 * on its current position. This is essential for spatial interactions with the
 * biochemical environment.
 * 
 * In cancer research, this enables modeling:
 * - Local variations in oxygen and nutrients within tumors
 * - Spatial gradients of signaling molecules
 * - Regional differences in treatment efficacy
 */
void Celula::actualizar_voxel_del_microambiente(){

	if(!microambiente->mgrilla.es_valida_la_posicion(posicion.x, posicion.y, posicion.z)){
    std::cout << "Posicion invalida: " << posicion.x << " " << posicion.y << " " << posicion.z << " " << id << "\n ";
		voxel_del_microambiente =-1;
		es_activa =false;
		return;
	}

	voxel_del_microambiente=microambiente->indice_del_voxel_mas_cercano(posicion);

}

/**
 * @brief Sets internal constants for secretion/consumption calculations
 * @details Prepares the cell-specific constants used for calculating the cell's
 * impact on the microenvironment. This includes factors accounting for:
 * - Cell volume relative to voxel volume
 * - Secretion and consumption rates
 * - Saturation levels
 * 
 * In cancer research, this models how individual cells modify their local
 * biochemical environment, a critical aspect of tumor metabolism and
 * microenvironment evolution.
 * 
 * @param dt Time step size for the update
 */
void Celula::set_constantes_de_consumo_interno( double dt ){

	double constante_interna_para_discretizar_la_aproximacion_delta = dt * fenotipo.volumen.total / ( (microambiente->voxeles(voxel_del_microambiente).volumen)); // needs a fix

	// temp1 = dt*(V_cell/V_voxel)*S*T
	temp_celula_fuente_sumidero_solver1.assign( (fenotipo.secrecion.tasas_de_secrecion).size() , 0.0 );
	temp_celula_fuente_sumidero_solver1 += fenotipo.secrecion.tasas_de_secrecion;
	temp_celula_fuente_sumidero_solver1 *= fenotipo.secrecion.densidades_de_saturacion;
	temp_celula_fuente_sumidero_solver1 *= constante_interna_para_discretizar_la_aproximacion_delta;

	// temp2 = 1 + dt*(V_cell/V_voxel)*( S + U )
	temp_celula_fuente_sumidero_solver2.assign( (fenotipo.secrecion.tasas_de_secrecion).size() , 1.0 );
	axpy( &(temp_celula_fuente_sumidero_solver2) , constante_interna_para_discretizar_la_aproximacion_delta , fenotipo.secrecion.tasas_de_secrecion );
	axpy( &(temp_celula_fuente_sumidero_solver2) , constante_interna_para_discretizar_la_aproximacion_delta , fenotipo.secrecion.tasas_de_consumo );

	// temp para el siguiente export
	temp_celula_fuente_sumidero_exportacion_solver1 = fenotipo.secrecion.tasas_de_exportacion_neta;
	temp_celula_fuente_sumidero_exportacion_solver1 *= dt; // amount exported in dt of time

	temp_celula_fuente_sumidero_exportacion_solver2 = temp_celula_fuente_sumidero_exportacion_solver1;
	temp_celula_fuente_sumidero_exportacion_solver2 /= ( (microambiente->voxeles(voxel_del_microambiente)).volumen ) ;

	fenotipo.volumen.cambio_el_volumen = false;

	return;


}

/**
 * @brief Gets the cell's microenvironment
 * @return Pointer to the microenvironment
 */
Microambiente* Celula::get_microambiente( void ){

	return microambiente;

}

/**
 * @brief Gets the cell's voxel index in the microenvironment
 * @return Voxel index
 */
int Celula::get_indice_del_voxel_del_microambiente( void )
{
	return voxel_del_microambiente;
}

/**
 * @brief Gets the nearest density vector from the microenvironment
 * @details Provides access to the substrate concentrations in the cell's voxel
 * This is essential for cells to sense their biochemical environment.
 * 
 * In cancer research, this enables modeling:
 * - Cell response to hypoxia
 * - Nutrient-dependent proliferation
 * - Adaptation to microenvironmental conditions
 * 
 * @return Reference to the vector of substrate densities
 */
std::vector<double>& Celula::vector_de_densidades_mas_cercano( void )
{
	return microambiente->vector_de_densidades_mas_cercano( voxel_del_microambiente );
}


/**
 * @brief Gets the nearest gradient vector for a specific substrate
 * @details Provides access to the spatial gradient of a substrate at the cell's location
 * This is used for processes like chemotaxis.
 * 
 * In cancer research, this enables modeling:
 * - Directional cell migration in response to chemical gradients
 * - Cell polarization based on microenvironmental cues
 * 
 * @param indice_del_sustrato Index of the substrate in the microenvironment
 * @return Reference to the gradient vector
 */
std::vector<double>& Celula::gradiente_mas_cercano( int indice_del_sustrato ){

    return microambiente->vector_de_gradientes(voxel_del_microambiente)[indice_del_sustrato];

}


/**
 * @brief Simulates the cell's secretion and consumption processes
 * @details Updates the microenvironment based on the cell's secretion and consumption rates.
 * This method physically changes the substrate concentrations in the cell's voxel.
 * 
 * In cancer research, this models:
 * - Oxygen consumption by tumor cells
 * - Growth factor and cytokine secretion
 * - pH changes due to metabolic byproducts
 * - Nutrient depletion in densely packed tumors
 * 
 * @param dt Time step size for the update
 */
void Celula::simular_secrecion_y_consumo( double dt ){

	if(!es_activa)
	{ return; }

	if( fenotipo.volumen.cambio_el_volumen )
	{
		set_constantes_de_consumo_interno(dt);
		fenotipo.volumen.cambio_el_volumen = false;
	}


	// //Debug  immunostimulatory facor secretion
	// std::ofstream file("out/immuno_factor.csv", std::ios::app);
	// if (file.is_open()) {
	// file  << pg->tiempo_total << "," 
	// 	<< (*microambiente)(voxel_del_microambiente)[1]  << ","<< temp_celula_fuente_sumidero_solver1[1] <<","<< temp_celula_fuente_sumidero_solver2[1]<< "," << temp_celula_fuente_sumidero_exportacion_solver2[1]<< "\n";
	// }


	(*microambiente)(voxel_del_microambiente) += temp_celula_fuente_sumidero_solver1;
	(*microambiente)(voxel_del_microambiente) /= temp_celula_fuente_sumidero_solver2;

	(*microambiente)(voxel_del_microambiente) += temp_celula_fuente_sumidero_exportacion_solver2;


	


	return;


}

/**
 * @brief Initializes a cell with default parameters
 * @details Sets up the cell's cycle model, secretion/consumption rates, 
 * volume change rates, death pathways, and other key properties.
 * 
 * In cancer research, this method establishes:
 * - Baseline metabolic parameters (oxygen consumption)
 * - Immunostimulatory factor secretion for immune response modeling
 * - Oncoprotein expression levels (variable between cells)
 * - Proliferation parameters and oxygen sensitivity
 * - Apoptosis and necrosis pathways
 * 
 * This method creates heterogeneity in the cell population through
 * stochastic parameter assignment, a key aspect of realistic tumor modeling.
 */
void Celula::inicializar_celula(){
	fenotipo.ciclo.sync_con_ciclo_modelo(vida);
	// fenotipo.ciclo.sync_con_ciclo_modelo(Ki67);
	static int oxigeno_ID = microambiente->encontrar_indice_de_densidad( "oxigeno" ); // 0
	fenotipo.secrecion.tasas_de_secrecion[oxigeno_ID] = pg->tasas_de_secrecion;
	fenotipo.secrecion.tasas_de_consumo[oxigeno_ID] = pg->tasas_de_consumo;
	fenotipo.secrecion.densidades_de_saturacion[oxigeno_ID] = pg->densidades_de_saturacion;
	//FACTOR INMUNOESTIMULATORIO
    if( pg->activar_respuesta_inmune == true ){
    static int imm_ID = microambiente->encontrar_indice_de_densidad( "immunostimulatory factor" ); // 0
	fenotipo.secrecion.tasas_de_secrecion[imm_ID] = 10.0;
	fenotipo.secrecion.densidades_de_saturacion[imm_ID] = 1.0;
	}
	nombre= pg->c_nombre;


    fenotipo.secrecion.oncoproteina = rng->NormalRandom_CM(pg->imm_mean, pg->imm_sd);
	// fenotipo.secrecion.oncoproteina =1.;//Debug
    if( fenotipo.secrecion.oncoproteina < 0.0 ){
    fenotipo.secrecion.oncoproteina = 0.0; }

	parametros.o2_saturacion_para_la_proliferacion=pg->o2_saturacion_para_la_proliferacion;
	parametros.o2_referencia=pg->o2_referencia;
    
    fenotipo.volumen.citoplasma_tasa_de_cambio = 0.13/60.0;
	fenotipo.volumen.nucleo_tasa_de_cambio = 0.22/60.0;
	fenotipo.volumen.fluido_tasa_de_cambio = 1.3/60.0;
    
    fenotipo.ciclo.tasa_aleatoria= 1/(rng->NormalRandom_CM(38.6, 3.7)*60);
    fenotipo.ciclo.actualizar_mis_tasas_de_transicion(0,0) = fenotipo.ciclo.tasa_aleatoria;

	fenotipo.muerte.agregar_ciclo_de_muerte(0.0, &necrosis, necrosis_parametros);
	fenotipo.muerte.agregar_ciclo_de_muerte(0.0, &apoptosis, apoptosis_parametros);


	actualizar_voxel_del_microambiente();

	return;
}

/**
 * @brief Displays detailed information about the cell
 * @details Prints the cell's type, name, physical properties (radius, volume),
 * position, voxel indices, cycle model details, and secretion/consumption rates.
 * 
 * In cancer research, this method is valuable for:
 * - Debugging simulation behavior
 * - Tracking individual cell characteristics
 * - Analyzing heterogeneity in the tumor cell population
 * - Validating cell parameter initialization
 * 
 * @param os Output stream to write the information to
 */
void Celula::mostrar_informacion_de_la_celula(std::ostream& os){

	os << " tipo:" << tipo << " nombre: " << nombre << std::endl;
	os << " radio:" << fenotipo.geometria.radio << " volumen: " << fenotipo.volumen.total << std::endl;
	os << " posicion: " << posicion << std::endl;
	os << " voxel_M: " << voxel_del_microambiente << std::endl;
	os << " voxel_C: " << voxel << std::endl;
	//Resumen del ciclo
	os << " cycle model: " << fenotipo.ciclo.pCiclo_Modelo->nombre
				<< " (codigo=" << fenotipo.ciclo.pCiclo_Modelo->codigo << ")" << std::endl;
	fenotipo.ciclo.pCiclo_Modelo->mostrar_ciclo(std::cout);

	//Secreciones
	for( unsigned int i = 0 ; i < fenotipo.secrecion.tasas_de_secrecion.size() ; i++ )
	{
		os << " tasa de secrecion =" << fenotipo.secrecion.tasas_de_secrecion[i] << std::endl;
		os << " tasa de consumo =" << fenotipo.secrecion.tasas_de_consumo[i] << std::endl;
		os << " densidades_de_saturacion =" << fenotipo.secrecion.densidades_de_saturacion[i] << std::endl;}

		return;

}


/**
 * @brief Initializes a healthy cell with non-proliferative parameters
 * @details Sets up a cell with healthy cell parameters, including turning off
 * proliferation and setting appropriate oxygen consumption rates. The cell is
 * set to type 1 (healthy) and given the name "Sana".
 * 
 * In cancer research, this method enables:
 * - Modeling the normal tissue surrounding tumors
 * - Establishing baseline healthy cell behavior for comparison
 * - Creating heterogeneous tissue environments with both normal and cancer cells
 * - Studying cancer cell interactions with healthy tissue
 * 
 * This is crucial for realistic tumor microenvironment simulation.
 */
void Celula::inicializar_celula_sana(){
	fenotipo.ciclo.sync_con_ciclo_modelo(vida);

////////Apagar la proliferación
	int indice_inicio_ciclo = vida.encontrar_indice_de_la_fase(c->viva);
	int indice_fin_ciclo = vida.encontrar_indice_de_la_fase(c->viva);
	fenotipo.ciclo.pCiclo_Modelo->tasa_de_transicion(indice_inicio_ciclo, indice_fin_ciclo) = 0.0;
	fenotipo.ciclo.tasas_de_transicion = fenotipo.ciclo.pCiclo_Modelo->tasas_de_transicion;
////////

	static int oxigeno_ID = microambiente->encontrar_indice_de_densidad( "oxigeno" ); // 0
	fenotipo.secrecion.tasas_de_secrecion[oxigeno_ID] = 0.0;
	fenotipo.secrecion.tasas_de_consumo[oxigeno_ID] = 1.0;
	nombre= "Sana";
	tipo = 1;

	//fenotipo.mecanica.fuerza_de_adhesion_cc = 2.0;

    fenotipo.muerte.agregar_ciclo_de_muerte(0.0, &necrosis, necrosis_parametros);

	actualizar_voxel_del_microambiente();

	return;
}

/**
 * @brief Constructor for the Linfocito (Lymphocyte) class
 * @details Initializes a lymphocyte immune cell with specialized parameters for:
 * - Migration behavior (motility)
 * - Cell-cell adhesion and interaction
 * - Cytotoxic activity against tumor cells
 * - Limited lifespan
 * 
 * In cancer research, this class is crucial for modeling:
 * - Anti-tumor immune responses
 * - Immunotherapy mechanisms
 * - Tumor immune evasion
 * - Immune cell trafficking and infiltration
 * 
 * The lymphocyte maintains specialized parameters for detecting and killing
 * tumor cells based on their oncoprotein expression levels.
 */
Linfocito::Linfocito():Celula(){

	nombre= "Linfocito";
	tipo = 2;
    pg->numero_id += 1;
    id=pg->numero_id;
	fenotipo.ciclo.sync_con_ciclo_modelo(vida);

////////Apagar la proliferación
	int indice_inicio_ciclo = vida.encontrar_indice_de_la_fase(c->viva);
	int indice_fin_ciclo = vida.encontrar_indice_de_la_fase(c->viva);
    fenotipo.ciclo.actualizar_mis_tasas_de_transicion(indice_inicio_ciclo, indice_fin_ciclo) = 0.0;
//	fenotipo.ciclo.pCiclo_Modelo->tasa_de_transicion(indice_inicio_ciclo, indice_fin_ciclo) = 0.0;
//	fenotipo.ciclo.tasas_de_transicion = fenotipo.ciclo.pCiclo_Modelo->tasas_de_transicion;
////////

	static int oxigeno_ID = get_microambiente()->encontrar_indice_de_densidad( "oxigeno" ); // 0
	fenotipo.secrecion.tasas_de_secrecion[oxigeno_ID] = 0.0;
	fenotipo.secrecion.tasas_de_consumo[oxigeno_ID] = 1.0;

    // set apoptosis to survive 10 days (on average): 1.0 / (10.0 * 24.0 * 60.0 )
	fenotipo.muerte.agregar_ciclo_de_muerte(1.0 / (c->dt_ciclo * 10.0 * 24.0 * 60.0 ), &apoptosis, apoptosis_parametros);

	tasa_de_asesinato = 0.06667; // 1/min // how often it tries to kill
    tiempo_de_adhesion = 60.0; // min
    tasa_de_adhesion = 0.2; // 1/min
    constante_elastica = 0.01; // 1/min
    distancia_de_adhesion_maxima = 18.0; // micrones
    distancia_de_adhesion_minima = 14.0; // micrones
    saturacion_de_oncoproteina = 2.0;
    limite_de_oncoproteina = 0.5;
    diferencia_de_oncoproteina = saturacion_de_oncoproteina - limite_de_oncoproteina;
    diferencia_de_adhesion = distancia_de_adhesion_maxima - distancia_de_adhesion_minima;
    hora_de_la_ultima_mecanica = pg->tiempo_total;
    hora_del_ultimo_ciclo = pg->tiempo_total;



	motilidad.es_movil = true;
	motilidad.tiempo_de_persistencia = 10.0;
	motilidad.velocidad_de_migracion = 5.0;
	motilidad.bias_de_la_migracion = 0.5;


	fenotipo.mecanica.fuerza_de_adhesion_cc = 0.0;
	fenotipo.mecanica.fuerza_de_adhesion_co = 0.0;
	fenotipo.mecanica.fuerza_de_adhesion_mb = 0.0;

	fenotipo.mecanica.fuerza_de_repulsion_cc *= 5.0;//5.0
	fenotipo.mecanica.fuerza_de_repulsion_co = fenotipo.mecanica.fuerza_de_repulsion_cc;



	actualizar_voxel_del_microambiente();

	return;


}

/**
 * @brief Sets the motility state of the lymphocyte
 * @details Controls whether the lymphocyte can actively migrate or not.
 * This is used when lymphocytes are attached to target cells or
 * during other states where movement should be restricted.
 * 
 * In cancer research, this enables modeling:
 * - Immune cell arrest during engagement with tumor cells
 * - Transition between scanning and effector functions
 * 
 * @param valor Boolean value indicating whether the cell should be motile
 */
void Linfocito::es_movil(bool valor){
    motilidad.es_movil = valor;
    return;
}

/**
 * @brief Updates the lymphocyte's motility vector
 * @details Calculates a new migration direction based on:
 * - Random motion component
 * - Chemotactic bias (toward immunostimulatory factors)
 * - Current motility state
 * - Time persistence
 * 
 * In cancer research, this models:
 * - Immune cell chemotaxis toward tumor sites
 * - Scanning behavior within the tumor microenvironment
 * - Lymphocyte trafficking and infiltration dynamics
 * 
 * @param dt Time step size
 * @param celulas_en_mi_voxel Vector of cell pointers sharing the same voxel
 */
void Linfocito::actualizar_vector_de_motilidad( double dt, std::vector<Celula*> celulas_en_mi_voxel ){

	if( motilidad.es_movil == false )
	{
		motilidad.vector_de_motilidad.x = 0.0;
		motilidad.vector_de_motilidad.y = 0.0;
		motilidad.vector_de_motilidad.z = 0.0;
		return;
	}
	// std::cout <<pg->tiempo_total <<"velocity: " << velocidad.x <<","<< velocidad.y <<","<< velocidad.z <<",norma "<< std::sqrt(velocidad.x*velocidad.x+velocidad.y*velocidad.y+velocidad.z*velocidad.z) << std::endl;//Debug

	// if( true )//Debug
	if( rng->RandomNumber() < dt / motilidad.tiempo_de_persistencia || motilidad.tiempo_de_persistencia < dt )//Debug needs tu be uncommented
	{
	//std::cout << "Actualizo vector motilidarks \n";
	//std::cin.get();
		// choose a uniformly random unit vector
		double temp_angle = 6.28318530717959*rng->RandomNumber();
		double temp_phi = 3.1415926535897932384626433832795*rng->RandomNumber();

		double sin_phi = sin(temp_phi);
		double cos_phi = cos(temp_phi);

		Vector randvec;
		randvec.x = sin_phi;
		randvec.y = sin_phi;
		randvec.z = sin_phi;

		randvec.x *= cos( temp_angle ); // cos(theta)*sin(phi)
		randvec.y *= sin( temp_angle ); // sin(theta)*sin(phi)
		randvec.z = cos_phi; //  cos(phi)

		motilidad_de_linfocito( dt, celulas_en_mi_voxel  ); //////////////////////////PROGRAMAR
		motilidad.vector_de_motilidad = motilidad.bias_de_la_migracion_direccion;
		motilidad.vector_de_motilidad = motilidad.vector_de_motilidad * motilidad.bias_de_la_migracion;
		double uno_menos_bias = 1.0 - motilidad.bias_de_la_migracion;

		axpy( &(motilidad.vector_de_motilidad), uno_menos_bias, randvec );
		v->normalizame( &(motilidad.vector_de_motilidad) );
		motilidad.vector_de_motilidad = motilidad.vector_de_motilidad * motilidad.velocidad_de_migracion;
		velocidad = velocidad + motilidad.vector_de_motilidad;
// std::cout <<pg->tiempo_total <<"velocity: after motility" << velocidad.x <<","<< velocidad.y <<","<< velocidad.z <<",norma "<< std::sqrt(velocidad.x*velocidad.x+velocidad.y*velocidad.y+velocidad.z*velocidad.z) << std::endl;//Debug
		
	}
	return;

}

/**
 * @brief Determines the lymphocyte's chemotactic response
 * @details Sets the direction of biased migration based on:
 * - Gradients of immunostimulatory factors
 * - Presence of nearby cells
 * - Cell type in the current voxel
 * 
 * In cancer research, this models:
 * - Directed migration toward inflammatory signals
 * - Immune cell arrest during target cell engagement
 * - Lymphocyte behavior in different tumor regions
 * 
 * @param dt Time step size
 * @param celulas_en_mi_voxel Vector of cell pointers sharing the same voxel
 */
void Linfocito::motilidad_de_linfocito( double dt, std::vector<Celula*> celulas_en_mi_voxel  )
{



	static int indice_del_factor_inmune = get_microambiente()->encontrar_indice_de_densidad("immunostimulatory factor");

	if( celulas_en_mi_voxel.size() == 1 )
	{


		motilidad.es_movil = true;
		motilidad.bias_de_la_migracion_direccion.x = gradiente_mas_cercano(indice_del_factor_inmune)[0];
		motilidad.bias_de_la_migracion_direccion.y = gradiente_mas_cercano(indice_del_factor_inmune)[1];
		motilidad.bias_de_la_migracion_direccion.z = gradiente_mas_cercano(indice_del_factor_inmune)[2];
		v->normalizame( &( motilidad.bias_de_la_migracion_direccion ) );
	}
	else if(celulas_en_mi_voxel.size() > 1)
	{

        bool todos_linfos = true;
        for(unsigned int i=1; i < celulas_en_mi_voxel.size(); i++){
            if(celulas_en_mi_voxel[i]->tipo == 1){
                todos_linfos = false;
                motilidad.es_movil = false;
                continue;
            }
        }
            if(todos_linfos == true){
                motilidad.es_movil = true;
                motilidad.bias_de_la_migracion_direccion.x = gradiente_mas_cercano(indice_del_factor_inmune)[0];
                motilidad.bias_de_la_migracion_direccion.y = gradiente_mas_cercano(indice_del_factor_inmune)[1];
                motilidad.bias_de_la_migracion_direccion.z = gradiente_mas_cercano(indice_del_factor_inmune)[2];
                v->normalizame( &( motilidad.bias_de_la_migracion_direccion ) );

            }
	}

	return;
}


/**
 * @brief Controls lymphocyte behavior including attachment, detachment, and killing
 * @details Main method governing lymphocyte interaction with target cells, including:
 * - Force-based movement toward attached cells
 * - Attempts to induce apoptosis in tumor cells
 * - Detachment after killing or timeout
 * - Seeking new target cells
 * 
 * In cancer research, this models:
 * - Cytotoxic T-cell killing mechanisms
 * - Immune synapse formation and dynamics
 * - Immune cell serial killing behavior
 * - Tumor cell-immune cell interactions
 * 
 * @param dt Time step size
 * @param celulas_en_mi_voxel Vector of cell pointers sharing the same voxel
 */
void Linfocito::avanzar_linfocito( double dt, std::vector<Celula*> celulas_en_mi_voxel )
{


	if( fenotipo.muerte.muerta == true ){
		return;
	}

	// std::cout << adherida<< motilidad.es_movisl << std::endl;//Debug

	//Debug CAR-T positions
	std::ofstream file("out/posiciones_cart.csv", std::ios::app);
	if (file.is_open()) {
	file  << pg->tiempo_total << "," 
		<< posicion.x << "," << posicion.y << "," << posicion.z << ", norma " << std::sqrt(posicion.x*posicion.x+posicion.y*posicion.y+posicion.z*posicion.z) << "\n";
	}


	if( celulas_en_mi_voxel.size() > 1 )
	{
		for(long unsigned int i=0; i < celulas_en_mi_voxel.size(); i++ ){
            if(celulas_en_mi_voxel[i] != this && celulas_en_mi_voxel[i]->tipo != 2){

            Vector desplazamiento = celulas_en_mi_voxel[i]->posicion - posicion;
            
         	//  std::cout <<pg->tiempo_total <<"Movement towards tumor cell: " << desplazamiento.x*constante_elastica << std::endl;//Debug
			
			
			axpy( &(velocidad) , constante_elastica , desplazamiento );
			//  std::cout <<pg->tiempo_total <<"Movement towards tumor cell: " << velocidad.x <<","<< velocidad.y <<","<< velocidad.z <<",norma "<< std::sqrt(velocidad.x*velocidad.x+velocidad.y*velocidad.y+velocidad.z*velocidad.z) << std::endl;//Debug
			}
		}



		bool soltarme = false;

		if( adherida && intento_de_apoptosis( celula_adherida, dt ) )
		{
			// std::cout << "Linfocito " << id << " desencadena apoptosis en celula " << celula_adherida->id << std::endl;//Debug
			desencadenar_apoptosis( celula_adherida );

			soltarme = true;

		}




		if( adherida && rng->RandomNumber() < dt / ( tiempo_de_adhesion + 1e-15 ) )
		{ soltarme = true; }



		if( soltarme )
		{

			soltar_celula( celula_adherida );

			motilidad.es_movil = true;
		}







	if( adherida==false){

	// Debug uncomment
	bool check = chequear_vecinos_para_adherirse( celulas_en_mi_voxel , dt );

		if (check == true){
		motilidad.es_movil = false;
		return;
		}
	}
	motilidad.es_movil = true;
	}

	return;
}

/**
 * @brief Determines if a lymphocyte should induce apoptosis in a target cell
 * @details Makes a probabilistic decision to kill a target cell based on:
 * - Target cell's oncoprotein expression level
 * - Lymphocyte's killing rate
 * - Time step size
 * 
 * In cancer research, this models:
 * - Immune recognition of cancer cells
 * - Stochastic nature of cytotoxic activity
 * - Tumor immunogenicity and immune evasion
 * - Tumor antigen threshold effects
 * 
 * @param celula_adherida Pointer to the attached target cell
 * @param dt Time step size
 * @return Boolean indicating whether killing should occur
 */
bool Linfocito::intento_de_apoptosis( Celula* celula_adherida, double dt )
{
	//static int oncoproteina_i = celula_adherida->fenotipo.secrecion.oncoproteina;
	//static int indice_del_ciclo_apoptosis = celula_adherida->fenotipo.muerte.encontrar_indice_del_ciclo_de_muerte("apoptosis");
	//static int kill_rate_index = pAttacker->custom_data.find_variable_index( "kill rate" );



	//static double oncoprotein_saturation = parameters.doubles("oncoprotein_saturation"); // 2.0;
	//static double oncoprotein_threshold = parameters.doubles("oncoprotein_threshold"); // 0.5; // 0.1;
	//static double diferencia_de_oncoproteina = saturacion_de_oncoproteina - limite_de_oncoproteina;



	if( celula_adherida->fenotipo.secrecion.oncoproteina < limite_de_oncoproteina )
	{ return false; }


	double escala = celula_adherida->fenotipo.secrecion.oncoproteina;
	escala -= limite_de_oncoproteina;
	escala /= diferencia_de_oncoproteina;
	if( escala > 1.0 )
	{ escala = 1.0; }



	if( rng->RandomNumber() < tasa_de_asesinato * escala * dt )
	{
		return true;
	}
	return false;
}


/**
 * @brief Initiates apoptosis in a target cell
 * @details Activates the apoptosis death pathway in the target cell.
 * 
 * In cancer research, this models:
 * - Cytotoxic T-cell induced tumor cell death
 * - Immune-mediated tumor control
 * - Programmed cell death mechanisms
 * 
 * @param celula_adherida Pointer to the target cell
 * @return Boolean indicating success of apoptosis initiation
 */
bool Linfocito::desencadenar_apoptosis( Celula* celula_adherida )
{
	static int indice_del_ciclo_apoptosis = celula_adherida->fenotipo.muerte.encontrar_indice_del_ciclo_de_muerte("Apoptosis");



	if( celula_adherida->fenotipo.muerte.muerta == true ){

        return false; }

	celula_adherida->comenzar_muerte(indice_del_ciclo_apoptosis);


	return true;
}

/**
 * @brief Creates an attachment between the lymphocyte and a target cell
 * @details Establishes a bidirectional link between the lymphocyte and its target,
 * setting the adherida (attached) flag in both cells and cross-referencing them.
 * 
 * In cancer research, this models:
 * - Immune synapse formation
 * - T-cell engagement with tumor cells
 * - Physical coupling preceding cytotoxic action
 * 
 * @param celula_objetivo Pointer to the target cell to attach to
 */
void Linfocito::adherir_celula(Celula* celula_objetivo){

    if(adherida==false && celula_objetivo->adherida == false){
        celula_adherida = celula_objetivo;
        celula_objetivo->adherida = true;
        adherida = true;
        celula_objetivo->celula_adherida = this;
    }

 return;
}

/**
 * @brief Breaks the attachment between the lymphocyte and a target cell
 * @details Releases the bidirectional link between the lymphocyte and its target,
 * clearing the adherida (attached) flag and references in both cells.
 * 
 * In cancer research, this models:
 * - Immune synapse dissolution
 * - T-cell detachment after killing or unsuccessful engagement
 * - Serial killing behavior of cytotoxic T-cells
 * 
 * @param celula_objetivo Pointer to the target cell to detach from
 */
void Linfocito::soltar_celula(Celula* celula_objetivo){

    if(adherida==true && celula_objetivo->adherida == true){
        celula_objetivo->adherida = false;
        celula_objetivo->celula_adherida = NULL;
        celula_adherida = NULL;
        adherida = false;
        
    }

 return;
}

/**
 * @brief Scans neighboring cells for potential targets to attach to
 * @details Evaluates cells in the same voxel as candidates for attachment,
 * based on oncoprotein expression levels and proximity.
 * 
 * In cancer research, this models:
 * - Immune cell surveillance
 * - Target recognition based on tumor-associated antigens
 * - Local scanning behavior of immune cells
 * 
 * @param celulas_en_mi_voxel Vector of cell pointers sharing the same voxel
 * @param dt Time step size
 * @return Boolean indicating whether attachment was successful
 */
bool Linfocito::chequear_vecinos_para_adherirse( std::vector<Celula*> celulas_en_mi_voxel , double dt )
{
    unsigned int i = 0;
	while( i < celulas_en_mi_voxel.size() && adherida == false)
	{



		if( celulas_en_mi_voxel[i] != this )
		{
			if( intentar_adherirse( celulas_en_mi_voxel[i] , dt ) )
			{

				return true;
			}
		}
		i= i+1;

	}

	return false;
}

/**
 * @brief Attempts to adhere to a specific target cell
 * @details Determines whether attachment should occur based on:
 * - Target cell's oncoprotein expression level
 * - Distance between the cells
 * - Random chance modulated by adhesion rate
 * 
 * In cancer research, this models:
 * - Antigen-dependent recognition of tumor cells
 * - Spatial constraints in immune-tumor interaction
 * - Probabilistic nature of immune synapse formation
 * 
 * @param celula_objetivo Pointer to the potential target cell
 * @param dt Time step size
 * @return Boolean indicating whether attachment was successful
 */
bool Linfocito::intentar_adherirse( Celula* celula_objetivo , double dt )
{


	if( celula_objetivo->fenotipo.secrecion.oncoproteina > limite_de_oncoproteina && celula_objetivo->fenotipo.muerte.muerta == false && celula_objetivo->adherida == false )
	{
		Vector desplazamiento = celula_objetivo->posicion - posicion;
		double escala_de_distancia = norma(desplazamiento);

		if( escala_de_distancia > distancia_de_adhesion_maxima ){

            return false;
        }

		double escala = celula_objetivo->fenotipo.secrecion.oncoproteina;
		escala -= limite_de_oncoproteina;
		escala /= diferencia_de_oncoproteina;
		if( escala > 1.0 )
		{ escala = 1.0; }


		escala_de_distancia *= -1.0;
		escala_de_distancia += distancia_de_adhesion_maxima;
		escala_de_distancia /= diferencia_de_adhesion;
		if( escala_de_distancia > 1.0 )
		{ escala_de_distancia = 1.0; }


		if( rng->RandomNumber() < tasa_de_adhesion * escala * dt * escala_de_distancia )
		{

            adherir_celula(celula_objetivo);

            return true;
		}


	}


	return false;
}