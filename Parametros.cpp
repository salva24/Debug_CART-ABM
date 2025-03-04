/**
 * @file Parametros.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of oxygen parameter definitions for cancer cell simulation
 *
 * @details Implements the Parametros class that defines oxygen thresholds for cancer cell
 * behaviors, particularly proliferation and necrosis. Default values are set based on
 * experimental observations of cancer cell responses to varying oxygen levels.
 * 
 * Inbound Dependencies: Parametros.h, <string>
 *
 * Outbound Dependencies: None
 * 
 * @license: Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#include "Parametros.h"

#include <string>

/**
 * @brief Constructor that initializes oxygen threshold parameters
 * @details Sets default values for oxygen-dependent behaviors based on experimental data from
 * cancer research. These parameters model how cancer cells respond to varying oxygen levels
 * in the tumor microenvironment, including proliferation thresholds and necrosis triggers.
 * 
 * Key values are set to model:
 * 1. Severe hypoxia (<2.5 mmHg): Maximum necrosis rate
 * 2. Hypoxia (2.5-5 mmHg): Increasing necrosis rate
 * 3. Restricted proliferation (5-10 mmHg): Cells viable but non-proliferating
 * 4. Optimal proliferation (>38 mmHg): Cells proliferate at maximum rate
 * 
 * @return None
 */
Parametros::Parametros(){

	//o2_hypoxia_limite = 15.0;
	//o2_hypoxia_respuesta = 8.0;
	//o2_hypoxia_saturacion = 4.0;

	o2_necrosis_limite = 5.0;//5.0;
	o2_necrosis_max = 2.5;//2.5;

	o2_limite_de_proliferacion = 10.0;//5.0;
	o2_saturacion_para_la_proliferacion = 38.0; //38.0;
	o2_referencia = 38.0; //38.0;

	tasa_necrosis_max = 1.0 / (6.0 * 60.0);

	//pReferencia_a_un_fenotipo_vivo = &fenotipo;

	return;


}
