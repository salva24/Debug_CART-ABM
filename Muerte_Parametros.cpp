/**
 * @file Muerte_Parametros.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of cell death parameters for cancer simulation
 *
 * @details This file implements the Muerte_parametros (Death Parameters) class, which
 * defines the default parameter values governing cell death processes in the cancer
 * simulation. These parameters control volume changes, fluid exchange rates, and
 * structural degradation during different death pathways.
 * 
 * The default values are calibrated based on experimental observations of cancer cell
 * death, particularly referencing MCF-7 breast cancer cells, to provide realistic
 * modeling of cell death processes in tumors.
 * 
 * Inbound dependencies Muerte_Parametros.h
 *
 * Outbound dependencies None
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#include "Muerte_Parametros.h"

/**
 * @brief Default constructor with cancer-relevant death parameters
 * @details Initializes death parameters with default values calibrated for
 * cancer cell simulation, with specific reference to MCF-7 breast cancer cells.
 * 
 * The parameters include:
 * - Time units: minutes, for consistency with other simulation time scales
 * - Fluid exchange rates: Different rates for pre-lysis (3.0/60.0 min⁻¹) and 
 *   post-lysis (0.05/60.0 min⁻¹) phases, reflecting the faster initial volume
 *   changes followed by slower equilibration
 * - Cytoplasmic degradation rate (1.0/60.0 min⁻¹): Relatively rapid breakdown
 *   of cytoplasmic contents during cell death
 * - Nuclear degradation rate (0.35/60.0 min⁻¹): Slower degradation of nuclear
 *   material, consistent with observations of persistent nuclear fragments in
 *   dying cancer cells
 * - Calcification rate (0.0 by default): Can be increased to model microcalcifications
 *   observed in certain cancer types
 * - Rupture volume threshold (2.0x normal volume): Defines when necrotic cells
 *   rupture and release their contents into the tumor microenvironment
 * 
 * These values collectively model the physical manifestation of cell death in
 * cancer, with particular relevance to necrotic regions in poorly vascularized
 * tumors and treatment-induced cell death.
 */
Muerte_parametros::Muerte_parametros(){
	
		tiempo_unidades = "min";
		

		tasa_de_cambio_fluido_no_lisado = 3.0/60.0;
		tasa_de_cambio_fluido_lisado = 0.05/60.0;
		
		citoplasma_tasa_de_cambio = 1.0/60.0;
	    nucleo_tasa_de_cambio = 0.35/60.0;
	    
	    tasa_de_calcificacion = 0.0;
	    
	    volumen_de_ruptura_relativo = 2.0;

	
}
