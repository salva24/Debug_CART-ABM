/**
 * @file Muerte_Parametros.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Cell death parameters for cancer simulation
 *
 * @details This file defines the Muerte_parametros (Death Parameters) class, which
 * encapsulates the parameters governing cell death processes in the cancer simulation.
 * These parameters control how cells change their physical properties during death,
 * including volume changes, fluid exchange, and structural degradation.
 * 
 * In cancer research, accurate modeling of cell death is critical for:
 * - Simulating tumor necrotic core formation in poorly vascularized regions
 * - Modeling treatment response and cell death mechanisms (apoptosis vs. necrosis)
 * - Representing calcification processes in certain cancer types
 * - Understanding the impact of cell debris on the tumor microenvironment
 * - Modeling immune-mediated cancer cell killing
 * 
 * Inbound dependencies string
 *
 * Outbound dependencies Muerte.h, Ciclo.h
 *
 * @usage Used by death pathways to control the physical manifestation of cell death
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __MUERTE_PARAMETROS_H__
#define __MUERTE_PARAMETROS_H__

#include <string>

/**
 * @class Muerte_parametros
 * @brief Parameters governing cell death processes in cancer simulation
 * 
 * @details This class defines the rates and thresholds that control how cells
 * physically change during death processes. These parameters determine fluid
 * exchange rates, volume changes, and structural degradation during different
 * death pathways.
 * 
 * In cancer research, these parameters enable modeling of:
 * - Different death mechanisms (apoptosis vs. necrosis) with distinct physical manifestations
 * - Formation of necrotic regions in poorly vascularized tumor areas
 * - Calcification processes observed in certain cancer types
 * - Cell lysis and release of intracellular contents into the microenvironment
 * - Morphological changes during programmed cell death
 */
class Muerte_parametros{
	
	public:
		/**
		 * @brief Time units for rate parameters
		 * @details Specifies the time units (typically "min") for all rate parameters.
		 * Ensures consistent time scales across the simulation.
		 */
		std::string tiempo_unidades;
		
		/**
		 * @brief Rate of fluid change in non-lysed cells
		 * @details Controls how quickly dying cells exchange fluid with their environment
		 * before membrane rupture. In cancer, this represents early volume changes during
		 * apoptosis or controlled cell death.
		 */
		double tasa_de_cambio_fluido_no_lisado;
		
		/**
		 * @brief Rate of fluid change in lysed cells
		 * @details Controls fluid exchange after membrane rupture. In cancer modeling,
		 * this represents the slower equilibration phase after cell lysis, affecting
		 * how quickly cell debris is processed in the tumor microenvironment.
		 */
		double tasa_de_cambio_fluido_lisado;
		
		/**
		 * @brief Rate of cytoplasmic degradation
		 * @details Controls how quickly the cytoplasm breaks down during cell death.
		 * In cancer research, this parameter affects how long dead cells persist in
		 * the tumor microenvironment before complete degradation.
		 */
		double citoplasma_tasa_de_cambio;
		
		/**
		 * @brief Rate of nuclear degradation
		 * @details Controls how quickly the nucleus breaks down during cell death.
		 * In cancer modeling, the typically slower nuclear degradation represents
		 * the persistence of nuclear fragments that can be observed in histological
		 * samples of tumors.
		 */
	    double nucleo_tasa_de_cambio;
	    
	    /**
	     * @brief Rate of cell calcification
	     * @details Controls the rate at which dying cells undergo calcification.
	     * In cancer research, this models microcalcifications observed in certain
	     * cancer types, particularly breast cancer, which are important diagnostic
	     * markers in mammography.
	     */
	    double tasa_de_calcificacion;
	    
	    /**
	     * @brief Relative volume threshold for cell rupture
	     * @details Defines the volume ratio (relative to normal) at which cell membrane
	     * ruptures during death. In cancer modeling, this threshold determines when
	     * necrotic cells release their contents, potentially triggering inflammatory
	     * responses in the tumor microenvironment.
	     */
	    double volumen_de_ruptura_relativo;
	    
	    /**
	     * @brief Default constructor with cancer-relevant death parameters
	     * @details Initializes death parameters with default values calibrated for
	     * cancer cell simulation. Values are based on experimental observations of
	     * cancer cell death processes.
	     */
	    Muerte_parametros();
	
};


#endif
