/**
 * @file Volumen.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of the cell volume class for cancer modeling
 *
 * @details This file implements the Volumen class methods, establishing default
 * volume values based on MCF-7 breast cancer cells. It provides functionality
 * for volume manipulation (division and multiplication) needed during
 * cell growth, division, and other volume-changing processes.
 * 
 * @inbound_dependencies Volumen.h
 * 
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#include "Volumen.h"

/**
 * @brief Constructor that initializes volume parameters based on MCF-7 breast cancer cells
 * 
 * @details Initializes all volume components with reference values from MCF-7 breast cancer cells.
 * MCF-7 is a widely used breast cancer cell line in cancer research, serving as a standard
 * model for estrogen receptor-positive breast adenocarcinoma.
 * 
 * Key cancer research relevance:
 * - The default total volume (2494 cubic microns) represents typical MCF-7 cell measurements
 * - The fluid fraction (0.75) is characteristic of cancer cells, which often have higher
 *   water content than normal cells, affecting drug diffusion
 * - The nuclear volume (540 cubic microns) reflects the enlarged nuclei typical of cancer cells
 * - The cytoplasmic-to-nucleus ratio (~3.6) is lower than in normal cells, consistent with 
 *   the high nuclear-to-cytoplasmic ratio characteristic of cancer cells
 * 
 * The rate parameters model how quickly different cell components can change in volume,
 * which is particularly relevant for cancer cells that often exhibit dysregulated growth.
 * 
 * @return void
 */
Volumen::Volumen(){
	
	//Parámetros de referencia en micrómetros

	fraccion_de_fluido=0.75;
	
	total=2494.0;
	fluido = fraccion_de_fluido * total;
	solido = total- fluido;
	
	
	nuclear=540.0;
	nuclear_fluido = fraccion_de_fluido * nuclear;
	nuclear_solido = nuclear - nuclear_fluido;
	    
	citoplasmatico = total - nuclear;
    citoplasmatico_fluido = fraccion_de_fluido * citoplasmatico;
	citoplasmatico_solido = citoplasmatico - citoplasmatico_fluido;
	    
	fraccion_calcificada = 0.0;
	    
	relacion_citoplasma_nucleo = citoplasmatico / ( 1e-16 + nuclear);
	    
	citoplasma_tasa_de_cambio = 0.27/60.0;
	nucleo_tasa_de_cambio = 0.33/60.0;
	fluido_tasa_de_cambio = 3.0 / 60.0;
	    
	tasa_de_calcificacion = 0.0;
	    
	target_citoplasma_solido = citoplasmatico_solido;
	target_nucleo_solido = nuclear_solido;
	target_fraccion_fluido = fraccion_de_fluido;
	target_relacion_citoplasma_nucleo = relacion_citoplasma_nucleo;
	    
	volumen_de_ruptura_relativo = 2.0;
	volumen_de_ruptura = volumen_de_ruptura_relativo * total;
	    
	cambio_el_volumen = true;	    
	
	// Commented out code represents previous implementation approach
	// total = 2494;
	// tasa_de_cambio = -log(1-0.99);
	// target = 2494;

	
	return;
}

/**
 * @brief Multiplies all volume components by a specified factor
 * 
 * @details Scales all volume components uniformly by the specified factor. This method
 * is used to model uniform cell growth or shrinkage, which occurs during various
 * cellular processes. In cancer research, this is particularly relevant for modeling:
 * - Gradual cell growth during the cell cycle
 * - Cell swelling during certain stress responses
 * - Cell shrinkage during apoptosis and other cell death processes
 * - Response to osmotic challenges in the tumor microenvironment
 * 
 * All components (total, solid, fluid, nuclear, cytoplasmic) are scaled proportionally,
 * preserving their relative ratios while changing absolute values.
 * 
 * @param numero The multiplication factor to apply to all volume components
 * @return void
 */
void Volumen::multiplicar(double numero){
	
	total *= numero;
	solido *= numero;
	fluido *= numero;
	
	nuclear *= numero;
	nuclear_fluido *= numero;
	nuclear_solido *= numero;
	    
	citoplasmatico *= numero;
	citoplasmatico_fluido *= numero;
	citoplasmatico_solido *= numero;
	    
	volumen_de_ruptura *= numero;
	    
	target_citoplasma_solido *= numero;
	target_nucleo_solido *= numero;
	
	
	// Commented out code represents previous implementation approach
	// total *= numero;
	// target *= numero;
	
	return;
	
}

/**
 * @brief Divides all volume components in half
 * 
 * @details Models cell division by halving all volume components. In cancer research,
 * this is particularly relevant as:
 * - Cancer cells often undergo rapid and frequent divisions
 * - Volume division is a key step in the cell cycle completion
 * - Abnormal volume distribution during division can lead to aneuploid cells,
 *   a common feature in cancer
 * 
 * This method is a specialized case of multiplicar() with a factor of 0.5,
 * used specifically for modeling mitotic division.
 * 
 * @return void
 */
void Volumen::dividir(void){
	
	multiplicar(0.5);
	return;
}

