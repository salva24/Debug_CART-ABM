/**
 * @file Secrecion.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of secretion and consumption behaviors for cancer cells
 *
 * @details This file implements the Secrecion class methods that model how cancer cells
 * interact with their microenvironment through secretion and consumption of various
 * biochemical substances. These interactions are critical for understanding:
 * - Cancer cell influence on their surrounding tissue
 * - Tumor-induced angiogenesis and immune suppression
 * - Metabolic reprogramming in cancer
 * - Formation of ecological niches within tumors
 * 
 * @inbound_dependencies "Secrecion.h"
 *
 * @outbound_dependencies None
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#include "Secrecion.h"

/**
 * @brief Constructor initializes all secretion and consumption vectors to empty
 * @details Cancer cells often have substrate-specific secretion and consumption profiles
 * that differ from normal cells. This constructor creates the data structure framework
 * that will later be populated with cancer-specific values depending on phenotype.
 */
Secrecion::Secrecion()
{
	tasas_de_secrecion.resize( 0 , 0.0 );
	tasas_de_consumo.resize( 0 , 0.0 );
	densidades_de_saturacion.resize( 0 , 0.0 );
	tasas_de_exportacion_neta.resize( 0 , 0.0 );
	oncoproteina = 0.0;
	return;
}



/**
 * @brief Sets all secretion rates to zero
 * @details In cancer research, this models several scenarios:
 * - Drug treatments that inhibit secretory pathways
 * - Hypoxic conditions that may temporarily suppress non-essential secretion
 * - Quiescent cancer states with minimal microenvironmental interaction
 * - Experimental conditions where secretion is blocked to study tumor behavior
 * 
 * This method affects both direct secretion rates and net export rates.
 */
void Secrecion::set_todas_las_secreciones_a_cero( void )
{
	for( unsigned int i=0; i < tasas_de_secrecion.size(); i++ )
	{
		tasas_de_secrecion[i] = 0.0;
		tasas_de_exportacion_neta[i] = 0.0;
	}
	return;
}


/**
 * @brief Sets all consumption rates to zero
 * @details In cancer research, this models several scenarios:
 * - Metabolic inhibitor treatments that block consumption pathways
 * - Cell death processes where cells cease to consume resources
 * - Dormant cancer states with minimal metabolic activity
 * - Experimental conditions testing the effects of altered metabolism
 * 
 * Unlike normal cells, cancer cells often have high consumption rates for
 * glucose (Warburg effect) and glutamine, making consumption an important
 * therapeutic target.
 */
void Secrecion::set_todos_los_consumos_a_cero( void )
{
	for( unsigned int i=0; i < tasas_de_consumo.size(); i++ )
	{ tasas_de_consumo[i] = 0.0; }
	return;
}


/**
 * @brief Scales all secretion rates by a multiplicative factor
 * @param factor The multiplication factor to apply to secretion rates
 * @details This method models physiological and pathological changes in cancer secretion:
 * - Hypoxia-induced upregulation of angiogenic factors (factor > 1)
 * - Drug-induced suppression of secretory pathways (factor < 1)
 * - Phenotypic changes during disease progression (variable factor)
 * - Circadian rhythm effects on secretion (sinusoidal factor variation)
 * 
 * In cancer research, these changes in secretion are critical for understanding
 * how tumors interact with their microenvironment during different stages and
 * under different conditions.
 */
void Secrecion::multiplicar_las_secreciones_por_un_factor( double factor )
{
	for( unsigned int i=0; i < tasas_de_secrecion.size(); i++ )
	{
		tasas_de_secrecion[i] *= factor;
		tasas_de_exportacion_neta[i] *= factor;
	}
	return;
}


/**
 * @brief Scales all consumption rates by a multiplicative factor
 * @param factor The multiplication factor to apply to consumption rates
 * @details This method models physiological and pathological changes in cancer metabolism:
 * - Increased glycolysis in aggressive cancers (factor > 1 for glucose)
 * - Metabolic inhibitor effects (factor < 1)
 * - Adaptation to nutrient scarcity (variable factor)
 * - Changes in metabolic demands during different cell cycle phases
 * 
 * Cancer metabolism is a key research area and therapeutic target, as
 * cancer cells often exhibit altered consumption patterns compared to normal cells,
 * particularly through the Warburg effect (aerobic glycolysis) and glutamine addiction.
 */
void Secrecion::multiplicar_los_consumos_por_un_factor( double factor )
{
	for( unsigned int i=0; i < tasas_de_consumo.size(); i++ )
	{ tasas_de_consumo[i] *= factor; }
	return;
}

