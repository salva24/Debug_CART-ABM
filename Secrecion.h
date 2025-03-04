/**
 * @file Secrecion.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Cell secretion and consumption model for cancer cell-microenvironment interactions
 *
 * @details This file defines the Secrecion class that models how cancer cells secrete and
 * consume various substances in the tumor microenvironment, including:
 * - Extracellular signaling molecules (growth factors, cytokines)
 * - Metabolites (lactate, waste products)
 * - Matrix-degrading enzymes (MMPs)
 * - Oncoproteins that may affect neighboring cells
 * 
 * @inbound_dependencies <vector>
 *
 * @outbound_dependencies Celula.h
 *
 * @usage Used by Celula class to model secretion/consumption behaviors
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __SECRECION_H__
#define __SECRECION_H__


#include <vector>

/**
 * @brief Models how cancer cells secrete substances into and consume substances from their microenvironment
 * @details The Secrecion class encapsulates the secretion and consumption behaviors of 
 * cancer cells, which are critical for:
 * - Tumor microenvironment remodeling (acidification, ECM degradation)
 * - Autocrine/paracrine signaling to promote proliferation and survival
 * - Angiogenesis induction through secretion of pro-angiogenic factors
 * - Immunosuppression through secretion of immunomodulatory factors
 * - Metabolic adaptation through altered consumption patterns
 * 
 * In cancer research, altered secretion profiles are characteristic of various cancer types
 * and can serve as biomarkers, therapeutic targets, and indicators of disease progression.
 */
class Secrecion
{
 private:
 public:

	/**
	 * @brief Secretion rates for each substrate in the microenvironment
	 * @details In cancer cells, secretion rates for certain factors are often upregulated,
	 * including VEGF (promotes angiogenesis), MMPs (promote invasion), and TGF-β 
	 * (promotes invasion in advanced stages). Units: substrate/min
	 */
	std::vector<double> tasas_de_secrecion;
	
	/**
	 * @brief Consumption rates for each substrate in the microenvironment
	 * @details Cancer cells typically show altered consumption patterns, particularly
	 * increased glucose consumption (Warburg effect) and glutamine addiction.
	 * Units: substrate/min
	 */
	std::vector<double> tasas_de_consumo;
	
	/**
	 * @brief Saturation densities limiting secretion rates
	 * @details Represents negative feedback mechanisms where high environmental concentrations
	 * inhibit further secretion. Important for modeling homeostatic-like behavior in tumors.
	 * Units: substrate/volume
	 */
	std::vector<double> densidades_de_saturacion;
	
	/**
	 * @brief Net export rates (secretion minus uptake)
	 * @details Net impact of the cell on each substrate's concentration in the microenvironment.
	 * Positive values indicate net addition to the environment, negative values indicate net removal.
	 * Units: substrate/min
	 */
	std::vector<double> tasas_de_exportacion_neta;
	
	/**
	 * @brief Oncoprotein concentration
	 * @details Models the production of oncogenic proteins that may be exported to affect
	 * neighboring cells. Examples include mutant p53 that can be transmitted via exosomes
	 * or oncogenic WNT signaling factors. Units: concentration
	 */
	double oncoproteina;

	/**
	 * @brief Default constructor
	 * @details Initializes empty substrate vectors and sets oncoprotein level to zero.
	 * These will be populated based on the microenvironment configuration.
	 */
	Secrecion();

	/**
	 * @brief Advances the secretion and consumption processes by time step dt
	 * @param dt Time step duration
	 * @details Updates the microenvironment based on the cell's current secretion and
	 * consumption rates. Critical for modeling how cancer cells modify their environment
	 * over time, creating gradients that affect nearby cells.
	 */
	void avanzar( double dt );

	/**
	 * @brief Sets all secretion rates to zero
	 * @details Used when modeling cancer treatments that block secretion pathways
	 * or when modeling quiescent cancer cells with reduced metabolic activity.
	 */
	void set_todas_las_secreciones_a_cero( void );
	
	/**
	 * @brief Sets all consumption rates to zero
	 * @details Used when modeling cancer treatments that block metabolic pathways
	 * or when modeling cell death where consumption ceases.
	 */
	void set_todos_los_consumos_a_cero( void );
	
	/**
	 * @brief Scales all secretion rates by a factor
	 * @param factor Multiplication factor for secretion rates
	 * @details Used to model altered secretion in different cancer stages or
	 * in response to treatments. For example, hypoxia can increase VEGF secretion 
	 * by an order of magnitude.
	 */
	void multiplicar_las_secreciones_por_un_factor( double factor );
	
	/**
	 * @brief Scales all consumption rates by a factor
	 * @param factor Multiplication factor for consumption rates
	 * @details Used to model altered metabolism in different cancer stages or
	 * in response to treatments. For example, some therapies target the increased
	 * metabolic demands of cancer cells.
	 */
	void multiplicar_los_consumos_por_un_factor( double factor );
};

#endif
