/**
 * @file Motilidad.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Cell motility model for cancer cell migration and invasion
 *
 * @details This file defines the Motilidad (Motility) class that models how cancer
 * cells move through their microenvironment. Cell motility is a critical aspect of:
 * - Cancer invasion into surrounding tissues
 * - Metastatic spread to distant organs
 * - Immune cell infiltration into tumors
 * - Cellular responses to microenvironmental gradients
 * 
 * The model supports random migration, directional bias, and chemotaxis, allowing
 * simulation of various motility patterns observed in different cancer types.
 * 
 * Inbound dependencies "Vector.h", <vector>
 *
 * Outbound dependencies Celula.h
 *
 * @usage Used by Celula class to model cell movement behaviors
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __MOTILIDAD_H__
#define __MOTILIDAD_H__


#include "Vector.h"
#include<vector>

extern Vector *v;

/**
 * @brief Models how cancer cells move through their microenvironment
 * @details The Motilidad class encapsulates cell movement behaviors critical for
 * cancer progression, particularly invasion and metastasis. It includes:
 * - Basic motility properties (speed, persistence)
 * - Directional bias capabilities (modeling haptotaxis/durotaxis)
 * - Chemotaxis properties (response to chemical gradients)
 * 
 * In cancer research, altered motility is associated with:
 * - Epithelial-to-mesenchymal transition (EMT)
 * - Metastatic potential
 * - Therapeutic resistance through cell escape
 * - Collective vs. individual cell migration modes
 */
class Motilidad{
 public:
 	/**
	 * @brief Indicates whether the cell can move
	 * @details In cancer, motility activation is often a key step in progression
	 * from in situ to invasive disease. This boolean flag allows modeling of:
	 * - Motile invasive cancer cells vs. stationary cancer in situ
	 * - EMT-induced motility switches
	 * - Treatment-induced changes in motility status
	 */
	bool es_movil;

	/**
	 * @brief Time duration for which a cell maintains its current direction
	 * @details Persistence time affects the pattern of cell migration:
	 * - Longer times create more ballistic (straight-line) movement
	 * - Shorter times create more diffusive (random) movement
	 * - Different cancer types show characteristic persistence times
	 * Units: minutes
	 */
	double tiempo_de_persistencia;
	
	/**
	 * @brief Cell migration speed
	 * @details The speed at which cells migrate varies significantly across cancer types:
	 * - Melanoma: Often highly motile (0.1-1.0 μm/min)
	 * - Glioblastoma: Rapidly invasive in brain tissue (~0.5 μm/min)
	 * - Breast cancer: Variable speed depending on subtype
	 * Units: μm/minute
	 */
	double velocidad_de_migracion;
	
	/**
	 * @brief Direction vector for biased migration
	 * @details Cancer cells often move preferentially in certain directions due to:
	 * - Haptotaxis (movement along adhesion molecule gradients)
	 * - Durotaxis (movement along stiffness gradients)
	 * - Mechanical cues from surrounding ECM
	 * - Pre-existing anatomical structures ("highways" for invasion)
	 */
	Vector bias_de_la_migracion_direccion;
	
	/**
	 * @brief Strength of directional bias
	 * @details Controls how strongly the bias direction influences cell movement:
	 * - 0.0: No bias (pure random motion)
	 * - 1.0: Maximum bias (deterministic motion in bias direction)
	 * This allows modeling of cancer cells' varying responses to environmental cues.
	 */
	double bias_de_la_migracion;
	
	/**
	 * @brief Current motility vector
	 * @details The actual direction of cell movement, resulting from:
	 * - Random persistence
	 * - Directional bias
	 * - Chemotactic response
	 * This vector drives the physical displacement of the cell.
	 */
	Vector vector_de_motilidad;

	/**
	 * @brief Index of the chemotactic substrate the cell responds to
	 * @details Cancer cells can migrate along gradients of various factors:
	 * - Growth factors (EGF, HGF) driving migration toward blood vessels
	 * - Oxygen gradients (hypoxia-driven migration)
	 * - Nutrient gradients
	 * This index points to the specific substrate in the microenvironment.
	 */
	int quimiotaxis_indice;
	
	/**
	 * @brief Direction of chemotactic response
	 * @details Determines whether cells move up or down the gradient:
	 * - 1: Positive chemotaxis (toward higher concentration)
	 * - -1: Negative chemotaxis (toward lower concentration)
	 * Cancer cells often show positive chemotaxis to growth factors and
	 * negative chemotaxis away from hypoxic regions.
	 */
	int quimiotaxis_direccion;

	/**
	 * @brief Default constructor
	 * @details Initializes a cell with default motility properties:
	 * - Non-motile by default (es_movil = false)
	 * - Standard persistence time (1.0 minute)
	 * - Standard migration speed (1.0 μm/minute) when active
	 * - No initial bias or motility direction
	 * - Set for positive chemotaxis to first substrate by default
	 */
	Motilidad();


};

#endif
