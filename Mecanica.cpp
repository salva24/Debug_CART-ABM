/**
 * @file Mecanica.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of cell mechanical properties for cancer simulation
 *
 * @details This file implements the Mecanica (Mechanics) class, which provides the
 * mechanical interaction properties used by cells in the cancer simulation. These
 * properties govern physical interactions between cells, with the extracellular matrix,
 * and with domain boundaries.
 * 
 * The default parameter values are calibrated based on reference data for cancer cells,
 * with specific emphasis on the balance of adhesive and repulsive forces that regulate
 * tumor organization, growth patterns, and invasion potential.
 * 
 * Inbound dependencies Mecanica.h
 *
 * Outbound dependencies None
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#include "Mecanica.h"

/**
 * @brief Constructor implementing default mechanical properties for cancer cells
 * @details Initializes mechanical properties with values calibrated for cancer simulation:
 * 
 * - Cell-cell interactions:
 *   - Adhesion (0.4): Moderate value reflecting reduced cell-cell adhesion in cancer cells
 *     compared to normal epithelial cells, modeling partial loss of contact inhibition
 *   - Repulsion (10.0): Strong repulsive force preventing excessive cell overlap,
 *     ensuring physical realism in tumor density modeling
 * 
 * - Cell-ECM interactions:
 *   - Adhesion (0.0): Default zero value allows for cell type-specific ECM adhesion
 *     to be defined separately, reflecting variable integrin expression in cancer
 *   - Repulsion (10.0): Strong repulsive force preventing cell penetration into ECM
 *     objects, important for modeling physical barriers to invasion
 * 
 * - Cell-boundary interactions:
 *   - Adhesion (4.0): Strong boundary adhesion useful for modeling attachment to
 *     basement membranes or blood vessel walls, critical in cancer invasion and metastasis
 *   - Repulsion (10.0): Strong repulsive force preventing boundary penetration,
 *     modeling anatomical barriers that constrain tumor growth
 * 
 * - Maximum adhesion distance (1.25): Set to 1.25 times the cell radius, allowing
 *   adhesive interactions to occur before direct contact, modeling extended reach
 *   through cell surface molecules and cellular protrusions
 * 
 * These default parameters produce realistic tumor spheroid formation while allowing
 * for invasion in the right conditions, balancing cohesion with invasive potential.
 * 
 * @return void
 */
Mecanica::Mecanica()
{
	fuerza_de_adhesion_cc = 0.4;
	fuerza_de_repulsion_cc = 10.0;
    fuerza_de_adhesion_co = 0.0;
    fuerza_de_repulsion_co = 10.0;
 	fuerza_de_adhesion_mb = 4.0;
	fuerza_de_repulsion_mb = 10.0;

	distancia_de_adhesion_maxima_relativa = 1.25; //Relativa a la célula. La distancia total es 1.25*Rcelula1 + 1.25*Rcelula2

	return;
}
