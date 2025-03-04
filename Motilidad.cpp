/**
 * @file Motilidad.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of cancer cell motility model
 *
 * @details This file implements the Motilidad class methods that model how cancer cells
 * migrate through their microenvironment. Cancer cell motility is a key factor in:
 * - Tumor invasion into surrounding tissue
 * - Metastatic dissemination
 * - Resistance to treatments through escape mechanisms
 * - Interaction with immune cells
 * 
 * The implementation provides default motility parameters based on common
 * cancer cell types studied in vitro and in vivo.
 * 
 * Inbound dependencies "Motilidad.h"
 *
 * Outbound dependencies None
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#include "Motilidad.h"

/**
 * @brief Constructor initializes default motility parameters
 * @details Creates a motility object with parameters similar to observed cancer cell behavior:
 * - Initially non-motile (es_movil = false)
 * - 1.0 minute persistence time (typical for epithelial cancer cells)
 * - 1.0 μm/minute migration speed (moderate migration rate, variable by cancer type)
 * - No directional bias (unbiased random migration)
 * - Zero motility vector (no initial movement)
 * - Chemotactic response to first substrate (index 0, often oxygen)
 * - Positive chemotaxis direction (cells move toward higher concentrations)
 * 
 * In cancer models, these parameters would be modified to reflect:
 * - Highly invasive phenotypes (higher speeds, longer persistence)
 * - EMT-induced motility (activation of movement)
 * - Gradient-directed invasion (chemotaxis toward nutrients/away from stress)
 * - Treatment effects on migration (reduced speed/persistence)
 */
Motilidad::Motilidad()
{
	es_movil = false;

	tiempo_de_persistencia = 1.0;
	velocidad_de_migracion = 1.0;

	bias_de_la_migracion_direccion.x=0.0;
	bias_de_la_migracion_direccion.y=0.0;
	bias_de_la_migracion_direccion.z=0.0;
	bias_de_la_migracion = 0.0;

	vector_de_motilidad.x=0.0;
	vector_de_motilidad.y=0.0;
	vector_de_motilidad.z=0.0;

	quimiotaxis_indice = 0;
	quimiotaxis_direccion = 1;

	return;
}
