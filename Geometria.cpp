/**
 * @file Geometria.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of cell geometry functions for cancer simulations
 *
 * @details This file implements the Geometria class, which manages the geometric properties
 *          of cells, including radius, nuclear radius, and surface area calculations.
 *          These properties are critical for accurately representing cancer cell morphology
 *          and physical interactions in the simulation.
 * 
 * Inbound dependencies math.h, Geometria.h, Volumen.h (indirect)
 *
 * Outbound dependencies None
 * 
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#include <math.h>
#include "Geometria.h"

/**
 * @brief Constructor for Geometria class with default values based on MCF-7 breast cancer cells
 * @details Initializes cell geometry parameters with values derived from experimental measurements
 *          of MCF-7 breast cancer cells. These reference values provide a biologically realistic
 *          starting point for cancer cell simulations.
 * 
 * @cancer_research Used to establish baseline morphological characteristics for simulated cancer cells,
 *                  allowing for realistic representation of cell size and shape in tumor models.
 */
Geometria::Geometria(){




	radio = 8.412710547954228;
	radio_nuclear = 5.051670902881889;
	area_superficial = 889.3685284131693;

	polaridad = 0.0;
	return;
}

/**
 * @brief Updates the cell radius based on its total volume
 * @details Calculates the radius of a cell assuming a spherical shape, using the formula
 *          r = (V / (4π/3))^(1/3) where V is the total cell volume.
 * 
 * @param volumen Reference to the cell's volume properties
 * 
 * @cancer_research Cell radius is a critical parameter for modeling physical interactions between
 *                  cancer cells and determining spatial organization within tumors. Changes in cell
 *                  volume and radius are associated with cancer progression and response to therapy.
 */
void Geometria::actualizar_radio(Volumen& volumen){

	static double cuatro_tercios_de_pi = 4.188790204786391;
	radio = volumen.total;
	radio /= cuatro_tercios_de_pi;
	radio = pow( radio, 0.333333333333333333333333333333333333333);
	return;

}

/**
 * @brief Updates the nuclear radius based on nuclear volume
 * @details Calculates the radius of the cell nucleus assuming a spherical shape, using the
 *          formula r = (V / (4π/3))^(1/3) where V is the nuclear volume.
 * 
 * @param volumen Reference to the cell's volume properties
 * 
 * @cancer_research Nuclear radius is an important indicator in cancer diagnosis and progression.
 *                  Cancer cells often show altered nuclear morphology with increased nuclear size,
 *                  irregular shape, and altered nuclear-to-cytoplasmic ratio, all of which are
 *                  captured by this model.
 */
void Geometria::actualizar_radio_nuclear(Volumen& volumen){

	static double cuatro_tercios_de_pi = 4.188790204786391;
	radio_nuclear = volumen.nuclear;
	radio_nuclear /= cuatro_tercios_de_pi;
	radio_nuclear = pow( radio_nuclear, 0.333333333333333333333333333333333333333);
	return;

}

/**
 * @brief Updates the cell surface area based on total volume
 * @details Calculates the surface area of a cell assuming a spherical shape, using a derived
 *          constant that relates volume to surface area: SA = V^(2/3) / constant
 * 
 * @param volumen Reference to the cell's volume properties
 * 
 * @cancer_research Surface area affects cellular interactions, drug uptake, and mechanical properties.
 *                  Cancer cells often have altered surface area due to morphological changes, which
 *                  impacts their ability to interact with the microenvironment, migrate, and respond
 *                  to treatments.
 */
void Geometria::actualizar_area_superficial(Volumen& volumen){

	// 4pi / (4pi/3)^(2/3)
	static double constante = 4.835975862049409;
	area_superficial = pow( volumen.total , 0.666666666666667 );
	area_superficial /= constante;

	return;

}

/**
 * @brief Updates all geometric properties of the cell based on its volume
 * @details Calls individual update methods for radius and nuclear radius, then calculates
 *          surface area using a different method based on the relationship between volume and radius:
 *          SA = 3V/r, which simplifies to 4πr² for a sphere.
 * 
 * @param volumen Reference to the cell's volume properties
 * 
 * @cancer_research Complete geometric updating is essential for modeling how cancer cells change
 *                  morphologically during growth, division, and in response to treatments. This method
 *                  ensures all geometric properties remain consistent as the cell volume changes,
 *                  providing realistic simulations of tumor growth dynamics and cellular interactions.
 */
void Geometria::actualizar(Volumen& volumen){

	actualizar_radio(volumen);
	actualizar_radio_nuclear(volumen);

	// surface area = 4*pi*r^2 = (4/3)*pi*r^3 / (r/3)
	area_superficial = volumen.total;
	area_superficial /= radio;
	area_superficial *= 3.0;

	return;
}

