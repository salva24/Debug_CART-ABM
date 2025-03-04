/**
 * @file Geometria.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Header file defining the cell geometry class for cancer simulations
 *
 * @details This file defines the Geometria class, which manages geometric properties
 *          of cells such as radius, nuclear radius, surface area, and polarity.
 *          These properties are essential for accurately representing cancer cell
 *          morphology and physical interactions in simulations.
 * 
 * Inbound dependencies Volumen.h
 *
 * Outbound dependencies Used by Celula.h
 * 
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#include "Volumen.h"

/**
 * @brief Class for managing cell geometric properties in cancer simulations
 * @details The Geometria class calculates and maintains geometric properties of cells
 *          based on their volume. It handles the conversion between volume measurements
 *          and geometric properties like radius and surface area, which are critical for
 *          physical interactions between cells in the simulation.
 * 
 * @cancer_research Cell geometry is a fundamental aspect of cancer cell behavior and interactions.
 *                  Cancer cells often exhibit altered morphology compared to normal cells, including
 *                  changes in cell size, nuclear size, and nuclear-to-cytoplasmic ratio. These
 *                  alterations affect how cells interact with each other and respond to their
 *                  microenvironment, influencing tumor growth, invasion, and metastasis.
 */
class Geometria{
	public:
	/**
	 * @brief Cell radius in micrometers
	 * @details Represents the effective radius of the entire cell, calculated from total volume
	 * @cancer_research Changes in cell size are observed during cancer progression and in response to treatments
	 */
	double radio;
	
	/**
	 * @brief Nuclear radius in micrometers
	 * @details Represents the effective radius of the cell nucleus, calculated from nuclear volume
	 * @cancer_research Nuclear enlargement is a common feature in cancer cells and is used in cancer diagnosis
	 */
	double radio_nuclear;
	
	/**
	 * @brief Cell surface area in square micrometers
	 * @details Represents the total surface area of the cell, calculated from volume and radius
	 * @cancer_research Surface area affects cellular interactions, drug uptake, and mechanical properties
	 */
	double area_superficial;
	
	/**
	 * @brief Cell polarity factor (0.0-1.0)
	 * @details Represents the degree of polarization in the cell, which can affect migration direction
	 * @cancer_research Loss of cell polarity is a hallmark of epithelial cancers and promotes invasiveness
	 */
	double polaridad;

	/**
	 * @brief Constructor with default values based on MCF-7 breast cancer cells
	 * @cancer_research Initializes with empirically derived values from breast cancer cell measurements
	 */
	Geometria(); // done
	
	/**
	 * @brief Updates the cell radius based on total volume
	 * @param volumen Reference to the cell's volume properties
	 * @cancer_research Maintains accurate cell size for tumor growth and interaction modeling
	 */
	void actualizar_radio(Volumen& volumen); // done
	
	/**
	 * @brief Updates the nuclear radius based on nuclear volume
	 * @param volumen Reference to the cell's volume properties
	 * @cancer_research Nuclear size changes are important in cancer progression and diagnosis
	 */
	void actualizar_radio_nuclear(Volumen& volumen); // done
	
	/**
	 * @brief Updates the cell surface area based on total volume
	 * @param volumen Reference to the cell's volume properties
	 * @cancer_research Surface area affects interactions with drugs, immune cells, and other cells
	 */
	void actualizar_area_superficial(Volumen& volumen); // done
	
	/**
	 * @brief Updates all geometric properties based on volume
	 * @param volumen Reference to the cell's volume properties
	 * @cancer_research Ensures consistent cell morphology throughout simulation of tumor development
	 */
	void actualizar(Volumen& volumen); // done
};

