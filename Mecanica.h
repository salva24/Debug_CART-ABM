/**
 * @file Mecanica.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Cell mechanics properties definition for cancer simulation
 *
 * @details This file defines the Mecanica (Mechanics) class, which encapsulates
 * the mechanical properties of cells that govern physical interactions with other cells,
 * the extracellular matrix, and boundary surfaces. These properties determine how
 * cells adhere to and repel each other, forming the basis for tissue architecture.
 * 
 * In cancer research, mechanical properties are critically important as they influence:
 * - Tumor compactness and organization
 * - Invasive capacity and metastatic potential
 * - Response to mechanical stresses in the tumor microenvironment
 * - Physical barriers to drug delivery
 * - Mechanical cues that can trigger phenotypic changes
 * 
 * Inbound dependencies None
 *
 * Outbound dependencies Fenotipo.h, Celula.h
 *
 * @usage Instantiated as a member variable within each Fenotipo object
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __MECANICA_H__
#define __MECANICA_H__

/**
 * @class Mecanica
 * @brief Cell mechanical properties for modeling physical interactions in cancer
 * 
 * @details The Mecanica class defines the parameters that govern how cells
 * physically interact with their environment, including other cells, the extracellular
 * matrix (ECM), and domain boundaries. These properties determine force generation
 * during cell-cell contact, cell-ECM adhesion, and boundary interactions.
 * 
 * In cancer research, altered mechanical properties are key characteristics that:
 * - Enable tumor cells to overcome contact inhibition
 * - Facilitate invasion through altered adhesion properties
 * - Create mechanical feedback that influences tumor growth patterns
 * - Affect tissue mechanics that can impact drug delivery and efficacy
 * - Contribute to metastatic potential through modified ECM interactions
 */
class Mecanica
{
 public:
    /**
     * @brief Repulsive force coefficient between cells
     * @details In cancer, this parameter affects tumor density and organization.
     * Higher values model more rigid cells that resist compression, while lower
     * values allow greater cell packing, affecting tumor compactness and interstitial
     * pressure, which impacts drug penetration and hypoxia.
     */
    double fuerza_de_repulsion_cc;
    
    /**
     * @brief Adhesive force coefficient between cells
     * @details Represents cell-cell adhesion strength, typically altered in cancer cells
     * through changes in cadherin expression. Reduced values model cancer cells with 
     * decreased cell-cell adhesion (a hallmark of EMT), associated with increased
     * invasiveness and metastatic potential.
     */
    double fuerza_de_adhesion_cc;
    
    /**
     * @brief Repulsive force coefficient between cells and ECM objects
     * @details Models cell stiffness response to ECM components. In cancer research,
     * this parameter affects how tumor cells navigate through the stroma and interact
     * with structural tissue components, influencing invasion patterns.
     */
    double fuerza_de_repulsion_co;
    
    /**
     * @brief Adhesive force coefficient between cells and ECM objects
     * @details Represents binding strength to ECM through integrins and other adhesion
     * molecules. Cancer cells often show altered ECM adhesion, with different patterns
     * enabling migration along ECM tracks or through tissue boundaries. Key parameter
     * for modeling invasion and metastasis.
     */
    double fuerza_de_adhesion_co;
    
    /**
     * @brief Repulsive force coefficient between cells and domain boundaries
     * @details Models cell response to physical boundaries, important for simulating
     * confined tumor growth or interactions with anatomical barriers. In cancer research,
     * this affects modeling of tumor growth against physical constraints.
     */
    double fuerza_de_repulsion_mb;
    
    /**
     * @brief Adhesive force coefficient between cells and domain boundaries
     * @details Models cell adhesion to physical boundaries, important for simulating
     * basement membrane attachment or blood vessel wall interactions. Particularly
     * relevant for modeling intravasation/extravasation during metastasis.
     */
    double fuerza_de_adhesion_mb;
    
    /**
     * @brief Maximum relative distance for adhesion interactions
     * @details Defines the maximum distance (as a multiple of cell radius) at which
     * adhesion forces are active. In cancer modeling, this parameter affects cell-cell
     * interaction ranges and the formation of clusters. Increased values can represent
     * extended reach through cellular protrusions, important in cancer cell migration.
     */
	double distancia_de_adhesion_maxima_relativa;


    /**
     * @brief Default constructor with cancer-relevant mechanical parameters
     * @details Initializes mechanical properties with default values calibrated
     * for cancer cell simulation. Different cancer types can be modeled by 
     * adjusting these baseline mechanical properties.
     */
	Mecanica();


};

#endif
