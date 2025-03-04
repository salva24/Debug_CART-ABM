/**
 * @file Volumen.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Cell volume class that models volumetric properties of cancer cells
 *
 * @details This class defines the volume components of a cell, including nuclear and cytoplasmic
 * volumes, fluid and solid fractions, and target values for growth and division.
 * The volume properties are particularly relevant for cancer research, as cancer cells
 * often exhibit altered volume regulation, nuclear-cytoplasmic ratios, and fluid content.
 * 
 * @inbound_dependencies math.h
 *
 * @outbound_dependencies Used by Geometria.h, Celula.h, Ciclo_Modelo.h
 * 
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __VOLUMEN_H__
#define __VOLUMEN_H__

#include <math.h> 

/**
 * @class Volumen
 * @brief Models the volumetric properties of a cell with specific cancer research relevance
 * 
 * @details The Volumen class represents the various volume components of a cell, which
 * are critical factors in cancer biology. Cancer cells typically show altered volume
 * regulation, increased nuclear-cytoplasmic ratios, and changes in fluid content
 * compared to normal cells. This class models these properties to support realistic
 * cancer cell simulation, with default values based on MCF-7 breast cancer cells.
 * 
 * Cancer research applications:
 * - Nuclear enlargement is a key diagnostic criterion in cancer histopathology
 * - Nuclear-cytoplasmic ratio changes are associated with cancer progression
 * - Cellular volume affects drug uptake and response to therapy
 * - Volume regulation is often dysregulated in cancer cells
 * - Cell rupture due to excessive volume relates to certain cell death mechanisms
 */
class Volumen{
	public:

	    /**
	     * @brief Total cell volume in cubic microns
	     * @details Critical parameter for cancer cell modeling as tumor cells often show
	     * increased volume compared to normal cells of the same tissue type
	     */
	    double total;
	    
	    /**
	     * @brief Total solid (non-fluid) component of the cell volume
	     * @details Represents proteins, organelles, and structural components of the cell
	     */
	    double solido;
	    
	    /**
	     * @brief Total fluid component of the cell volume
	     * @details Cancer cells may have altered fluid content affecting drug diffusion
	     */
	    double fluido;
	    
	    /**
	     * @brief Fraction of cell volume that is fluid (0.0-1.0)
	     * @details Important for modeling drug diffusion and osmotic effects in cancer cells
	     */
	    double fraccion_de_fluido;
	    
	    /**
	     * @brief Nuclear volume in cubic microns
	     * @details Nuclear enlargement is a hallmark of many cancer types and used in diagnosis
	     */
	    double nuclear;
	    
	    /**
	     * @brief Fluid component of the nuclear volume
	     * @details Altered nuclear fluid content may affect chromatin organization and gene expression
	     */
	    double nuclear_fluido;
	    
	    /**
	     * @brief Solid component of the nuclear volume
	     * @details Represents chromatin, nucleoli and other nuclear structures, often increased in cancer
	     */
	    double nuclear_solido;
	    
	    /**
	     * @brief Cytoplasmic volume in cubic microns
	     * @details The cytoplasmic volume is total volume minus nuclear volume
	     */
	    double citoplasmatico;
	    
	    /**
	     * @brief Fluid component of the cytoplasmic volume
	     * @details Cytoplasmic fluid content affects molecular diffusion rates within cancer cells
	     */
	    double citoplasmatico_fluido;
	    
	    /**
	     * @brief Solid component of the cytoplasmic volume
	     * @details Represents cytoskeletal elements, organelles, and proteins
	     */
	    double citoplasmatico_solido;
	    
	    /**
	     * @brief Fraction of cell volume that is calcified (0.0-1.0)
	     * @details Models calcifications that can occur in certain cancer types (e.g., breast cancer)
	     */
	    double fraccion_calcificada;
	    
	    /**
	     * @brief Ratio of cytoplasmic to nuclear volume
	     * @details A key diagnostic parameter in cancer pathology; high N/C ratio often indicates malignancy
	     */
	    double relacion_citoplasma_nucleo;
	    
	    /**
	     * @brief Volume threshold at which the cell will rupture
	     * @details Important for modeling certain forms of cell death like necrosis in tumors
	     */
	    double volumen_de_ruptura;
	    
	    /**
	     * @brief Rate of change of cytoplasmic volume
	     * @details Models cytoplasmic growth rate, often dysregulated in cancer cells
	     */
	    double citoplasma_tasa_de_cambio;
	    
	    /**
	     * @brief Rate of change of nuclear volume
	     * @details Models nuclear growth rate, which may be accelerated in cancer cells
	     */
	    double nucleo_tasa_de_cambio;
	    
	    /**
	     * @brief Rate of change of fluid volume
	     * @details Models fluid uptake rate, affecting cell swelling in tumor microenvironments
	     */
	    double fluido_tasa_de_cambio;
	    
	    /**
	     * @brief Rate of calcification
	     * @details Models formation of microcalcifications seen in certain cancer types
	     */
	    double tasa_de_calcificacion;
	    
	    /**
	     * @brief Target value for solid cytoplasmic component
	     * @details Used to model directed growth toward a specific cytoplasmic composition
	     */
	    double target_citoplasma_solido;
	    
	    /**
	     * @brief Target value for solid nuclear component
	     * @details Used to model directed growth toward a specific nuclear composition
	     */
	    double target_nucleo_solido;
	    
	    /**
	     * @brief Target fluid fraction
	     * @details Models the cell's regulatory target for fluid content, often altered in tumors
	     */
	    double target_fraccion_fluido;
	    
	    /**
	     * @brief Target cytoplasm-to-nucleus ratio
	     * @details Models the cell's regulatory target for cytoplasm/nucleus balance
	     */
	    double target_relacion_citoplasma_nucleo;
	    
	    /**
	     * @brief Relative volume at which rupture occurs (multiple of initial volume)
	     * @details Used to calculate absolute rupture volume based on initial cell size
	     */
	    double volumen_de_ruptura_relativo;
	    
	    /**
	     * @brief Flag indicating if volume has changed
	     * @details Used to trigger recalculation of dependent parameters like cell geometry
	     */
	    bool cambio_el_volumen;
	    
	    /**
	     * @brief Constructor that initializes volume parameters
	     * @details Sets default volume values based on MCF-7 breast cancer cells
	     */
	    Volumen();
	    
	    /**
	     * @brief Divides all volume components in half
	     * @details Models volume redistribution during cell division
	     */
	    void dividir(void);
	    
	    /**
	     * @brief Multiplies all volume components by a factor
	     * @details Used to model uniform scaling of cell volume, such as during growth or shrinkage
	     */
	    void multiplicar(double);
};

#endif
