/**
 * @file: Tejido.h
 *
 * @author: Luciana Melina Luque
 *
 * @details
 * Defines the tissue structure that organizes cells and the microenvironment.
 * This file implements the Tejido (Tissue) class, which serves as a container and coordinator
 * for cells and their environment, enabling organized simulation of tissue-level phenomena.
 * 
 * Inbound Dependencies:
 * - Random.h - Stochastic behavior generation
 * - Microambiente.h - Environment for cell interactions
 * - Contenedor_de_Celulas.h - Cell container for management
 * - Parametros_globales.h - Global simulation parameters
 * 
 * Outbound Dependencies:
 * - Main.cpp - Uses tissue for simulation
 * 
 * Usage:
 * Tejido tejido;
 * tejido.inicializar_tejido();
 * tejido.geometria_del_tumor();
 * tejido.introducir_linfocitos(n);
 * 
 * @license: Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#ifndef __TEJIDO_H__
#define __TEJIDO_H__

#include "Random.h"
#include "Microambiente.h"
#include "Contenedor_de_Celulas.h"
#include "Parametros_globales.h"

extern RNG *rng;
extern Parametros_globales *pg;
extern std::vector<Celula*> todas_las_celulas;
extern std::vector<Celula*> celulas_listas_para_dividirse;
extern std::vector<Celula*> celulas_para_registrar_en_voxeles;


/**
 * @class Tejido
 * @brief Models a tissue structure containing cells and microenvironment
 * 
 * The Tejido (Tissue) class serves as the primary organizational structure in the simulation,
 * containing both the microenvironment and the cell container. It provides methods for
 * initializing tissue structure, creating spatial arrangements of cells, introducing
 * immune cells, and calculating tumor geometry.
 * 
 * Cancer Research Context:
 * In cancer research, the tissue organization is critical for studying tumor growth,
 * infiltration, and interaction with surrounding healthy tissue. This class enables
 * realistic spatial modeling of tumor development within host tissue.
 */
class Tejido{
	public:

	Microambiente microambiente;
	Contenedor_de_Celulas cdc;
    double radio_del_tumor;
    double volumen_del_tumor;
    double volumen_del_tumor2;
    int celulas_tumorales;
    int celulas_muertas;

    /**
     * @brief Default constructor for the Tejido class
     */
	Tejido();

    /**
     * @brief Initializes the tissue structure with cells and microenvironment
     * 
     * Sets up the microenvironment, cell container, and creates the initial
     * spherical arrangement of cells for the simulation.
     */
	void inicializar_tejido();
    
    /**
     * @brief Creates planes of healthy cells along the Z-axis
     * 
     * @param cantidad Number of planes to create
     * @param posicion_en_z_del_primer_plano Z-position of the first plane
     * @return Vector of cell positions
     */
	std::vector<Vector> crear_planos_de_celulas_sanas_en_Z(int cantidad, double posicion_en_z_del_primer_plano);
    
    /**
     * @brief Creates a spherical arrangement of cells
     * 
     * @param radio_de_la_esfera Radius of the sphere
     * @return Vector of cell positions
     */
	std::vector<Vector> crear_esfera_de_celulas(double radio_de_la_esfera);
    
    /**
     * @brief Introduces lymphocytes around the tumor periphery
     * 
     * Adds immune cells (lymphocytes) in a shell around the tumor at a
     * specified distance from the tumor edge.
     * 
     * @param cantidad Number of lymphocytes to introduce
     */
	void introducir_linfocitos(int cantidad);
    
    /**
     * @brief Introduces lymphocytes randomly in the domain
     * 
     * Adds immune cells (lymphocytes) at random positions in the simulation
     * domain, outside a minimum distance from the tumor.
     * 
     * @param cantidad Number of lymphocytes to introduce
     */
    void introducir_linfocitos_aleatorios(int cantidad);
    
    /**
     * @brief Calculates geometric properties of the tumor
     * 
     * Computes the tumor radius, volume, cell count, and other geometric
     * metrics to track tumor progression over time.
     */
    void geometria_del_tumor();
};


#endif

