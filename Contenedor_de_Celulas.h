/**
 * @file: Contenedor_de_Celulas.h
 *
 * @author: Luciana Melina Luque
 *
 * @details Defines the cell container system that manages collections of cells in the simulation.
 * This class handles cell organization, interaction, mechanics, and lifecycle events (division, death)
 * within the spatial domain of the simulation.
 * 
 * Inbound Dependencies:
 * - Grillado.h - Cartesian grid for spatial organization
 * - Celula.h - Cell class definition
 * - vector - Standard container for cell collections
 * - fstream - I/O operations
 * 
 * Outbound Dependencies:
 * - Tejido.h - Uses cell container for tissue organization
 * 
 * Usage:
 * Contenedor_de_Celulas cdc;
 * cdc.inicializar(x1, x2, y1, y2, z1, z2, dx, dy, dz);
 * cdc.registrar_celula(celula);
 * cdc.actualizar_todas_las_celulas(tiempo, dt_diff, dt_mec, dt_ciclo);
 * 
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#ifndef __CONTENEDOR_DE_CELULAS_H__
#define __CONTENEDOR_DE_CELULAS_H__

#include "Grillado.h"
#include "Celula.h"
#include <vector>
#include <fstream>


/**
 * @class Contenedor_de_Celulas
 * @brief Manages collections of cells within the spatial domain
 * 
 * The cell container class organizes cells in a spatial grid (voxels),
 * handles cell interactions, updates cell states, and manages mechanical
 * forces between cells. It also coordinates cell lifecycle events like
 * division and death.
 * 
 * Cancer Research Context:
 * In cancer modeling, efficient cell organization and interaction handling
 * is critical for simulating tumor growth, cell-cell interactions, and
 * spatial organization. This container enables scaling to large numbers
 * of cells while maintaining computational efficiency.
 */
class Contenedor_de_Celulas{
	private:

	public:
		Grillado_Cartesiano grillado;
		std::vector<std::vector<Celula*> > celulas_en_voxel;
		std::vector<std::vector<Celula*> > celulas_fuera_del_domino;
		int num_de_divisiones_en_este_paso;
		int num_de_muertes_en_este_paso;
		int num_de_celulas;
        double tiempo_desde_la_ultima_mecanica;
        double hora_de_la_ultima_mecanica;
		Celula* celula;
		
		/**
		 * @brief Default constructor
		 * 
		 * Initializes a new cell container with default values.
		 */
		Contenedor_de_Celulas();
		
		/**
		 * @brief Initializes the spatial domain of the cell container
		 * 
		 * Sets up the Cartesian grid with the specified dimensions and
		 * creates the data structures needed to track cells within voxels.
		 * 
		 * @param x_ini Minimum x-coordinate
		 * @param x_fin Maximum x-coordinate
		 * @param y_ini Minimum y-coordinate
		 * @param y_fin Maximum y-coordinate
		 * @param z_ini Minimum z-coordinate
		 * @param z_fin Maximum z-coordinate
		 * @param dx Size of voxel in x-direction
		 * @param dy Size of voxel in y-direction
		 * @param dz Size of voxel in z-direction
		 */
		void inicializar(double x_ini, double x_fin, double y_ini, double y_fin, double z_ini, double z_fin , double dx, double dy, double dz);
		
		/**
		 * @brief Registers a cell in the container
		 * 
		 * Places the cell in the appropriate voxel based on its position.
		 * If the position is invalid, the cell is marked as outside the domain.
		 * 
		 * @param celula Pointer to the cell to register
		 */
		void registrar_celula( Celula* celula );
		
		/**
		 * @brief Adds a cell to a specific voxel
		 * 
		 * @param celula Pointer to the cell to add
		 * @param indice_de_voxel Index of the voxel to add the cell to
		 */
		void agregar_celula_a_voxel(Celula* celula, int indice_de_voxel);
		
		/**
		 * @brief Updates all cells for a simulation step
		 * 
		 * This is a key method that orchestrates the complete update sequence:
		 * 1. Updates cell secretion and consumption
		 * 2. Calculates and applies mechanical forces between cells
		 * 3. Updates cell positions
		 * 4. Advances cell phenotypes (cycle, death, etc.)
		 * 5. Processes division of cells ready to divide
		 * 6. Registers new cells in voxels
		 * 7. Removes cells marked for death
		 * 
		 * @param tiempo_total Current simulation time
		 * @param dt_difusion Time step for diffusion
		 * @param dt_mecanico Time step for mechanics
		 * @param dt_ciclo Time step for cell cycle
		 */
		void actualizar_todas_las_celulas(double tiempo_total, double dt_difusion, double dt_mecanico, double dt_ciclo);
		
		/**
		 * @brief Checks if a voxel contains any cells
		 * 
		 * @param indice_de_voxel Index of the voxel to check
		 * @return true if the voxel contains at least one cell, false otherwise
		 */
		bool contiene_alguna_celula(int indice_de_voxel);
		
		/**
		 * @brief Removes a cell from a voxel
		 * 
		 * Used when cells move between voxels or when they die.
		 * 
		 * @param celula Pointer to the cell to remove
		 * @param indice_de_voxel Index of the voxel containing the cell
		 */
		void sacar_celula_de_voxel(Celula* celula, int indice_de_voxel);
		
		/**
		 * @brief Updates voxel assignments for all cells
		 * 
		 * Checks each cell's position and moves it to the correct voxel
		 * if it has moved to a different voxel. This is typically called
		 * periodically rather than after every position update for efficiency.
		 */
		void actualizar_voxeles_de_celulas();
		
		/**
		 * @brief Calculates mechanical forces between two cells
		 * 
		 * Computes repulsive and adhesive forces between cells based on
		 * their distance, size, and mechanical properties. Updates velocity
		 * vectors for both cells accordingly.
		 * 
		 * @param celula Pointer to the first cell
		 * @param otra_celula Pointer to the second cell
		 */
		void agregar_potenciales_cdc(Celula* celula, Celula* otra_celula);
};

/** Global vectors for tracking cells throughout the simulation */
extern std::vector<Celula*> todas_las_celulas;
/** Vector of cells ready to divide in the current step */
extern std::vector<Celula*> celulas_listas_para_dividirse;
/** Vector of cells to be removed in the current step */
extern std::vector<Celula*> celulas_listas_para_remover;
/** Vector of cells to be registered in voxels */
extern std::vector<Celula*> celulas_para_registrar_en_voxeles;

#endif
