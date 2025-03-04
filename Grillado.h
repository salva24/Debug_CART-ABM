/**
 * @file Grillado.h
 *
 * @author Luciana Melina Luque
 *
 * @details Implements the spatial discretization grid system used for the cancer simulation.
 * 
 * This file defines the Cartesian grid structure that provides spatial organization for the cancer
 * simulation microenvironment. It enables calculations of substrate diffusion, gradient formation,
 * and spatial organization of cells within the tumor microenvironment.
 * 
 * Inbound Dependencies: Vector.h, Voxel.h, Parametros_globales.h, Macros.h
 *
 * Outbound Dependencies: Microambiente.h, Contenedor_de_Celulas.h
 *
 * Usage: Used by the microenvironment to implement the finite volume method for diffusion and
 *        by the cell container for spatial organization of cells.
 *
 * @license: Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __GRILLADO_H__
#define __GRILLADO_H__

#include <iostream>
#include "Macros.h"
#include "Voxel.h"
#include "Parametros_globales.h"

extern Parametros_globales *pg;

/**
 * @class Grillado_General
 * @brief Base class for grid systems used in cancer microenvironment modeling.
 * 
 * This class implements a general grid system that divides the simulation domain into voxels
 * for efficient spatial organization. In cancer modeling, this spatial discretization is 
 * essential for accurately representing tumor heterogeneity and the formation of
 * microenvironmental gradients that influence cancer cell behavior.
 */
class Grillado_General
{
 private:
	friend std::ostream& operator<<(std::ostream& os, const Grillado_General& grillado);


 public:
	/**
	 * @brief Constructor for the Grillado_General class.
	 * 
	 * Initializes a default grid with a small spatial domain. In cancer modeling,
	 * appropriate grid initialization ensures proper resolution of tumor microenvironments.
	 */
	Grillado_General();

	/**
	 * @brief The bounding box coordinates of the grid.
	 * 
	 * Stores the minimum and maximum coordinates in each dimension.
	 * Format: [x_min, y_min, z_min, x_max, y_max, z_max]
	 * Essential for defining the spatial extent of the tumor microenvironment.
	 */
	std::vector<double> caja;
	
	/**
	 * @brief Vector containing all voxels in the grid.
	 * 
	 * Each voxel represents a discrete spatial region within the tumor microenvironment
	 * where biochemical substrates can diffuse and cancer cells can reside.
	 */
	std::vector<Voxel> voxeles;
	
	/**
	 * @brief Connectivity information between voxels.
	 * 
	 * Stores indices of connected voxels for each voxel in the grid.
	 * In cancer modeling, voxel connectivity is essential for simulating
	 * diffusion of signaling molecules, nutrients, and oxygen that influence
	 * cancer cell behavior and tumor growth.
	 */
	std::vector< std::vector<int> > indices_de_voxeles_conectados;

	/**
	 * @brief Finds the index of the voxel nearest to a given position.
	 * 
	 * @param posicion The 3D position to query.
	 * @return The index of the nearest voxel.
	 * 
	 * In cancer modeling, this function allows cells to interact with the
	 * appropriate microenvironmental conditions at their location, which influences
	 * cancer cell phenotypes, division rates, and death rates.
	 */
	int indice_del_voxel_mas_cercano( Vector& posicion );
	
	/**
	 * @brief Checks if a position is valid within the grid domain.
	 * 
	 * @param x X-coordinate to check.
	 * @param y Y-coordinate to check.
	 * @param z Z-coordinate to check.
	 * @return Boolean indicating if the position is valid.
	 * 
	 * Ensures cells remain within the defined tumor microenvironment boundaries
	 * during cancer growth simulations.
	 */
	bool es_valida_la_posicion(double x, double y, double z);

	/**
	 * @brief Connects two voxels in the grid.
	 * 
	 * @param i Index of the first voxel.
	 * @param j Index of the second voxel.
	 * 
	 * Sets up bidirectional connectivity between voxels, which is
	 * essential for substrate diffusion calculations in cancer microenvironments.
	 */
	void conectar_voxeles(int i,int j);
	
	/**
	 * @brief Units of measurement for the grid.
	 * 
	 * Typically set to "micrometers" for cancer modeling to match
	 * cellular and tissue scales.
	 */
	std::string unidades;
	
	/**
	 * @brief Displays general information about the grid.
	 * 
	 * @param os Output stream to write information to.
	 * 
	 * Provides information about the grid dimensions and properties
	 * that define the cancer simulation spatial domain.
	 */
	void mostrar_informacion_general( std::ostream& os);
};


/**
 * @class Grillado_Cartesiano
 * @brief Implements a Cartesian grid for cancer microenvironment modeling.
 * 
 * Extends the general grid with specific Cartesian organization.
 * In cancer modeling, Cartesian grids allow for efficient representation
 * of tumor tissues and the surrounding microenvironment, enabling
 * accurate simulation of substrate gradients (oxygen, nutrients, growth factors)
 * that influence cancer cell behavior and tumor evolution.
 */
class Grillado_Cartesiano : public Grillado_General
{
 private:

 public:
	/**
	 * @brief Coordinate values along the x-axis.
	 * 
	 * Defines the x-positions of voxel centers, which is important
	 * for mapping the spatial distribution of tumor cells.
	 */
	std::vector<double> coordenadas_x;
	
	/**
	 * @brief Coordinate values along the y-axis.
	 * 
	 * Defines the y-positions of voxel centers, which is important
	 * for mapping the spatial distribution of tumor cells.
	 */
	std::vector<double> coordenadas_y;
	
	/**
	 * @brief Coordinate values along the z-axis.
	 * 
	 * Defines the z-positions of voxel centers, which is important
	 * for mapping the spatial distribution of tumor cells.
	 */
	std::vector<double> coordenadas_z;
	
	/**
	 * @brief Moore neighborhood connectivity for each voxel.
	 * 
	 * Stores indices of voxels in the Moore neighborhood (26 neighbors in 3D).
	 * In cancer modeling, Moore neighborhoods allow for more complex diffusion patterns
	 * and cell-cell interactions that better represent tumor microenvironments.
	 */
	std::vector< std::vector<int> > indices_de_voxeles_conectados_tipo_moore;
	
	/**
	 * @brief Creates Moore neighborhood connections between voxels.
	 * 
	 * Sets up 26-neighbor connectivity pattern for each voxel.
	 * In cancer modeling, this enables more realistic diffusion of biochemical
	 * factors that influence cancer cell behavior.
	 */
	void crear_vecindario_moore(void);
	
	/**
	 * @brief Creates optimized Moore neighborhood connections.
	 * 
	 * Implements a more efficient version of Moore neighborhood creation.
	 * Optimization is important for large-scale cancer simulations to
	 * reduce computational overhead.
	 */
	void crear_vecindario_moore_optimizado(void);    
	
	/**
	 * @brief Creates periodic Moore neighborhood connections.
	 * 
	 * Implements periodic boundary conditions for Moore neighborhoods.
	 * In cancer modeling, periodic boundaries can simulate larger tissue
	 * contexts without computational expense of modeling the entire domain.
	 */
	void crear_vecindario_moore_periodico(void);
	
	/**
	 * @brief Creates optimized periodic Moore neighborhood connections.
	 * 
	 * Combines optimization with periodic boundary conditions.
	 * Important for efficient simulation of large-scale cancer tissues.
	 */
    void crear_vecindario_moore_periodico_optimizado();
	
	/**
	 * @brief Calculates the linear index from 3D Cartesian indices.
	 * 
	 * @param i Index in x-direction.
	 * @param j Index in y-direction.
	 * @param k Index in z-direction.
	 * @return The linear voxel index.
	 * 
	 * Fundamental for accessing voxel data in the cancer microenvironment,
	 * especially for substrate concentration and gradient calculations.
	 */
	unsigned int indice_de_voxel( unsigned int i, unsigned int j, unsigned int k );
	
	/**
	 * @brief Converts a linear index to 3D Cartesian indices.
	 * 
	 * @param n The linear voxel index.
	 * @return Vector of 3D indices [i,j,k].
	 * 
	 * Enables mapping between linear storage and 3D spatial representation
	 * of the tumor microenvironment.
	 */
	std::vector<unsigned int> indices_cartesianos( unsigned int n );

	/**
	 * @brief Spacing between grid points in x-direction.
	 * 
	 * Defines spatial resolution in x-direction, which affects
	 * the accuracy of cancer microenvironment gradient representation.
	 */
	double dx;
	
	/**
	 * @brief Spacing between grid points in y-direction.
	 * 
	 * Defines spatial resolution in y-direction, which affects
	 * the accuracy of cancer microenvironment gradient representation.
	 */
	double dy;
	
	/**
	 * @brief Spacing between grid points in z-direction.
	 * 
	 * Defines spatial resolution in z-direction, which affects
	 * the accuracy of cancer microenvironment gradient representation.
	 */
	double dz;

	/**
	 * @brief Voxel volume (dx*dy*dz).
	 * 
	 * The volume of each voxel, used in diffusion calculations.
	 * Critical for accurate modeling of substrate dynamics in the
	 * tumor microenvironment.
	 */
	double dV;
	
	/**
	 * @brief Surface area of voxel faces.
	 * 
	 * Used in flux calculations for diffusion.
	 * Important for modeling nutrient and signaling molecule
	 * movement through the cancer tissue.
	 */
	double dS;

	/**
	 * @brief Surface area of voxel faces in xy-plane.
	 * 
	 * Used in directional flux calculations.
	 * Important for modeling anisotropic diffusion in tumor tissues.
	 */
	double dS_xy;
	
	/**
	 * @brief Surface area of voxel faces in yz-plane.
	 * 
	 * Used in directional flux calculations.
	 * Important for modeling anisotropic diffusion in tumor tissues.
	 */
	double dS_yz;
	
	/**
	 * @brief Surface area of voxel faces in xz-plane.
	 * 
	 * Used in directional flux calculations.
	 * Important for modeling anisotropic diffusion in tumor tissues.
	 */
	double dS_xz;

	/**
	 * @brief Default constructor.
	 * 
	 * Creates a minimal Cartesian grid with one voxel.
	 * Serves as the starting point for more complex cancer tissue
	 * geometries.
	 */
	Grillado_Cartesiano();

	/**
	 * @brief Constructor with grid dimensions.
	 * 
	 * @param x_nodes Number of grid points in x-direction.
	 * @param y_nodes Number of grid points in y-direction.
	 * @param z_nodes Number of grid points in z-direction.
	 * 
	 * Creates a Cartesian grid with specified number of nodes in each dimension.
	 * Used to initialize cancer microenvironments with appropriate spatial resolution.
	 */
	Grillado_Cartesiano( int x_nodes, int y_nodes, int z_nodes );


	/**
	 * @brief Resizes the grid with specified number of nodes.
	 * 
	 * @param x_nodes Number of grid points in x-direction.
	 * @param y_nodes Number of grid points in y-direction.
	 * @param z_nodes Number of grid points in z-direction.
	 * 
	 * Reinitializes the grid with new dimensions while preserving the domain size.
	 * Allows adjustment of spatial resolution for cancer simulations.
	 */
	void redimensionar( int x_nodes, int y_nodes, int z_nodes );
	
	/**
	 * @brief Resizes the grid with specified domain and node counts.
	 * 
	 * @param x_start Minimum x-coordinate.
	 * @param x_end Maximum x-coordinate.
	 * @param y_start Minimum y-coordinate.
	 * @param y_end Maximum y-coordinate.
	 * @param z_start Minimum z-coordinate.
	 * @param z_end Maximum z-coordinate.
	 * @param x_nodes Number of grid points in x-direction.
	 * @param y_nodes Number of grid points in y-direction.
	 * @param z_nodes Number of grid points in z-direction.
	 * 
	 * Redefines both the spatial extent and resolution of the grid.
	 * Enables customization of the cancer microenvironment spatial domain.
	 */
	void redimensionar( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , int x_nodes, int y_nodes, int z_nodes );
	
	/**
	 * @brief Resizes the grid with specified domain and spacing.
	 * 
	 * @param x_start Minimum x-coordinate.
	 * @param x_end Maximum x-coordinate.
	 * @param y_start Minimum y-coordinate.
	 * @param y_end Maximum y-coordinate.
	 * @param z_start Minimum z-coordinate.
	 * @param z_end Maximum z-coordinate.
	 * @param dx Grid spacing in x-direction.
	 * @param dy Grid spacing in y-direction.
	 * @param dz Grid spacing in z-direction.
	 * 
	 * Redefines the spatial extent and explicitly sets grid spacing.
	 * Useful for matching cancer microenvironment scale to experimental data.
	 */
	void redimensionar( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx, double dy, double dz );

	/**
	 * @brief Finds the nearest voxel to a given position.
	 * 
	 * @param posicion The 3D position to query.
	 * @return The index of the nearest voxel.
	 * 
	 * Allows cells to interact with local microenvironment conditions.
	 * Essential for modeling how cancer cells respond to spatial heterogeneity
	 * in substrate concentrations.
	 */
	int indice_del_voxel_mas_cercano( Vector& posicion );
	
	/**
	 * @brief Gets the Cartesian indices of the voxel nearest to a position.
	 * 
	 * @param posicion The 3D position to query.
	 * @return Vector containing the i,j,k indices.
	 * 
	 * Enables spatial mapping between continuous cell positions and
	 * discretized substrate fields in cancer simulations.
	 */
	Vector indices_cartesianos_mas_cercanos(Vector& posicion );
	
	/**
	 * @brief Gets a reference to the voxel nearest to a position.
	 * 
	 * @param posicion The 3D position to query.
	 * @return Reference to the nearest voxel.
	 * 
	 * Provides direct access to voxel properties at a given location.
	 * Used by cancer cells to sense and respond to their microenvironment.
	 */
	Voxel& voxel_mas_cercano( Vector& posicion );
	
	/**
	 * @brief Gets the center coordinates of a voxel.
	 * 
	 * @param indice_de_voxel The index of the voxel.
	 * @return Vector containing the center coordinates.
	 * 
	 * Used for determining precise spatial relationships between cancer cells
	 * and voxel-based substrate distributions.
	 */
	Vector get_centro_voxel(int indice_de_voxel);

	/**
	 * @brief Displays information about the Cartesian grid.
	 * 
	 * @param os Output stream to write information to.
	 * 
	 * Provides detailed information about the grid configuration.
	 * Useful for validating cancer microenvironment spatial setup.
	 */
	void mostrar_informacion_cartesiano( std::ostream& os );
	
	/**
	 * @brief Displays the Moore neighborhood connectivity.
	 * 
	 * @param os Output stream to write information to.
	 * 
	 * Shows voxel connectivity patterns, which influence
	 * substrate diffusion patterns in the cancer microenvironment.
	 */
	void mostrar_vecindarios_moore(std::ostream& os);

};


#endif
