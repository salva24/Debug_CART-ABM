/**
 * @file: Microambiente.h
 *
 * @author: Luciana Melina Luque
 *
 * @description: Defines the microenvironment where cells exist and interact.
 * This file implements the diffusive microenvironment that models the biochemical
 * substrates (like oxygen, nutrients, signaling molecules) using a finite volume method.
 * 
 * Inbound Dependencies:
 * - Grillado.h - Cartesian grid for spatial discretization
 * - Parametros_globales.h - Global simulation parameters
 * 
 * Outbound Dependencies:
 * - Celula.h - Cells interact with the microenvironment
 * - Tejido.h - Tissue organization uses the microenvironment
 * 
 * Usage:
 * The microenvironment handles diffusion and tracks substrate concentrations:
 * Microambiente m;
 * m.inicializar_microambiente();
 * m.simular_difusion_decaimiento(dt);
 * 
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#ifndef __MICROAMBIENTE_H__
#define __MICROAMBIENTE_H__

#include "Grillado.h"
#include "Parametros_globales.h"

extern Vector *v;
extern Parametros_globales *pg;

typedef std::vector<double> gradiente;

/**
 * @class Microambiente
 * @brief Models the biochemical environment where cells exist
 * 
 * Implements a multi-substrate diffusive microenvironment using a finite volume method.
 * Tracks concentrations of various biochemical factors (like oxygen, nutrients, and
 * signaling molecules) and models their diffusion, decay, and interaction with cells.
 * The microenvironment is discretized into voxels on a Cartesian grid.
 * 
 * Cancer Research Context:
 * The microenvironment is critical in cancer research, as it models tumor 
 * microenvironmental conditions like hypoxia (low oxygen), nutrient gradients,
 * and signaling molecule distributions that significantly affect cancer cell
 * behavior, proliferation patterns, and treatment response.
 */
class Microambiente{

	private:
		/*! Para uso interno y acelerar los solvers */
		std::vector< std::vector<double> > temp_vectores_densidad1;
		std::vector< std::vector<double> > temp_vectores_densidad2;

		/*! para uso interno en los solvers de fuente/sumidero del bulk */
		std::vector< std::vector<double> > temp_bulk_fuente_sumidero_solver1;
		std::vector< std::vector<double> > temp_bulk_fuente_sumidero_solver2;
		std::vector< std::vector<double> > temp_bulk_fuente_sumidero_solver3;
		bool setup_del_solver_bulk_fuente_sumidero_hecho;

		/*! Almacena los punteros a las densidades.
		Se accede mediante las funciones operador(). */
		std::vector< std::vector<double> >* p_vectores_densidad;

		std::vector< std::vector<gradiente> > vectores_gradiente;
		std::vector<bool> vector_gradiente_calculado;

		/*! Útiles para los solvers.
		Redimensionar cuando se agreguen/saquen sustratos*/
		std::vector<double> uno;
		std::vector<double> cero;
		std::vector<double> un_medio;
		std::vector<double> un_tercio;

		/*! Cosas de Thomas para uso interno en los solvers de difución.*/
		std::vector< std::vector<double> > thomas_temp1;
		std::vector< std::vector<double> > thomas_temp2;
		std::vector<double> thomas_constante1x;
		std::vector<double> thomas_constante1y;
		std::vector<double> thomas_constante1z;
		std::vector<double> thomas_neg_constante1x;
		std::vector<double> thomas_neg_constante1y;
		std::vector<double> thomas_neg_constante1z;
		bool setup_de_thomas_hecho;
		int thomas_salto_en_i;
		int thomas_salto_en_j;
		int thomas_salto_en_k;
		std::vector<double> thomas_constante1;
		std::vector<double> thomas_constante1a;
		std::vector<double> thomas_constante2;
		std::vector<double> thomas_constante3;
		std::vector<double> thomas_constante3a;
		std::vector< std::vector<double> > thomas_denomx;
		std::vector< std::vector<double> > thomas_cx;
		std::vector< std::vector<double> > thomas_denomy;
		std::vector< std::vector<double> > thomas_cy;
		std::vector< std::vector<double> > thomas_denomz;
		std::vector< std::vector<double> > thomas_cz;
		bool setup_del_solver_de_difusion_hecho;

		std::vector< std::vector<double> > vector_valores_de_dirichlet;
		std::vector<bool> vector_activacion_dirichlet;


		/* Los vectores de activación se pueden especificar voxel por voxel*/
		std::vector< std::vector<bool> > dirichlet_vectorES_activacion;

	public:


		/*Grillado para las cantidades difusivas*/
		Grillado_Cartesiano mgrilla;

		//Contenedor_de_Celulas * contenedor_de_celulas;
		std::string unidades_espaciales;
		std::string unidades_temporales;
		std::string nombre;

		/*Cantidades difusivas*/
		std::vector< std::string > densidades_nombres;
		std::vector< std::string > densidades_unidades;

		// coefficients
		std::vector< double > coeficientes_de_difusion;
		std::vector< double > tasas_de_decaimiento;

		std::vector< std::vector<double> > densidades_objetivo_de_suministro_por_tasas_de_suministro;
		std::vector< std::vector<double> > tasas_de_suministro;
		std::vector< std::vector<double> > tasas_de_consumo;

		std::vector<int> voxeles_del_vaso_sanguineo;

		/**
		 * @brief Updates all rates in the microenvironment
		 */
		void actualizar_tasas( void );

		/**
		 * @brief Default constructor
		 */
		Microambiente();

		/**
		 * @brief Gets the number of tracked densities/substrates
		 * 
		 * @return Number of densities
		 */
		unsigned int numero_de_densidades( void );
		
		/**
		 * @brief Gets the total number of voxels in the domain
		 * 
		 * @return Number of voxels
		 */
		unsigned int numero_de_voxeles( void );

		/**
		 * @brief Resizes the spatial domain of the microenvironment
		 * 
		 * @param x_ini Initial x coordinate
		 * @param x_fin Final x coordinate
		 * @param y_ini Initial y coordinate
		 * @param y_fin Final y coordinate
		 * @param z_ini Initial z coordinate
		 * @param z_fin Final z coordinate
		 * @param dx_nuevo New x spacing
		 * @param dy_nuevo New y spacing
		 * @param dz_nuevo New z spacing
		 */
		void redimensionar_espacio( double x_ini, double x_fin, double y_ini, double y_fin, double z_ini, double z_fin , double dx_nuevo , double dy_nuevo , double dz_nuevo );

		/**
		 * @brief Resizes the number of tracked densities/substrates
		 * 
		 * @param nuevo_tamano New number of densities
		 */
		void redimensionar_densidades( int nuevo_tamano );
		
		/**
		 * @brief Adds a new density/substrate with default properties
		 */
		void agregar_densidad( void );
		
		/**
		 * @brief Adds a new density/substrate with specified properties
		 * 
		 * @param nombre Name of the density
		 * @param unidades Units of the density
		 * @param coeficiente_de_difusion Diffusion coefficient
		 * @param tasa_de_decaimiento Decay rate
		 */
		void agregar_densidad (std::string nombre , std::string unidades, double coeficiente_de_difusion, double tasa_de_decaimiento );

		/**
		 * @brief Sets properties for an existing density/substrate
		 * 
		 * @param indice Index of the density
		 * @param nombre Name of the density
		 * @param unidades Units of the density
		 * @param coeficiente_de_difusion Diffusion coefficient
		 * @param tasa_de_decaimiento Decay rate
		 */
		void set_densidad(int indice, std::string nombre , std::string unidades , double coeficiente_de_difusion , double tasa_de_decaimiento );
		
		/**
		 * @brief Sets name and units for an existing density/substrate
		 * 
		 * @param indice Index of the density
		 * @param nombre Name of the density
		 * @param unidades Units of the density
		 */
		void set_densidad(int indice, std::string nombre , std::string unidades );

		/**
		 * @brief Finds the index of a density by its name
		 * 
		 * @param nombre Name of the density
		 * @return Index of the density, or -1 if not found
		 */
		int encontrar_indice_de_densidad(std::string nombre);
		
		/**
		 * @brief Converts 3D Cartesian indices to a voxel index
		 * 
		 * @param i x-index
		 * @param j y-index
		 * @param k z-index
		 * @return Voxel index
		 */
		int indice_de_voxel( int i, int j, int k );
		std::vector<unsigned int> indices_cartesianos( int n );

		/**
		 * @brief Finds the voxel index closest to a given position
		 * 
		 * @param posicion 3D position
		 * @return Closest voxel index
		 */
		int indice_del_voxel_mas_cercano( Vector& posicion);
		
		/**
		 * @brief Gets the center position of a voxel
		 * 
		 * @param indice_del_voxel Voxel index
		 * @return Center position as a Vector
		 */
		Vector centro_del_voxel(int indice_del_voxel);
		
		/**
		 * @brief Gets Cartesian indices closest to a position
		 * 
		 * @param posicion 3D position
		 * @return Vector containing the closest i,j,k indices
		 */
		Vector indices_cartesianos_mas_cercanos(Vector& posicion );
		Voxel& voxel_mas_cercano( Vector& posicion );
		Voxel& voxeles( int indice_de_voxel );
		std::vector<double>& vector_de_densidades_mas_cercano( Vector& posicion );
		std::vector<double>& vector_de_densidades_mas_cercano( int indice_de_voxel );

		/*! acceder al vector de densidades en [ X(i),Y(j),Z(k) ] */
		std::vector<double>& operator()( int i, int j, int k );
		/*! acceder al vector de densidades en el voxel (n) */
		std::vector<double>& operator()( int n );
		/*! acceder al vector de densidades en  [ X(i),Y(j),Z(k) ] */
		std::vector<double>& vector_de_densidades( int i, int j, int k );
		/*! acceder al vector de densidades en  el voxel (n) */
		std::vector<double>& vector_de_densidades( int n );


		std::vector<gradiente>& vector_de_gradientes(int i, int j, int k);
		std::vector<gradiente>& vector_de_gradientes(int n );

		std::vector<gradiente>& vector_de_gradiente_mas_cercano( Vector& posicion );

		void calcular_todos_los_vectores_de_gradientes( void );
		
		/**
		 * @brief Calculates gradient for a specific density
		 * 
		 * @param n Index of the density
		 */
		void calcular_vector_de_gradiente( int n );
		
		/**
		 * @brief Resets all gradient calculations
		 */
		void resetear_todos_los_vectores_de_gradientes( void );

		// Solvers

		/**
		 * @brief Simulates diffusion and decay of substrates
		 * 
		 * @param dt Time step
		 */
		void simular_difusion_decaimiento( double dt );

		/**
		 * @brief Simulates bulk supply and consumption
		 * 
		 * @param dt Time step
		 */
		void simular_suministro_y_consumo_del_bulk( double dt );

		/**
		 * @brief Simulates cell-based supply and consumption
		 * 
		 * @param dt Time step
		 */
		void simular_suministro_y_consumo_de_las_celulas( double dt );

		/**
		 * @brief Outputs information about the microenvironment
		 * 
		 * @param os Output stream
		 */
		void mostrar_informacion( std::ostream& os );
		
		// Dirichlet boundary condition methods
		
		/**
		 * @brief Adds a Dirichlet node at a voxel
		 * 
		 * @param indice_de_voxel Voxel index
		 * @param valor Vector of density values
		 */
		void agregar_nodo_de_dirichlet( int indice_de_voxel, std::vector<double>& valor );
		
		/**
		 * @brief Updates a Dirichlet node with new values
		 * 
		 * @param indice_de_voxel Voxel index
		 * @param nuevo_valor Vector of new density values
		 */
		void actualizar_nodo_de_dirichlet( int indice_de_voxel , std::vector<double>& nuevo_valor );
		
		/**
		 * @brief Updates a specific substrate at a Dirichlet node
		 * 
		 * @param indice_de_voxel Voxel index
		 * @param indice_de_sustrato Substrate index
		 * @param nuevo_valor New density value
		 */
		void actualizar_nodo_de_dirichlet( int indice_de_voxel , int indice_de_sustrato, double nuevo_valor );
		
		/**
		 * @brief Removes a Dirichlet node
		 * 
		 * @param indice_de_voxel Voxel index
		 */
		void remover_nodo_de_dirichlet( int indice_de_voxel );
		
		/**
		 * @brief Applies Dirichlet conditions to the domain
		 */
		void aplicar_condiciones_de_dirichlet( void );

		/**
		 * @brief Sets activation status for a substrate's Dirichlet condition
		 * 
		 * @param indice_de_sustrato Substrate index
		 * @param nuevo_valor New activation status
		 */
		void set_activacion_de_sustrato_de_dirichlet( int indice_de_sustrato , bool nuevo_valor );

		/**
		 * @brief Sets activation status for a substrate at a specific voxel
		 * 
		 * @param indice_de_sustrato Substrate index
		 * @param indice Voxel index
		 * @param nuevo_valor New activation status
		 */
		void set_activacion_de_sustrato_de_dirichlet( int indice_de_sustrato , int indice, bool nuevo_valor );
		
		/**
		 * @brief Sets activation status for all substrates at a voxel
		 * 
		 * @param indice Voxel index
		 * @param nuevo_valor Vector of activation statuses
		 */
		void set_activacion_de_sustrato_de_dirichlet( int indice, std::vector<bool>& nuevo_valor );
		
		/**
		 * @brief Gets activation status for a substrate at a voxel
		 * 
		 * @param indice_de_sustrato Substrate index
		 * @param indice Voxel index
		 * @return Activation status
		 */
		bool get_activacion_de_sustrato_de_dirichlet( int indice_de_sustrato, int indice );

		/**
		 * @brief Checks if a voxel is a Dirichlet node
		 * 
		 * @param indice_de_voxel Voxel index
		 * @return Reference to Dirichlet status
		 */
		bool& es_nodo_de_dirichlet( int indice_de_voxel );

		/**
		 * @brief Solver for diffusion-decay using LOD (Locally One-Dimensional) method
		 * 
		 * @param dt Time step
		 */
		void solver_decaimiento_de_la_difusion__coeficientes_constantes_LOD_3D( double dt );

		/**
		 * @brief Initializes the microenvironment with default settings
		 */
		void inicializar_microambiente();

		/**
		 * @brief Creates a blood vessel in the specified region
		 * 
		 * @param xmin Minimum x-coordinate
		 * @param ymin Minimum y-coordinate
		 * @param zmin Minimum z-coordinate
		 * @param xmax Maximum x-coordinate
		 * @param ymax Maximum y-coordinate
		 * @param zmax Maximum z-coordinate
		 */
		void crear_vaso_sanguineo(int xmin, int ymin, int zmin, int xmax, int ymax, int zmax);
		
		/** @cond */
		// Operator for accessing densities
		/** @endcond */
};

/**
 * @brief Sets default properties for a microenvironment
 * 
 * @param M Pointer to the microenvironment to configure
 */
void set_microambiente_default( Microambiente* M );
Microambiente* get_microambiente_default( void );


#endif
