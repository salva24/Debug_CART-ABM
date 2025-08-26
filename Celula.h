/**
 * @file Celula.h
 *
 * @author: Luciana Melina Luque
 *
 * @details Defines the cell agents that are the core entities in the simulation.
 * This file implements both the base Celula (Cell) class and the derived Linfocito (Lymphocyte)
 * class, representing the primary agents in the agent-based model.
 * 
 * Inbound Dependencies:
 * - fstream - For I/O operations
 * - Fenotipo.h - Cell phenotype properties
 * - Vector.h - Spatial position and movement
 * - Random.h - Stochastic behaviors
 * - Microambiente.h - Environmental interactions
 * - Parametros_globales.h - Global simulation parameters
 * - Parametros.h - Cell-specific parameters
 * - Motilidad.h - Cell motility behaviors
 * 
 * Outbound Dependencies:
 * - Tejido.h - Uses cells to build tissue structure
 * - Contenedor_de_Celulas.h - Manages collections of cells
 * 
 * Usage:
 * Cells are created, managed and updated throughout the simulation:
 * Celula* celula = crear_celula();
 * celula->actualizar_posicion(dt);
 * 
 * License: Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#ifndef __CELULA_H__
#define __CELULA_H__

#include <fstream>
#include "Fenotipo.h"
#include "Vector.h"
#include "Random.h"
#include "Microambiente.h"
#include "Parametros_globales.h"
#include "Parametros.h"
#include "Motilidad.h"


extern RNG *rng;
extern Vector *v;
extern Constantes *c;
extern Parametros_globales *pg;
extern Ciclo_Modelo Ki67, vida, necrosis, apoptosis;
extern Muerte_parametros necrosis_parametros, apoptosis_parametros;


/**
 * @class Celula
 * @brief Base class for cell agents in the simulation.
 * 
 * Represents individual cells with physical properties, behavioral characteristics,
 * and interaction capabilities. Each cell has position, velocity, size, and various
 * phenotypic states. Cells can interact with the microenvironment, move according
 * to physical forces, divide, die, and secrete/consume substances.
 * 
 * Cancer Research Context:
 * This class is the foundation for modeling cancer cells, capturing their behavior,
 * proliferation patterns, and response to environmental factors like oxygen.
 * This enables researchers to study tumor growth dynamics and heterogeneity.
 */
class Celula{
	private:
		Microambiente *microambiente;

	protected:
		std::vector<double> temp_celula_fuente_sumidero_solver1;
		std::vector<double> temp_celula_fuente_sumidero_solver2;
		std::vector<double> temp_celula_fuente_sumidero_exportacion_solver1;
		std::vector<double> temp_celula_fuente_sumidero_exportacion_solver2;
		Vector velocidad_anterior;
		bool es_activa;

	public:
		int voxel_del_microambiente;
		int voxel;
		std::string nombre;
	    int tipo;
		int id;
	    int indice;
	    int madre;
	    double tiempo_desde_el_ultimo_ciclo;
	    double tiempo_desde_la_ultima_mecanica;
	    double hora_del_ultimo_ciclo;
	    double hora_de_la_ultima_mecanica;
	    double hora_global;
	    double tasa;
	    bool adherida;
	    Fenotipo fenotipo;
	    Parametros parametros;
	    Vector posicion;
    	Vector velocidad;
    	Vector desplazamiento;
        Celula* celula_adherida;        

		/**
		 * @brief Default constructor for Cell
		 * 
		 * Initializes a cell with default properties
		 */
		Celula();
		
		/**
		 * @brief Updates the cell volume based on growth/shrinkage rates
		 * 
		 * @param fenotipo The cell's phenotype containing volume properties
		 * @param dt Time step for the update
		 */
		void actualizar_volumen( Fenotipo& fenotipo, double dt);
		
		/**
		 * @brief Updates cell and death parameters based on oxygen concentration
		 * 
		 * @param fenotipo The cell's phenotype
		 * @param dt Time step for the update
		 */
		void actualizar_parametros_de_celula_y_muerte_con_o2(Fenotipo& fenotipo, double dt);
		
		/**
		 * @brief Updates cell and death parameters based on oxygen and oncoprotein
		 * 
		 * @param fenotipo The cell's phenotype
		 * @param dt Time step for the update
		 */
		void actualizar_parametros_de_celula_y_muerte_con_o2_y_oncoproteina(Fenotipo& fenotipo, double dt);
		
		/**
		 * @brief Advances the cell's phenotypic functions over time
		 * 
		 * @param hora_global Current simulation time
		 * @param dt_ciclo Time step for the cycle update
		 */
		void avanzar_funciones_del_fenotipo(double hora_global, double dt_ciclo);
		
		/**
		 * @brief Advances the cell's phenotypic functions considering oxygen
		 * 
		 * @param hora_global Current simulation time
		 * @param dt_ciclo Time step for the cycle update
		 */
		void avanzar_funciones_del_fenotipo_con_O2(double hora_global, double dt_ciclo);
		
		/**
		 * @brief Advances phenotypic functions considering oxygen and oncoprotein
		 * 
		 * @param hora_global Current simulation time
		 * @param dt_ciclo Time step for the cycle update
		 */
		void avanzar_funciones_del_fenotipo_con_O2_y_oncoproteina(double hora_global, double dt_ciclo);
		
		/**
		 * @brief Creates a daughter cell through division
		 * 
		 * @return Pointer to the newly created daughter cell
		 */
		Celula* dividir();
		
		/**
		 * @brief Initiates cell death process
		 * 
		 * @param indice Index of the death type
		 */
		void morir(int indice);
		
		/**
		 * @brief Sets the cell's position using coordinates
		 * 
		 * @param x X coordinate
		 * @param y Y coordinate
		 * @param z Z coordinate
		 * @return Success status
		 */
		bool set_posicion(double x, double y ,double z);
		
		/**
		 * @brief Sets the cell's position using a Vector
		 * 
		 * @param posicion Position vector
		 * @return Success status
		 */
		bool set_posicion(Vector posicion);
		
		/**
		 * @brief Adds mechanical potentials from another cell
		 * 
		 * @param otra_celula Another cell to interact with
		 */
		void agregar_potenciales(Celula* otra_celula);
		
		/**
		 * @brief Adds mechanical potentials from microenvironment boundary
		 */
		void agregar_potenciales_mb();
		
		/**
		 * @brief Alternative implementation for boundary potentials
		 */
		void agregar_potenciales_mb_2();
		
		/**
		 * @brief Updates cell position based on velocity
		 * 
		 * @param dt Time step for the update
		 */
		void actualizar_posicion( double dt );
		
		/**
		 * @brief Starts the cell death process
		 * 
		 * @param indice_del_ciclo_de_muerte Index of the death cycle
		 */
		void comenzar_muerte( int indice_del_ciclo_de_muerte );

		/**
		 * @brief Registers the microenvironment with the cell
		 * 
		 * @param m Pointer to the microenvironment
		 */
		void registrar_microambiente( Microambiente* );
		
		/**
		 * @brief Updates the microenvironment voxel containing the cell
		 */
		void actualizar_voxel_del_microambiente();
		
		/**
		 * @brief Sets constants for internal consumption
		 * 
		 * @param dt Time step
		 */
		void set_constantes_de_consumo_interno( double dt );
		
		/**
		 * @brief Gets the microenvironment associated with this cell
		 * 
		 * @return Pointer to the microenvironment
		 */
		Microambiente* get_microambiente( void );
		
		/**
		 * @brief Simulates secretion and consumption of microenvironment substrates
		 * 
		 * @param dt Time step
		 */
		void simular_secrecion_y_consumo( double dt );
		
		/**
		 * @brief Gets the index of the microenvironment voxel
		 * 
		 * @return Voxel index
		 */
		int get_indice_del_voxel_del_microambiente();
		
		/**
		 * @brief Gets the closest density vector from the microenvironment
		 * 
		 * @return Reference to density vector
		 */
		std::vector<double>& vector_de_densidades_mas_cercano( void );
		
		/**
		 * @brief Gets the gradient of a specific substrate
		 * 
		 * @param substrate_index Index of the substrate
		 * @return Reference to gradient vector
		 */
		std::vector<double>& gradiente_mas_cercano( int substrate_index );


		/**
		 * @brief Initializes a cell with default properties
		 */
		void inicializar_celula();
		
		/**
		 * @brief Initializes a cell with healthy properties
		 */
		void inicializar_celula_sana();
		
		/**
		 * @brief Outputs detailed cell information
		 * 
		 * @param os Output stream
		 */
		void mostrar_informacion_de_la_celula(std::ostream& os);

		/**
		 * @brief Sets whether the cell is mobile (overridden in derived classes)
		 * 
		 * @param valor Mobility flag
		 */
        virtual void es_movil(bool valor) {return;};
        
        /**
         * @brief Updates the cell's motility vector (overridden in derived classes)
         * 
         * @param dt Time step
         * @param celulas_en_mi_voxel Cells in the same voxel
         */
		virtual void actualizar_vector_de_motilidad( double dt, std::vector<Celula*> celulas_en_mi_voxel ) {return;};
		
		/**
		 * @brief Advances lymphocyte-specific behaviors (overridden in derived classes)
		 * 
		 * @param dt Time step
		 * @param celulas_en_mi_voxel Cells in the same voxel
		 */
		virtual void avanzar_linfocito( double dt, std::vector<Celula*> celulas_en_mi_voxel ) {return;};
		
		/**
		 * @brief Virtual destructor for proper inheritance
		 */
        virtual ~Celula(){};

};

/**
 * @class Linfocito
 * @brief Specialized cell class representing a lymphocyte (immune cell)
 * 
 * Extends the base Celula class with immune-specific properties and behaviors,
 * such as directed motility, cell adhesion, and the ability to trigger apoptosis
 * in other cells (immune killing).
 * 
 * Cancer Research Context:
 * This class enables modeling immune responses against cancer cells, including
 * natural killer (NK) cell activity and T-cell responses. Lymphocytes can detect,
 * adhere to, and potentially eliminate cancer cells, allowing researchers to study
 * immune surveillance and cancer immunotherapy.
 */
class Linfocito: public Celula{

    public:
        Motilidad motilidad;
        double tasa_de_asesinato;
        double tiempo_de_adhesion;
        double tasa_de_adhesion;
        double constante_elastica;
        double distancia_de_adhesion_maxima;
        double distancia_de_adhesion_minima;
        double saturacion_de_oncoproteina;
        double limite_de_oncoproteina;
        double diferencia_de_oncoproteina;
        double diferencia_de_adhesion;

	//Debug
	public:
		static long double acumulator_probabilities;
		static long double GetAccumulatedProbabilities() { return acumulator_probabilities; }

		/**
		 * @brief Constructor for Lymphocyte
		 */
        Linfocito();
        
        /**
         * @brief Sets the lymphocyte's mobility
         * 
         * @param valor Mobility flag
         */
        void es_movil(bool valor);
        
        /**
         * @brief Updates the lymphocyte's motility vector
         * 
         * @param dt Time step
         * @param celulas_en_mi_voxel Cells in the same voxel
         */
        void actualizar_vector_de_motilidad( double dt, std::vector<Celula*> celulas_en_mi_voxel );
        
        /**
         * @brief Processes lymphocyte motility behavior
         * 
         * @param dt Time step
         * @param celulas_en_mi_voxel Cells in the same voxel
         */
        void motilidad_de_linfocito( double dt, std::vector<Celula*> celulas_en_mi_voxel  );
        
        /**
         * @brief Advances lymphocyte-specific behaviors
         * 
         * @param dt Time step
         * @param celulas_en_mi_voxel Cells in the same voxel
         */
        void avanzar_linfocito( double dt, std::vector<Celula*> celulas_en_mi_voxel );
        
        /**
         * @brief Attempts to trigger apoptosis in an attached cell
         * 
         * @param celula_adherida Attached cell
         * @param dt Time step
         * @return Success status
         */
        bool intento_de_apoptosis( Celula* celula_adherida, double dt );
        
        /**
         * @brief Triggers apoptosis in an attached cell
         * 
         * @param celula_adherida Attached cell
         * @return Success status
         */
        bool desencadenar_apoptosis( Celula* celula_adherida );
        
        /**
         * @brief Adheres to a target cell
         * 
         * @param celula_objetivo Target cell
         */
        void adherir_celula(Celula* celula_objetivo);
        
        /**
         * @brief Releases an attached cell
         * 
         * @param celula_objetivo Cell to release
         */
        void soltar_celula(Celula* celula_objetivo);
        
        /**
         * @brief Checks neighboring cells for potential adhesion
         * 
         * @param celulas_en_mi_voxel Cells in the same voxel
         * @param dt Time step
         * @return Success status
         */
        bool chequear_vecinos_para_adherirse( std::vector<Celula*> celulas_en_mi_voxel , double dt );
        
        /**
         * @brief Attempts to adhere to a target cell
         * 
         * @param celula_objetivo Target cell
         * @param dt Time step
         * @return Success status
         */
        bool intentar_adherirse( Celula* celula_objetivo , double dt );
        
        /**
         * @brief Virtual destructor
         */
	virtual ~Linfocito(){};

};

extern std::vector<Celula*> todas_las_celulas;
extern std::vector<Celula*> celulas_listas_para_dividirse;
extern std::vector<Celula*> celulas_listas_para_remover;
extern std::vector<Celula*> celulas_para_registrar_en_voxeles;

/**
 * @brief Factory function to create a new cell
 * 
 * @return Pointer to newly created cell
 */
Celula* crear_celula();

/**
 * @brief Factory function to create a new lymphocyte
 * 
 * @return Pointer to newly created lymphocyte
 */
Linfocito* crear_linfocito();

#endif
