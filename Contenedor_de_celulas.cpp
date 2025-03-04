/**
 * @file: Contenedor_de_celulas.cpp
 *
 * @author: Luciana Melina Luque
 *
 * @details Implementation of the cell container system that manages collections of cells
 * in the cancer simulation. This file implements the mechanics of cell organization, 
 * interaction, and lifecycle events.
 * 
 * Inbound Dependencies:
 * - Contenedor_de_Celulas.h - Header defining the cell container interface
 * 
 * Cancer Research Context:
 * The cell container implementation is critical for scaling cancer simulations to biologically
 * relevant cell populations. It enables efficient spatial organization, neighbor detection, and
 * cell-cell interaction which are all essential aspects of modeling tumor growth dynamics.
 * This implementation particularly supports:
 * - Modeling mechanical interactions between cancer cells
 * - Cell division and tumor growth dynamics
 * - Cell death processes
 * - Cell population dynamics in heterogeneous tumors
 * - Spatial organization of cancer and immune cells
 * 
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#include "Contenedor_de_Celulas.h"

//std::ofstream outfile2 ("Division.dat", std::ios::app );

/**
 * @brief Constructor for the cell container
 * 
 * Initializes a new cell container with default values for tracking mechanical
 * updates and division events. In cancer modeling, proper initialization ensures
 * accurate tracking of tumor growth rates and mechanical interactions.
 */
Contenedor_de_Celulas::Contenedor_de_Celulas(){



	num_de_divisiones_en_este_paso=0;
    tiempo_desde_la_ultima_mecanica=0.0;
    hora_de_la_ultima_mecanica=0.0;    


	return;
}

/**
 * @brief Initializes the spatial domain of the cell container
 * 
 * Sets up the Cartesian grid with the specified dimensions and
 * creates the data structures needed to track cells within voxels.
 * For cancer modeling, this defines the tissue domain within which
 * tumor growth can occur and determines the spatial resolution of
 * the simulation.
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
void Contenedor_de_Celulas::inicializar(double x_ini, double x_fin, double y_ini, double y_fin, double z_ini, double z_fin , double dx, double dy, double dz)
{
	grillado.redimensionar(x_ini, x_fin, y_ini, y_fin, z_ini, z_fin , dx, dy, dz);
	celulas_en_voxel.resize(grillado.voxeles.size());
    for(unsigned int v=0; v<grillado.voxeles.size(); v++){
        celulas_en_voxel[v].resize(0);
    }
    


	return;
}

/**
 * @brief Registers a cell in the container
 * 
 * Places the cell in the appropriate voxel based on its position.
 * If the position is invalid, the cell is marked as outside the domain.
 * This function is essential for maintaining spatial organization in
 * tumor simulations, especially when tracking cells that migrate to
 * different regions of the tissue.
 * 
 * @param celula Pointer to the cell to register
 */
void Contenedor_de_Celulas::registrar_celula( Celula* celula ){

	if( !grillado.es_valida_la_posicion(celula->posicion.x,celula->posicion.y,celula->posicion.z))
	{
    std::cout << "Posicion invalida: " << celula->posicion.x << " " << celula->posicion.y << " " << celula->posicion.z << " " << celula->id << "\n ";
		celula->voxel =-1;
		agregar_celula_a_voxel(celula, celula->voxel);
		return;
	}

	celula->voxel = grillado.indice_del_voxel_mas_cercano( celula->posicion );

	agregar_celula_a_voxel(celula, celula->voxel);

	return;

}

/**
 * @brief Adds a cell to a specific voxel
 * 
 * Updates the spatial data structure to include the cell in the specified voxel.
 * Efficient voxel-based organization is critical for cancer simulations as it
 * allows for quick identification of neighboring cells for interaction calculations.
 * 
 * @param celula Pointer to the cell to add
 * @param indice_de_voxel Index of the voxel to add the cell to
 */
void Contenedor_de_Celulas::agregar_celula_a_voxel(Celula* celula, int indice_de_voxel){

	celulas_en_voxel[indice_de_voxel].push_back(celula);

}

/**
 * @brief Calculates mechanical forces between two cells
 * 
 * Computes repulsive and adhesive forces between cells based on
 * their distance, size, and mechanical properties. Updates velocity
 * vectors for both cells accordingly.
 * 
 * In cancer research, these mechanical interactions are critical for modeling:
 * - Cell crowding effects in tumors
 * - Differential adhesion hypothesis in tumor organization
 * - Mechanical pressure effects on cancer cell behavior
 * - Tumor compaction and invasion mechanics
 * 
 * @param celula Pointer to the first cell
 * @param otra_celula Pointer to the second cell
 */
void Contenedor_de_Celulas::agregar_potenciales_cdc(Celula* celula, Celula* otra_celula){

	if( celula->id == otra_celula->id ){
		return;
	}

//	double distancia = 0;
//	for( int i = 0 ; i < 3 ; i++ )
//	{
//		desplazamiento[i] = posicion[i] - (*otra_celula).posicion[i];
//		distancia += desplazamiento[i] * desplazamiento[i];
//	}

	celula->desplazamiento.x = celula->posicion.x - otra_celula->posicion.x;
	celula->desplazamiento.y = celula->posicion.y - otra_celula->posicion.y;
	celula->desplazamiento.z = celula->posicion.z - otra_celula->posicion.z;

    //Agrego esto por las Condiciones Periódicas de Contorno (funciona también sin CPC)
    celula->desplazamiento.x = celula->desplazamiento.x - (2*pg->rango_en_X[1])*round(celula->desplazamiento.x/(2*pg->rango_en_X[1]));
	celula->desplazamiento.y = celula->desplazamiento.y - (2*pg->rango_en_Y[1])*round(celula->desplazamiento.y/(2*pg->rango_en_Y[1]));
	celula->desplazamiento.z = celula->desplazamiento.z - (2*pg->rango_en_Z[1])*round(celula->desplazamiento.z/(2*pg->rango_en_Z[1]));

    

    
	double distancia = celula->desplazamiento.x * celula->desplazamiento.x +
	celula->desplazamiento.y * celula->desplazamiento.y + celula->desplazamiento.z * celula->desplazamiento.z;



	distancia = std::max(sqrt(distancia), 0.00001);



	//Repulsive
	double R = celula->fenotipo.geometria.radio + otra_celula->fenotipo.geometria.radio;

	double temp_r;
	if( distancia > R )
	{
		temp_r=0;
	}
	else
	{
		// temp_r = 1 - distance/R;
		temp_r = -distancia; // -d
		temp_r /= R; // -d/R
		temp_r += 1.0; // 1-d/R
		temp_r *= temp_r; // (1-d/R)^2

	}

	if(celula->tipo == otra_celula->tipo){
	double repulsion_efectiva = sqrt( celula->fenotipo.mecanica.fuerza_de_repulsion_cc * otra_celula->fenotipo.mecanica.fuerza_de_repulsion_cc );
	temp_r *= repulsion_efectiva;
    }else{
	double repulsion_efectiva = sqrt( celula->fenotipo.mecanica.fuerza_de_repulsion_co * otra_celula->fenotipo.mecanica.fuerza_de_repulsion_co );
	temp_r *= repulsion_efectiva;
    }



	// Adhesive



	double distancia_de_interaccion_max = celula->fenotipo.mecanica.distancia_de_adhesion_maxima_relativa * celula->fenotipo.geometria.radio +
		otra_celula->fenotipo.mecanica.distancia_de_adhesion_maxima_relativa * otra_celula->fenotipo.geometria.radio;

	if(distancia < distancia_de_interaccion_max )
	{

		double temp_a = -distancia; // -d
		temp_a /= distancia_de_interaccion_max; // -d/S
		temp_a += 1.0; // 1 - d/S
		temp_a *= temp_a; // (1-d/S)^2



		if(celula->tipo == otra_celula->tipo){
		double adhesion_efectiva = sqrt( celula->fenotipo.mecanica.fuerza_de_adhesion_cc * otra_celula->fenotipo.mecanica.fuerza_de_adhesion_cc );
		temp_a *= adhesion_efectiva;
		}else{
		double adhesion_efectiva = sqrt( celula->fenotipo.mecanica.fuerza_de_adhesion_co * otra_celula->fenotipo.mecanica.fuerza_de_adhesion_co );
		temp_a *= adhesion_efectiva;
		}

		temp_r -= temp_a;
	}

	if( fabs(temp_r) < 1e-16 )
	{ return; }
	temp_r /= distancia;




	
    //axpy(&velocidad, temp_r, desplazamiento);
    celula->velocidad.x += temp_r * celula->desplazamiento.x;
	celula->velocidad.y += temp_r * celula->desplazamiento.y;
	celula->velocidad.z += temp_r * celula->desplazamiento.z;
    
    //Le agrego el potencial a la otra célula también
    otra_celula->velocidad.x -= temp_r * celula->desplazamiento.x;
	otra_celula->velocidad.y -= temp_r * celula->desplazamiento.y;
	otra_celula->velocidad.z -= temp_r * celula->desplazamiento.z; 
    

    return;
}

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
 * In cancer research, this method drives the entire simulation forward,
 * capturing key tumor behaviors such as:
 * - Tumor growth through cell division
 * - Cell death through various mechanisms
 * - Mechanical interactions within the tumor
 * - Phenotypic evolution of cancer cells
 * - Interaction with the microenvironment
 * 
 * @param tiempo_total Current simulation time
 * @param dt_difusion Time step for diffusion
 * @param dt_mecanico Time step for mechanics
 * @param dt_ciclo Time step for cell cycle
 */
void Contenedor_de_Celulas::actualizar_todas_las_celulas(double tiempo_total, double dt_difusion, double dt_mecanico, double dt_ciclo){



    //#pragma omp parallel for
	for (unsigned int l=0; l< todas_las_celulas.size(); l++ ){

		todas_las_celulas[l]->simular_secrecion_y_consumo(dt_difusion);

	}




    static double tolerancia_de_la_mecanica = 0.001 * dt_mecanico;
    
    tiempo_desde_la_ultima_mecanica = tiempo_total - hora_de_la_ultima_mecanica;
    if(fabs(tiempo_desde_la_ultima_mecanica - dt_mecanico) < tolerancia_de_la_mecanica){
        
        for(unsigned int i=0; i<grillado.voxeles.size(); i++){
            for(unsigned int j=0; j<celulas_en_voxel[i].size(); j++){
                for(unsigned int k=j+1; k<celulas_en_voxel[i].size(); k++){
                    agregar_potenciales_cdc(celulas_en_voxel[i][j], celulas_en_voxel[i][k]);
                }
                for(unsigned int n=0; n<grillado.indices_de_voxeles_conectados_tipo_moore[i].size(); n++){
                    for(unsigned int k=0; k<celulas_en_voxel[grillado.indices_de_voxeles_conectados_tipo_moore[i][n]].size(); k++){
                        agregar_potenciales_cdc(celulas_en_voxel[i][j], celulas_en_voxel[grillado.indices_de_voxeles_conectados_tipo_moore[i][n]][k]);
                    }
                }
                
            }

        }
        hora_de_la_ultima_mecanica = tiempo_total;
        

//        #pragma omp parallel for 
        for( unsigned int a=0; a < todas_las_celulas.size(); a++ ){
            //todas_las_celulas[a]->agregar_potenciales_mb_2();
            if(todas_las_celulas[a]->tipo == 2){

                todas_las_celulas[a]->actualizar_vector_de_motilidad(dt_mecanico, celulas_en_voxel[todas_las_celulas[a]->voxel]);
                todas_las_celulas[a]->avanzar_linfocito(dt_mecanico, celulas_en_voxel[todas_las_celulas[a]->voxel]);
            }
        }


        #pragma omp parallel for 
		for( unsigned int a=0; a < todas_las_celulas.size(); a++ ){
			todas_las_celulas[a]->actualizar_posicion(dt_mecanico);
            todas_las_celulas[a]->tiempo_desde_la_ultima_mecanica = tiempo_total - todas_las_celulas[a]->hora_de_la_ultima_mecanica;
            todas_las_celulas[a]->hora_de_la_ultima_mecanica = tiempo_total;
		}		        
    }


/*
	for( unsigned int a=0; a < todas_las_celulas.size(); a++ ){
//		    std::cout << "celula: " << todas_las_celulas[a]->id << "\n";

		todas_las_celulas[a]->tiempo_desde_la_ultima_mecanica = tiempo_total - todas_las_celulas[a]->hora_de_la_ultima_mecanica;

		//double tasa;
		//std::cout << "Celula: " << todas_las_celulas[a]->id << " Tiempo de la ultima mec: " << todas_las_celulas[a]->tiempo_desde_la_ultima_mecanica << std::endl;
		//std::cout<<" = "<<fabs(todas_las_celulas[a]->tiempo_desde_la_ultima_mecanica - dt_mecanico)<<"\n";
		if(fabs(todas_las_celulas[a]->tiempo_desde_la_ultima_mecanica - dt_mecanico) < tolerancia_de_la_mecanica){


		//std::cout << "Entre en la mecanica \n";

            todas_las_celulas[a]->agregar_potenciales_mb_2();

			for ( unsigned int b=0; b < celulas_en_voxel[todas_las_celulas[a]->voxel].size(); b++){

				todas_las_celulas[a]->agregar_potenciales(celulas_en_voxel[todas_las_celulas[a]->voxel][b]);
//				std::cout << "Agregue potenciales en la celulaA " << todas_las_celulas[a]->id << std::endl;
			}

			for( unsigned int c=0; c< grillado.indices_de_voxeles_conectados_tipo_moore[todas_las_celulas[a]->voxel].size(); c++){

				int voxel_vecino = grillado.indices_de_voxeles_conectados_tipo_moore[todas_las_celulas[a]->voxel][c];

				if(contiene_alguna_celula(voxel_vecino)){

					for( unsigned int d=0; d<celulas_en_voxel[voxel_vecino].size(); d++){

						todas_las_celulas[a]->agregar_potenciales(celulas_en_voxel[voxel_vecino][d]);
//						std::cout << "Agregue potenciales en la celulaB " << todas_las_celulas[a]->id << std::endl;
					}
				}
			}

			if(todas_las_celulas[a]->tipo == 2){
				//std::cout << "Entre al if de linfocitos \n";
                todas_las_celulas[a]->actualizar_vector_de_motilidad(dt_mecanico, celulas_en_voxel[todas_las_celulas[a]->voxel]);
               //std::cout << "Entre al if de linfocitos 1 \n";
                todas_las_celulas[a]->avanzar_linfocito(dt_mecanico, celulas_en_voxel[todas_las_celulas[a]->voxel]);
                //std::cout << "Entre al if de linfocitos 2 \n";
                //std::cin.get();
			}


			todas_las_celulas[a]->hora_de_la_ultima_mecanica = tiempo_total;

		}
	}
	
*/








/////////////////////////Avanzo en las funciones del fenotipo
	static double tolerancia_del_fenotipo = 0.001 * dt_ciclo;
    //#pragma omp parallel for 
	for( unsigned int i=0; i < todas_las_celulas.size(); i++ ){

        todas_las_celulas[i]->tiempo_desde_el_ultimo_ciclo = tiempo_total - todas_las_celulas[i]->hora_del_ultimo_ciclo;


        if(fabs(todas_las_celulas[i]->tiempo_desde_el_ultimo_ciclo - dt_ciclo) < tolerancia_del_fenotipo){




            todas_las_celulas[i]->avanzar_funciones_del_fenotipo_con_O2_y_oncoproteina(tiempo_total, dt_ciclo);

            todas_las_celulas[i]->hora_del_ultimo_ciclo = tiempo_total;
		}




	}




	for( unsigned int i=0; i < celulas_listas_para_dividirse.size(); i++ ){


		celulas_listas_para_dividirse[i]->dividir();

	}

///////////////////Registro las células en los vóxeles correspondientes.
	for( unsigned int i=0; i < celulas_para_registrar_en_voxeles.size(); i++ ){

		registrar_celula(celulas_para_registrar_en_voxeles[i]);

	}

//////////////////////Remuevo las celulas correspondientes

	for( unsigned int i=0; i < celulas_listas_para_remover.size(); i++ ){

        
 
		sacar_celula_de_voxel(celulas_listas_para_remover[i], celulas_listas_para_remover[i]->voxel);
		celulas_listas_para_remover[i]->morir(celulas_listas_para_remover[i]->indice);

	}


	num_de_divisiones_en_este_paso += celulas_listas_para_dividirse.size();
	num_de_muertes_en_este_paso += celulas_listas_para_remover.size();
	celulas_listas_para_dividirse.clear();
	celulas_listas_para_remover.clear();
	celulas_para_registrar_en_voxeles.clear();





	return;

}

/**
 * @brief Checks if a voxel contains any cells
 * 
 * Determines whether the specified voxel contains at least one cell.
 * This is used to optimize mechanical calculations by only processing
 * voxels that contain cells. In cancer simulations, this helps improve
 * performance when modeling sparse regions of the tumor.
 * 
 * @param indice_de_voxel Index of the voxel to check
 * @return true if the voxel contains at least one cell, false otherwise
 */
bool Contenedor_de_Celulas::contiene_alguna_celula(int indice_de_voxel){
	return celulas_en_voxel[indice_de_voxel].size()==0?false:true;
}

/**
 * @brief Removes a cell from a voxel
 * 
 * Used when cells move between voxels or when they die.
 * The method updates the spatial data structure by removing
 * the reference to the cell from the voxel. In cancer modeling,
 * this is essential for tracking cell migration, division, and death.
 * 
 * @param celula Pointer to the cell to remove
 * @param indice_de_voxel Index of the voxel containing the cell
 */
void Contenedor_de_Celulas::sacar_celula_de_voxel(Celula* celula, int indice_de_voxel){

	int indice_a_borrar = 0;
	while( celulas_en_voxel[indice_de_voxel][ indice_a_borrar ] != celula )
	{
		indice_a_borrar++;
	}



	celulas_en_voxel[celula->voxel][indice_a_borrar] = celulas_en_voxel[celula->voxel][celulas_en_voxel[celula->voxel].size()-1 ];

	celulas_en_voxel[celula->voxel].pop_back();
	return;

}

/**
 * @brief Updates voxel assignments for all cells
 * 
 * Checks each cell's position and moves it to the correct voxel
 * if it has moved to a different voxel. This is typically called
 * periodically rather than after every position update for efficiency.
 * 
 * In cancer modeling, this method ensures proper spatial organization
 * of cells, which is critical for accurate simulation of:
 * - Cell migration through tissue
 * - Tumor invasion into surrounding tissue
 * - Spatial heterogeneity within tumors
 * - Local interactions between cancer and immune cells
 */
void Contenedor_de_Celulas::actualizar_voxeles_de_celulas(){

	int nuevo_voxel = 0;
	for( unsigned int a=0; a<todas_las_celulas.size(); a++){

		if( !grillado.es_valida_la_posicion(todas_las_celulas[a]->posicion.x,todas_las_celulas[a]->posicion.y,todas_las_celulas[a]->posicion.z)){
			nuevo_voxel=-1;
		}else{
			nuevo_voxel= grillado.indice_del_voxel_mas_cercano( todas_las_celulas[a]->posicion );

		}

		if(todas_las_celulas[a]->voxel != nuevo_voxel){

			// #pragma omp critical{
			sacar_celula_de_voxel(todas_las_celulas[a], todas_las_celulas[a]->voxel);
			agregar_celula_a_voxel(todas_las_celulas[a], nuevo_voxel);

			todas_las_celulas[a]->voxel = nuevo_voxel;

		}


	}
	return;

}
