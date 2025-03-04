/**
 * @file Tejido.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of the Tejido (Tissue) class for cancer simulation
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 * 
 * @details
 * Inbound Dependencies:
 * - Tejido.h: Header file defining the Tejido class
 * - Microambiente.h: Microenvironment implementation for biochemical factors
 * - Contenedor_de_Celulas.h: Cell container for spatial organization
 * - Celula.h: Base cell class and specialized cell types
 * - Vector.h: Vector mathematics for spatial positioning
 * 
 * Outbound Dependencies:
 * - Main.cpp: Uses Tejido for simulation control
 * 
 * Cancer Research Context:
 * This file implements the tissue-level organization for cancer simulations,
 * providing methods to create and manage cellular structures that model tumor
 * growth and development. It handles the initialization of the microenvironment,
 * placement of cancer and normal cells, and creation of specific tissue geometries
 * relevant to cancer research. The implementation supports studies of tumor
 * morphology, invasion patterns, and interactions with immune cells.
 */
#include "Tejido.h"

/**
 * @brief Constructor for the Tejido class
 * 
 * Cancer Research Context:
 * Initializes a tissue object that serves as the foundation for cancer simulations.
 * The tissue integrates cellular components with the biochemical microenvironment,
 * providing the structural framework for modeling tumor development, invasion,
 * and interactions with surrounding normal tissue.
 */
Tejido::Tejido(){

//	std::cout << "entre a tejido \n";
//	microambiente;
//	cdc;
}

/**
 * @brief Initializes the tissue with cells and microenvironment
 * 
 * Cancer Research Context:
 * This method establishes the initial conditions for cancer simulation by:
 * 1. Setting up the biochemical microenvironment (oxygen, nutrients, signaling molecules)
 * 2. Initializing the spatial cell container for efficient cell interactions
 * 3. Creating and placing the initial cancer cell
 * 4. Generating a sphere of normal cells surrounding the tumor
 * 
 * This initialization creates a controlled environment for studying early tumor
 * development, cellular interactions, and the influence of spatial organization
 * on cancer progression.
 */
void Tejido::inicializar_tejido(){

    //Microambiente
	microambiente.inicializar_microambiente();

	//Contenedor de células
	cdc.inicializar(pg->rango_en_X[0],pg->rango_en_X[1],pg->rango_en_Y[0],pg->rango_en_Y[1],
	pg->rango_en_Z[0],pg->rango_en_Z[1], pg->c_dx, pg->c_dy, pg->c_dz);

    //Célula Cancerosa
	cdc.celula = crear_celula();
	cdc.celula->inicializar_celula();
    cdc.celula->fenotipo.secrecion.oncoproteina = rng->NormalRandom_CM( pg->imm_mean, pg->imm_sd );
	cdc.registrar_celula(cdc.celula);
	celulas_para_registrar_en_voxeles.clear();

	//Planos de células sanas
	std::vector<Vector> posiciones = crear_esfera_de_celulas(150);
	Celula* pCelula = NULL;


    for( long unsigned int i=0; i < posiciones.size(); i++ ){
        if(i==0){
            cdc.celula->set_posicion( posiciones[i] );
            cdc.registrar_celula(cdc.celula);
            celulas_para_registrar_en_voxeles.clear();            
        }else{

            pCelula = crear_celula();
            pCelula->inicializar_celula();
            pCelula->set_posicion( posiciones[i] );
            cdc.registrar_celula(pCelula);
            celulas_para_registrar_en_voxeles.clear();
            }

    }


	return;


}


/**
 * @brief Creates planes of healthy cells along the Z-axis
 * 
 * @param cantidad Number of planes to create
 * @param posicion_en_z_del_primer_plano Z-position of the first plane
 * @return std::vector<Vector> Vector of 3D positions for cells
 * 
 * Cancer Research Context:
 * This method generates organized layers of healthy cells that can represent
 * normal tissue surrounding a tumor. This arrangement is valuable for studying:
 * - Tumor-stroma interactions at the invasion front
 * - Mechanical constraints imposed by surrounding tissue
 * - Diffusion of nutrients and drugs through layered tissue
 * - Patterns of cancer invasion into structured normal tissue
 * 
 * The hexagonal arrangement within each plane mimics the efficient packing
 * often seen in epithelial tissues, providing a biologically relevant
 * environment for cancer growth and invasion studies.
 */
std::vector<Vector> Tejido::crear_planos_de_celulas_sanas_en_Z(int cantidad, double posicion_en_z_del_primer_plano){

	std::vector<Vector> posiciones;
    int xc=0,yc=0,zc=0;
	double espacio_en_x = cdc.celula->fenotipo.geometria.radio*sqrt(3);
	double espacio_en_y = cdc.celula->fenotipo.geometria.radio*2;
	double espacio_en_z = cdc.celula->fenotipo.geometria.radio*sqrt(3);

    Vector punto_temporal;

	for(double z=0; z<=cantidad; z++, zc++)
	{
		for(double x=pg->rango_en_X[0]; x<pg->rango_en_X[1]; x+=espacio_en_x, xc++)
		{
			for(double y=pg->rango_en_Y[0]; y<pg->rango_en_Y[1]; y+=espacio_en_y, yc++)
			{
				punto_temporal.x= x + (zc%2) * 0.5 * cdc.celula->fenotipo.geometria.radio;
				punto_temporal.y= y + (xc%2) * cdc.celula->fenotipo.geometria.radio;
				punto_temporal.z= posicion_en_z_del_primer_plano + z*espacio_en_z;

				posiciones.push_back(punto_temporal);
			}

		}
	}
	return posiciones;



}

/**
 * @brief Creates a sphere of cell positions with specified radius
 * 
 * @param radio_de_la_esfera Radius of the sphere to create
 * @return std::vector<Vector> Vector of 3D positions for cells
 * 
 * Cancer Research Context:
 * This method generates a spherical arrangement of cell positions that can be used
 * to create realistic 3D tumor geometries. Spherical tumor models are important in
 * cancer research for studying:
 * - Diffusion gradients of oxygen and nutrients from the periphery to the core
 * - Development of hypoxic and necrotic regions in the tumor center
 * - Spatial heterogeneity in proliferation rates and phenotypes
 * - Mechanical constraints on tumor growth and invasion
 * 
 * The method creates a hexagonal close-packed (HCP) lattice structure that
 * efficiently fills 3D space while maintaining realistic cell-cell spacing.
 */
std::vector<Vector> Tejido::crear_esfera_de_celulas(double radio_de_la_esfera){

	std::vector<Vector> posiciones;
	int xc=0,yc=0,zc=0;
	double espacio_en_x = cdc.celula->fenotipo.geometria.radio*sqrt(3);
	double espacio_en_y = cdc.celula->fenotipo.geometria.radio*2;
	double espacio_en_z = cdc.celula->fenotipo.geometria.radio*sqrt(3);

	Vector punto_temporal;


	for(double z=-radio_de_la_esfera;z<radio_de_la_esfera;z+=espacio_en_z, zc++)
	{
		for(double x=-radio_de_la_esfera;x<radio_de_la_esfera;x+=espacio_en_x, xc++)
		{
			for(double y=-radio_de_la_esfera;y<radio_de_la_esfera;y+=espacio_en_y, yc++)
			{
				punto_temporal.x=x + (zc%2) * 0.5 * cdc.celula->fenotipo.geometria.radio;
				punto_temporal.y=y + (xc%2) * cdc.celula->fenotipo.geometria.radio;
				punto_temporal.z=z;

				if(sqrt(punto_temporal.x*punto_temporal.x+punto_temporal.y*punto_temporal.y+punto_temporal.z*punto_temporal.z)< radio_de_la_esfera)
				{ posiciones.push_back(punto_temporal); }
			}

		}
	}
	return posiciones;


}


/**
 * @brief Introduces lymphocytes (immune cells) around the tumor periphery
 * 
 * @param cantidad Number of lymphocytes to introduce
 * 
 * Cancer Research Context:
 * This method models the immune response to cancer by introducing lymphocytes
 * in a ring-like region surrounding the tumor. This approach is critical for studying:
 * - Tumor-immune interactions at the invasion front
 * - Immune surveillance and its role in controlling tumor growth
 * - Immune cell infiltration patterns into solid tumors
 * - Potential for immunotherapy efficacy in different tumor microenvironments
 * 
 * The method calculates the current tumor radius and places lymphocytes at a
 * specified distance from the tumor boundary, mimicking the immune infiltration
 * patterns observed in many solid tumors.
 */
void Tejido::introducir_linfocitos(int cantidad){

	double radio_del_tumor_b = -9e9; // 250.0;
	double radio_temporal_b = 0.0;


	for( unsigned int i=0; i < todas_las_celulas.size() ; i++ )
	{
		radio_temporal_b = norm_squared( todas_las_celulas[i]->posicion );
		if( radio_temporal_b > radio_del_tumor_b )
		{ radio_del_tumor_b = radio_temporal_b; }
	}

	radio_del_tumor_b = sqrt( radio_del_tumor_b );



	double radio_interno = radio_del_tumor_b + 30.0;

	double radio_externo = radio_interno + 75.0;

	double radio_medio = 0.5*(radio_interno + radio_externo);
	double radio_std = 0.33*( radio_externo-radio_interno)/2.0;
    Celula* pLinfocito = NULL;

	for( int i=0 ;i < cantidad ; i++ )
	{
		double theta = rng->RandomNumber() * 6.283185307179586476925286766559;
		double phi = acos( 2.0*rng->RandomNumber() - 1.0 );

		double radio = rng->NormalRandom_CM( radio_medio, radio_std );


		pLinfocito = new Linfocito;
		todas_las_celulas.push_back(pLinfocito);
		celulas_para_registrar_en_voxeles.push_back(pLinfocito);
		pLinfocito->indice=todas_las_celulas.size()-1;
		//pLinfocito->id=todas_las_celulas.size();
    	//pLinfocito->inicializar_linfocito();
    	pLinfocito->set_posicion(radio*cos(theta)*sin(phi), radio*sin(theta)*sin(phi), radio*cos(phi));

    	cdc.registrar_celula(pLinfocito);
    	celulas_para_registrar_en_voxeles.clear();
	}

	return;



}


/**
 * @brief Introduces lymphocytes (immune cells) randomly throughout the domain
 * 
 * @param cantidad Number of lymphocytes to introduce
 * 
 * Cancer Research Context:
 * This method models a different pattern of immune cell distribution by randomly
 * placing lymphocytes throughout the simulation domain, but outside the tumor.
 * This approach is valuable for studying:
 * - Stochastic immune surveillance mechanisms
 * - Heterogeneous immune infiltration patterns
 * - Immune cell recruitment from distant sites
 * - Differences between localized and systemic immune responses
 * 
 * The random distribution creates a more realistic representation of the
 * variability in immune cell positioning observed in some tumor types,
 * particularly those with less structured immune infiltration patterns.
 */
void Tejido::introducir_linfocitos_aleatorios(int cantidad){

	double radio_del_tumor_b = -9e9; // 250.0;
	double radio_temporal_b = 0.0;
    Vector pos_temp;
    double radio_temp_2 = 0.0;


	for( unsigned int i=0; i < todas_las_celulas.size() ; i++ )
	{
        if( todas_las_celulas[i]->tipo == 0){
            radio_temporal_b = norm_squared( todas_las_celulas[i]->posicion );
            if( radio_temporal_b > radio_del_tumor_b )
            { radio_del_tumor_b = radio_temporal_b; }
        }
    }   

	radio_del_tumor_b = sqrt( radio_del_tumor_b );
	double radio_interno = radio_del_tumor_b + 50.0;

    


    Celula* pLinfocito = NULL;

	for( int i=0 ;i < cantidad ; i++ )
	{

        do{
            pos_temp.x = rng->RandomNumber(pg->rango_en_X[0],pg->rango_en_X[1]);
            pos_temp.y = rng->RandomNumber(pg->rango_en_Y[0],pg->rango_en_Y[1]);
            pos_temp.z = rng->RandomNumber(pg->rango_en_Z[0],pg->rango_en_Z[1]);
            radio_temp_2 = sqrt(pos_temp.x*pos_temp.x + pos_temp.y*pos_temp.y + pos_temp.z*pos_temp.z);
        }while(radio_interno > radio_temp_2);
        
        
        
        if(radio_interno < radio_temp_2){
            pLinfocito = new Linfocito;
            todas_las_celulas.push_back(pLinfocito);
            celulas_para_registrar_en_voxeles.push_back(pLinfocito);
            pLinfocito->indice=todas_las_celulas.size()-1;
            pLinfocito->set_posicion(pos_temp.x, pos_temp.y, pos_temp.z);
            cdc.registrar_celula(pLinfocito);
            celulas_para_registrar_en_voxeles.clear();
        }
	}

	return;



}

/**
 * @brief Calculates and updates tumor geometry metrics
 * 
 * Cancer Research Context:
 * This method analyzes the current state of the tumor to calculate key metrics
 * that are essential for quantitative cancer research:
 * - Tumor radius: Indicates the extent of tumor spread and invasion
 * - Tumor volume: Critical measure for assessing tumor growth kinetics
 * - Cell counts: Quantifies tumor cellularity and death rates
 * 
 * These metrics enable quantitative analysis of:
 * - Tumor growth rates under different conditions
 * - Effects of treatments on tumor size and composition
 * - Correlation between tumor geometry and microenvironmental factors
 * - Validation of computational models against experimental data
 * 
 * The method uses efficient computational approaches (calculating norm squared
 * before taking square root) to optimize performance for large cell populations.
 */
void Tejido::geometria_del_tumor(){
    
 	radio_del_tumor = -9e9; // 250.0;
	double radio_temporal = 0.0;
    volumen_del_tumor = 0.0;
    volumen_del_tumor2 = 0.0;
    celulas_tumorales = 0;
    celulas_muertas = 0;
    static double cuatrotercios = 1.33333333333333333333;
    static double pi = 3.1415926535897932384626433832795;


	for( unsigned int i=0; i < todas_las_celulas.size() ; i++ )
	{
        if( todas_las_celulas[i]->tipo == 0){
            radio_temporal = norm_squared( todas_las_celulas[i]->posicion );
            volumen_del_tumor2 += todas_las_celulas[i]->fenotipo.volumen.total;
            if( radio_temporal > radio_del_tumor ){
                radio_del_tumor = radio_temporal;                
            }
            celulas_tumorales += 1;
        }
        if(todas_las_celulas[i]->fenotipo.muerte.muerta){
            celulas_muertas +=1;
        }
	}

	radio_del_tumor = sqrt( radio_del_tumor );
    
    volumen_del_tumor = cuatrotercios * pi * radio_del_tumor * radio_del_tumor * radio_del_tumor;
    
    return;
    
}
