/*! @file Ciclos_estandares.cpp
 *
 * @brief Implementation of standard cell cycle models used in cancer simulation.
 *
 * @author Luciana Melina Luque
 *
 * @details
 * This file implements several standard cell cycle models that are fundamental
 * to cancer simulation: Ki67 (based on this proliferation marker), basic life cycle,
 * necrosis (cell death due to nutrient deprivation/hypoxia), and apoptosis
 * (programmed cell death). These models allow simulation of tumor growth dynamics
 * and response to microenvironmental conditions.
 *
 * Inbound Dependencies:
 *   - Ciclos_estandares.h - Function declarations
 *   - Volumen.h - Volume calculations for cells
 *   - Constantes.h - Phase constants
 *
 * Outbound Dependencies:
 *   - None
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */


#include "Ciclos_estandares.h"

/*! 
 * \brief Implementation of entry function for Ki67 positive phase
 *
 * When a cell enters the Ki67 positive phase (indicating proliferation),
 * this function doubles the target volumes for the cell's nucleus and cytoplasm.
 * This simulates the growth process that occurs before cell division,
 * which is a critical aspect of modeling tumor growth rates.
 *
 * \param volumen The cell's volume parameters to be modified
 * \param mp Death parameters (unused in this function)
 */
void Ki67_fase_positiva_funcion_de_entrada(Volumen& volumen, Muerte_parametros& mp)
{

	volumen.target_nucleo_solido *= 2.0;
	volumen.target_citoplasma_solido *= 2.0;

	return;
}

/*!
 * \brief Implementation of entry function for basic life cycle phase
 *
 * Similar to the Ki67 positive phase, this function doubles volume targets
 * when a cell enters the basic life phase. This is used for simpler
 * cell cycle models where we don't distinguish between different proliferative phases.
 *
 * \param volumen The cell's volume parameters to be modified
 * \param mp Death parameters (unused in this function)
 */
void fase_viva_funcion_de_entrada(Volumen& volumen, Muerte_parametros& mp){


	volumen.target_nucleo_solido *= 2.0;
	volumen.target_citoplasma_solido *= 2.0;

	return;

}

/*!
 * \brief Implementation of entry function for apoptosis phase
 *
 * When a cell undergoes apoptosis (programmed cell death), this function
 * sets volume targets to simulate cellular shrinkage. Apoptosis is a critical
 * process in cancer research as it's often dysregulated in tumors and is
 * a target for many cancer therapies. The function also sets rate parameters
 * for the volume changes based on death parameters.
 *
 * \param volumen The cell's volume parameters to be modified
 * \param mp Death parameters that control the rate of volume changes
 */
void standard_apoptosis_funcion_de_entrada( Volumen& volumen, Muerte_parametros& mp )
{

	volumen.target_fraccion_fluido = 0.0;
	volumen.target_citoplasma_solido = 0.0;
	volumen.target_nucleo_solido = 0.0;

	volumen.target_relacion_citoplasma_nucleo = 0.0;



	volumen.citoplasma_tasa_de_cambio = mp.citoplasma_tasa_de_cambio;

	volumen.nucleo_tasa_de_cambio = mp.nucleo_tasa_de_cambio;

	volumen.fluido_tasa_de_cambio = mp.tasa_de_cambio_fluido_lisado;

	volumen.tasa_de_calcificacion = mp.tasa_de_calcificacion;

	volumen.volumen_de_ruptura_relativo = mp.volumen_de_ruptura_relativo;
	volumen.volumen_de_ruptura = volumen.total * volumen.volumen_de_ruptura_relativo;

	return;
}

/*!
 * \brief Implementation of entry function for necrosis phase
 *
 * When a cell becomes necrotic (typically due to hypoxia or nutrient deprivation),
 * this function sets volume targets to simulate cellular swelling that occurs
 * during necrosis. In tumors, necrosis commonly occurs in poorly vascularized regions
 * and is an important indicator of aggressive cancer. The function also sets
 * various rate parameters for the volume changes.
 *
 * \param volumen The cell's volume parameters to be modified
 * \param mp Death parameters that control the rate of volume changes
 */
void standard_necrosis_funcion_de_entrada( Volumen& volumen, Muerte_parametros& mp ){

	volumen.target_fraccion_fluido = 1.0;
	volumen.target_citoplasma_solido = 0.0;
	volumen.target_nucleo_solido = 0.0;

volumen.target_relacion_citoplasma_nucleo = 0.0;


	volumen.citoplasma_tasa_de_cambio = mp.citoplasma_tasa_de_cambio;

	volumen.nucleo_tasa_de_cambio = mp.nucleo_tasa_de_cambio;

	volumen.fluido_tasa_de_cambio = mp.tasa_de_cambio_fluido_lisado;

	volumen.tasa_de_calcificacion = mp.tasa_de_calcificacion;


	volumen.volumen_de_ruptura_relativo = mp.volumen_de_ruptura_relativo;

	volumen.volumen_de_ruptura =
	volumen.total * volumen.volumen_de_ruptura_relativo;

	return;

}

/*!
 * \brief Implementation of entry function for lysis phase
 *
 * When a cell enters the lysis phase (final breakdown after death),
 * this function sets volume targets to simulate complete cellular
 * degradation. This phase represents the final fate of dead cells
 * within the tumor microenvironment before they are removed from the simulation.
 *
 * \param volumen The cell's volume parameters to be modified
 * \param mp Death parameters that control the rate of volume changes
 */
void standard_lysis_funcion_de_entrada( Volumen& volumen, Muerte_parametros& mp ){


	volumen.target_fraccion_fluido = 0.0;
	volumen.target_citoplasma_solido = 0.0;
	volumen.target_nucleo_solido = 0.0;


	volumen.citoplasma_tasa_de_cambio = mp.citoplasma_tasa_de_cambio;

	volumen.nucleo_tasa_de_cambio = mp.nucleo_tasa_de_cambio;

	volumen.fluido_tasa_de_cambio = mp.tasa_de_cambio_fluido_lisado;

	volumen.tasa_de_calcificacion = mp.tasa_de_calcificacion;


	volumen.volumen_de_ruptura_relativo = 9e99;
	volumen.volumen_de_ruptura =
	volumen.total * volumen.volumen_de_ruptura_relativo;

	return;

}

/*!
 * \brief Determines if a necrotic cell should transition to the lysed phase
 *
 * This function checks if a necrotic cell has reached the volume threshold for
 * bursting/rupture. In necrotic cell death, cells typically swell until they rupture,
 * which this function models. This process is relevant for understanding tumor
 * necrotic core formation and related inflammatory responses.
 *
 * \param volumen The cell's current volume parameters
 * \param mp Death parameters (unused in this function)
 * \return true if the cell should remain in current phase, false if it should transition
 */
bool standard_necrosis_funcion_de_arrest( Volumen& volumen, Muerte_parametros& mp ){

	if( volumen.total < volumen.volumen_de_ruptura ){

		return true;

	}

	return false;

}

/*!
 * \brief Creates and initializes the Ki67-based cell cycle model
 *
 * This function initializes the Ki67 cell cycle model with three phases:
 * Ki67 negative, Ki67 positive pre-mitotic, and Ki67 positive post-mitotic.
 * The model is based on the Ki67 proliferation marker, which is widely used
 * in cancer research to assess tumor growth rates and proliferation indices.
 * It sets phase transition rates based on experimental data from cancer cell lines.
 */
void crear_ciclo_ki67( void )
{

	Ki67.codigo = Constantes::ciclo_Ki67;
	Ki67.nombre="ki67";
	Ki67.unidades_tiempo ="min";

	Ki67.agregar_fase( Constantes::Ki67_negativa, "ki67-" );
	Ki67.agregar_fase( Constantes::Ki67_positiva_premitotica, "ki67+ (premitotica)");
	Ki67.agregar_fase( Constantes::Ki67_positiva_postmitotica, "ki67+ (postmitotica)");


	Ki67.fases[1].division_al_final_de_la_fase = true;


	Ki67.fases[1].actualizar_volumen = true;
	Ki67.fases[2].actualizar_volumen = true;


	Ki67.agregar_link( 0 , 1, NULL ); // - to + (pre-mitotic)
	Ki67.agregar_link( 1 , 2, NULL ); // + (pre-mitotic) to + (post-mitotic)
	Ki67.agregar_link( 2 , 0, NULL ); // + (post-mitotic) to -


	Ki67.fase_link(1,2).duracion_fija = true;
	Ki67.fase_link(2,0).duracion_fija = true;
	//Ki67.fase_link(0,1).duracion_fija = true; //Volver a false cuando termine de probar

	// Set transition rates based on experimental data
	Ki67.tasa_de_transicion(0,1) = 1.0/(3.62*60.0);  // Rate from - to + (pre-mitotic)
	Ki67.tasa_de_transicion(1,2) = 1.0/(13.0*60.0);  // Rate from + (pre-mitotic) to + (post-mitotic)
	Ki67.tasa_de_transicion(2,0) = 1.0/(2.5*60.0);   // Rate from + (post-mitotic) to -


	Ki67.fases[1].funcion_de_entrada = Ki67_fase_positiva_funcion_de_entrada;

	return;
}


/*!
 * \brief Creates and initializes the basic life cycle model
 *
 * This function initializes a simplified cell cycle model with just one phase
 * (living). This basic model is useful for simulations where detailed cell cycle
 * phases are not needed. The model's division rate is calibrated to match
 * the net birth rate of MCF10A cells, which are often used as normal breast
 * epithelial cell controls in breast cancer research.
 */
void crear_ciclo_vida( void )
{
	vida.codigo = Constantes::ciclo_vida;
	vida.nombre = "Vida";
	vida.unidades_tiempo = "min";

	vida.agregar_fase( Constantes::viva , "Viva" );

	vida.fases[0].actualizar_volumen = true;

	vida.fases[0].division_al_final_de_la_fase = true;

	vida.agregar_link( 0 , 0 , NULL );
    
    vida.fase_link(0,0).duracion_fija = true;    
    //vida.tasa_de_transicion(0,0) = 1 /(0.24*60.0);


    vida.tasa_de_transicion(0,0) = 0.02717 / 60.0;
    
	vida.fases[0].funcion_de_entrada = fase_viva_funcion_de_entrada;

	return;
}


/*!
 * \brief Creates and initializes the necrosis death cycle model
 *
 * This function initializes the necrosis model with its parameters and phases.
 * Necrosis is a type of cell death that occurs when cells are deprived of
 * oxygen or nutrients, which is common in the center of rapidly growing tumors.
 * The model includes a swelling phase followed by lysis, with parameters
 * calibrated to match experimental observations of necrotic cell death.
 */
void crear_ciclo_necrosis(void){

	necrosis_parametros.tiempo_unidades = "min";


	necrosis_parametros.citoplasma_tasa_de_cambio = 0.0032 /60.0;
	necrosis_parametros.nucleo_tasa_de_cambio = 0.013 / 60.0;

	necrosis_parametros.tasa_de_cambio_fluido_no_lisado = 0.67/ 60.0;
	necrosis_parametros.tasa_de_cambio_fluido_lisado = 0.050 /60.0;

	necrosis_parametros.tasa_de_calcificacion = 0.0042 / 60.0;

	necrosis_parametros.volumen_de_ruptura_relativo = 2.0;


	necrosis.nombre = "Necrosis";
	necrosis.codigo = Constantes::ciclo_de_muerte_necrosis;


	necrosis.agregar_fase(Constantes::necrotica_hinchada, "Necrotica (swelling)");
	necrosis.fases[0].funcion_de_entrada = standard_necrosis_funcion_de_entrada;


	necrosis.agregar_fase(Constantes::necrotica_lisada, "Necrotica (lysed)");
	necrosis.fases[1].funcion_de_entrada = standard_lysis_funcion_de_entrada;
	necrosis.fases[1].remover_al_final_de_la_fase = true;



	necrosis.agregar_fase(Constantes::debris, "Debris");

	necrosis.fases[0].actualizar_volumen = true; 
	necrosis.fases[1].actualizar_volumen = true;


	necrosis.agregar_link(0,1,standard_necrosis_funcion_de_arrest);
	necrosis.agregar_link(1,2,NULL);



	necrosis.tasa_de_transicion(0,1) = 9e9;
	necrosis.tasa_de_transicion(1,2) = 1.0/(60.0 * 24.0 * 60.0); // 1 day

	necrosis.fase_link(1,2).duracion_fija = true;

	return;

}


/*!
 * \brief Creates and initializes the apoptosis death cycle model
 *
 * This function initializes the apoptosis model with its parameters and phases.
 * Apoptosis is programmed cell death, characterized by cell shrinkage and fragmentation.
 * It's a crucial process in cancer biology, as many cancers develop mechanisms to
 * evade apoptosis, and many therapies aim to induce it. The model parameters are
 * calibrated to match experimental observations of apoptotic cell death.
 */
void crear_ciclo_apoptosis(void){

	apoptosis_parametros.tiempo_unidades = "min";


	apoptosis_parametros.citoplasma_tasa_de_cambio = 1.0 / 60.0;
	apoptosis_parametros.nucleo_tasa_de_cambio = 0.35 / 60.0;

	apoptosis_parametros.tasa_de_cambio_fluido_no_lisado = 3.0 / 60.0;
	apoptosis_parametros.tasa_de_cambio_fluido_lisado = 0.0;

	apoptosis_parametros.tasa_de_calcificacion = 0.0;

	apoptosis_parametros.volumen_de_ruptura_relativo = 2.0;


	apoptosis.nombre = "Apoptosis";
	apoptosis.codigo = Constantes::ciclo_de_muerte_apoptosis;




	apoptosis.agregar_fase( Constantes::apoptotica , "Apoptotica" );
	apoptosis.fases[0].funcion_de_entrada = standard_apoptosis_funcion_de_entrada;
	apoptosis.fases[0].remover_al_final_de_la_fase = true;


	apoptosis.agregar_fase( Constantes::debris, "Debris");



	apoptosis.agregar_link( 0, 1, NULL );
	apoptosis.tasa_de_transicion( 0, 1) = 1.0 / (8.6 * 60.0);  // 8.6 hours


	apoptosis.fase_link(0,1).duracion_fija = true;

	return;

}

/*!
 * \brief Creates and initializes all standard cell cycles
 *
 * This function creates and initializes the standard proliferative cell cycles
 * (Ki67 and basic life cycle). It only performs initialization if it hasn't been
 * done already. These cell cycle models are essential for simulating
 * tumor growth dynamics in cancer research.
 *
 * \return true if initialization was successful, false if already initialized
 */
bool crear_ciclo_celular_estandar( void )
{
	if( ciclo_celular_estandar_inicializado ){

		return false;
	}else{

		crear_ciclo_ki67();
		crear_ciclo_vida();
		return true;
	}
}

/*!
 * \brief Creates and initializes all standard death cycles
 *
 * This function creates and initializes the standard cell death cycles
 * (necrosis and apoptosis). It only performs initialization if it hasn't been
 * done already. These death cycle models are crucial for simulating
 * cell death processes in cancer research, including effects of hypoxia and
 * response to therapies.
 *
 * \return true if initialization was successful, false if already initialized
 */
bool crear_ciclo_de_muerte_estandar( void )
{

	if( ciclo_celular_de_muerte_inicializado){

		return false;
	}else{

		crear_ciclo_necrosis();
		crear_ciclo_apoptosis();
		return true;
	}

}
