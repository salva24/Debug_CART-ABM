/**
 * @file: Microambiente.cpp
 *
 * @author: Luciana Melina Luque
 *
 * @description: Implementation of the microenvironment for cancer simulation.
 * This file implements the multi-substrate diffusive microenvironment that models
 * biochemical factors like oxygen, nutrients, and signaling molecules using a
 * finite volume method on a Cartesian grid.
 * 
 * Cancer Research Context:
 * The tumor microenvironment is a critical aspect of cancer biology. This implementation
 * models diffusion of oxygen and other substrates, creating gradients that affect
 * cancer cell behavior. Key cancer-relevant aspects include:
 * - Hypoxia (low oxygen) regions that drive cancer cell adaptation and resistance
 * - Nutrient gradients that influence cancer cell proliferation and metabolism
 * - Signaling molecule distributions that affect cancer cell phenotypes
 * - Blood vessel modeling that simulates tumor angiogenesis and vascularization
 * 
 * Inbound Dependencies:
 * - Microambiente.h - Class definition and method declarations
 * - Grillado.h - Cartesian grid implementation for spatial discretization
 * - Parametros_globales.h - Global simulation parameters
 * 
 * Outbound Dependencies:
 * - Celula.h - Cells interact with the microenvironment
 * - Tejido.h - Tissue uses the microenvironment for cellular organization
 * 
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#include "Microambiente.h"


/**
 * @brief Global pointer to the default microenvironment instance
 * 
 * This variable stores a pointer to the default microenvironment, which
 * can be accessed throughout the simulation for convenience.
 * 
 * Cancer Research Context: 
 * Having a globally accessible microenvironment enables consistent 
 * modeling of the tumor microenvironment across different components.
 */
Microambiente* microambiente_default = NULL;

/**
 * @brief Sets the default microenvironment instance
 * 
 * @param M Pointer to the microenvironment to set as default
 * 
 * Cancer Research Context:
 * Establishes a single microenvironment that maintains consistency
 * across the simulation of tumor and its surrounding conditions.
 */
void set_microambiente_default( Microambiente* M ){

	microambiente_default = M;

}

/**
 * @brief Returns a pointer to the default microenvironment
 * 
 * @return Microambiente* Pointer to the default microenvironment instance
 * 
 * Cancer Research Context:
 * Provides global access to the shared tumor microenvironment, ensuring
 * that all cells interact with the same biochemical conditions.
 */
Microambiente* get_microambiente_default( void ){

	return microambiente_default;

}

/**
 * @brief Constructor for the Microambiente class
 * 
 * Initializes a new microenvironment with default properties, including
 * substrate vectors, diffusion coefficients, decay rates, and boundary 
 * conditions. Creates a minimal 1x1x1 grid initially, which will typically
 * be resized before simulation.
 * 
 * Cancer Research Context:
 * Creates the foundation for modeling biochemical conditions within and
 * around a tumor, which significantly influence cancer progression.
 */
Microambiente::Microambiente(){

	nombre = "unnamed";
	unidades_espaciales = "none";
	unidades_temporales = "none";

	setup_del_solver_bulk_fuente_sumidero_hecho = false;
	setup_de_thomas_hecho = false;
	setup_del_solver_de_difusion_hecho = false;




	mgrilla.redimensionar(1,1,1);

	uno.resize( 1 , 1.0 );
	cero.resize( 1 , 0.0 );

	temp_vectores_densidad1.resize( mgrilla.voxeles.size() , cero );
	temp_vectores_densidad2.resize( mgrilla.voxeles.size() , cero );
	p_vectores_densidad = &temp_vectores_densidad1;

	vectores_gradiente.resize( mgrilla.voxeles.size() );
	for( unsigned int k=0 ; k < mgrilla.voxeles.size() ; k++ )
	{
		vectores_gradiente[k].resize( 1 );
		(vectores_gradiente[k])[0].resize( 3, 0.0 );
	}
	vector_gradiente_calculado.resize( mgrilla.voxeles.size() , false );





	densidades_nombres.assign( 1 , "unnamed" );
	densidades_unidades.assign( 1 , "none" );

	coeficientes_de_difusion.assign( numero_de_densidades() , 0.0 );
	tasas_de_decaimiento.assign( numero_de_densidades() , 0.0 );

	un_medio = uno;
	un_medio *= 0.5;

	un_tercio = uno;
	un_tercio /= 3.0;

	vector_valores_de_dirichlet.assign( mgrilla.voxeles.size(), uno );
	vector_activacion_dirichlet.assign( 1 , true );

	dirichlet_vectorES_activacion.assign( 1 , vector_activacion_dirichlet );
/*
	pg->dirichlet_todo.assign( 1 , true );
	pg->dirichlet_xmin.assign( 1 , true );
	pg->dirichlet_xmax.assign( 1 , true );
	pg->dirichlet_ymin.assign( 1 , true );
	pg->dirichlet_ymax.assign( 1 , true );
	pg->dirichlet_zmin.assign( 1 , true );
	pg->dirichlet_zmax.assign( 1 , true );
	pg->dirichlet_vs.assign( 1 , true );
*/
    pg->dirichlet_xmin_valores.assign( 1 , 1.0 );
	pg->dirichlet_xmax_valores.assign( 1 , 1.0 );
	pg->dirichlet_ymin_valores.assign( 1 , 1.0 );
	pg->dirichlet_ymax_valores.assign( 1 , 1.0 );
	pg->dirichlet_zmin_valores.assign( 1 , 1.0 );
	pg->dirichlet_zmax_valores.assign( 1 , 1.0 );

	if(microambiente_default==NULL){
		microambiente_default=this;
	}

	return;

}

unsigned int Microambiente::numero_de_densidades( void ){

	return (*p_vectores_densidad)[0].size();
}

unsigned int Microambiente::numero_de_voxeles( void ){

	return mgrilla.voxeles.size();
}

/**
 * @brief Resizes the microenvironment's spatial domain
 * 
 * @param x_ini Initial x-coordinate of the domain
 * @param x_fin Final x-coordinate of the domain
 * @param y_ini Initial y-coordinate of the domain
 * @param y_fin Final y-coordinate of the domain
 * @param z_ini Initial z-coordinate of the domain
 * @param z_fin Final z-coordinate of the domain
 * @param dx_nuevo New voxel size in x-direction
 * @param dy_nuevo New voxel size in y-direction
 * @param dz_nuevo New voxel size in z-direction
 * 
 * This method resizes the spatial domain of the microenvironment and 
 * reinitializes all data structures that depend on the spatial discretization.
 * 
 * Cancer Research Context:
 * The spatial domain configuration is crucial for simulating realistic tumor 
 * microenvironments. Proper dimensioning allows for accurate representation of:
 * - Tumor boundaries and invasion zones
 * - Heterogeneous regions with varying substrate concentrations
 * - Local microenvironmental niches that drive cancer cell adaptations
 * The spatial resolution (voxel size) directly impacts the fidelity of diffusion 
 * gradient modeling, which is essential for simulating hypoxia-driven malignant 
 * progression.
 */
void Microambiente::redimensionar_espacio( double x_ini, double x_fin, double y_ini, double y_fin, double z_ini, double z_fin , double dx_nuevo , double dy_nuevo , double dz_nuevo ){

	mgrilla.redimensionar(x_ini, x_fin, y_ini, y_fin, z_ini, z_fin , dx_nuevo, dy_nuevo, dz_nuevo );

	temp_vectores_densidad1.assign( mgrilla.voxeles.size() , cero );
	temp_vectores_densidad2.assign( mgrilla.voxeles.size() , cero );

	vectores_gradiente.resize( mgrilla.voxeles.size() );
	for( unsigned int k=0 ; k < mgrilla.voxeles.size() ; k++ )
	{
		vectores_gradiente[k].resize( numero_de_densidades() );
		for( unsigned int i=0 ; i < numero_de_densidades() ; i++ )
		{
			(vectores_gradiente[k])[i].resize( 3, 0.0 );
		}
	}
	vector_gradiente_calculado.resize( mgrilla.voxeles.size() , false );

	vector_valores_de_dirichlet.assign( mgrilla.voxeles.size(), uno );

	dirichlet_vectorES_activacion.assign( mgrilla.voxeles.size() , vector_activacion_dirichlet );

	return;
}

/**
 * @brief Resizes the density vectors for all substrates in the microenvironment
 * 
 * @param nuevo_tamano New size for density vectors
 * 
 * This method resizes all data structures related to substrate densities when
 * the number of substrates changes.
 * 
 * Cancer Research Context:
 * The ability to model multiple substrates simultaneously is critical for cancer research
 * as it allows simulation of:
 * - Oxygen and nutrient gradients that drive metabolic adaptations in cancer cells
 * - Signaling molecules that affect cancer cell phenotype and behavior
 * - Drug distributions for treatment response modeling
 * - Waste product accumulation that creates acidic microenvironments
 * Each substrate can have unique diffusion and decay properties, enabling realistic
 * modeling of the complex biochemical landscape that influences cancer progression.
 */
void Microambiente::redimensionar_densidades( int nuevo_tamano ){

	cero.assign( nuevo_tamano, 0.0 );
	uno.assign( nuevo_tamano, 1.0 );

	temp_vectores_densidad1.assign( mgrilla.voxeles.size() , cero );
	temp_vectores_densidad2.assign( mgrilla.voxeles.size() , cero );

	for( unsigned int k=0 ; k < mgrilla.voxeles.size() ; k++ )
	{
		vectores_gradiente[k].resize( numero_de_densidades() );
		for( unsigned int i=0 ; i < numero_de_densidades() ; i++ )
		{
			(vectores_gradiente[k])[i].resize( 3, 0.0 );
		}
	}
	vector_gradiente_calculado.resize( mgrilla.voxeles.size() , false );

	coeficientes_de_difusion.assign( nuevo_tamano , 0.0 );
	tasas_de_decaimiento.assign( nuevo_tamano , 0.0 );

	densidades_nombres.assign( nuevo_tamano, "unnamed" );
	densidades_unidades.assign( nuevo_tamano , "none" );

	un_medio = uno;
	un_medio *= 0.5;

	un_tercio = uno;
	un_tercio /= 3.0;

	vector_valores_de_dirichlet.assign( mgrilla.voxeles.size(), uno );
	vector_activacion_dirichlet.assign( nuevo_tamano, true );

	dirichlet_vectorES_activacion.assign(mgrilla.voxeles.size(), vector_activacion_dirichlet );

	pg->vector_condicion_de_dirichlet.assign( nuevo_tamano , 1.0 );
	pg->vector_activacion_dirichlet.assign( nuevo_tamano, true );

	pg->vector_condiciones_iniciales.assign( nuevo_tamano , 1.0 );

	pg->dirichlet_todo.assign( nuevo_tamano , true );

	pg->dirichlet_xmin.assign( nuevo_tamano , false );
	pg->dirichlet_xmax.assign( nuevo_tamano , false );
	pg->dirichlet_ymin.assign( nuevo_tamano , false );
	pg->dirichlet_ymax.assign( nuevo_tamano , false );
	pg->dirichlet_zmin.assign( nuevo_tamano , false );
	pg->dirichlet_zmax.assign( nuevo_tamano , false );
	pg->dirichlet_vs.assign( nuevo_tamano , false );

	pg->dirichlet_xmin_valores.assign( nuevo_tamano , 1.0 );
	pg->dirichlet_xmax_valores.assign( nuevo_tamano , 1.0 );
	pg->dirichlet_ymin_valores.assign( nuevo_tamano , 1.0 );
	pg->dirichlet_ymax_valores.assign( nuevo_tamano , 1.0 );
	pg->dirichlet_zmin_valores.assign( nuevo_tamano , 1.0 );
	pg->dirichlet_zmax_valores.assign( nuevo_tamano , 1.0 );

	return;


}

/**
 * @brief Adds a new substrate to the microenvironment with specified properties
 * 
 * @param nombre Name of the substrate
 * @param unidades Units of measurement for the substrate
 * @param coeficiente_de_difusion Diffusion coefficient of the substrate
 * @param tasa_de_decaimiento Decay rate of the substrate
 * 
 * Cancer Research Context:
 * Different biochemical factors have unique diffusion and decay properties that
 * affect their distribution within tumors. For example, oxygen has a high diffusion
 * coefficient but is rapidly consumed, creating steep gradients in tumors. Growth
 * factors may diffuse more slowly but persist longer. These properties directly
 * influence tumor growth patterns, invasion, and responses to therapy.
 */
void Microambiente::agregar_densidad(std::string nombre , std::string unidades, double coeficiente_de_difusion, double tasa_de_decaimiento ){

	// update 1, 0
	cero.push_back( 0.0 );
	uno.push_back( 1.0 );

	// update units
	densidades_nombres.push_back( nombre );
	densidades_unidades.push_back( unidades );

	// update coefficients
	coeficientes_de_difusion.push_back( coeficiente_de_difusion );
	tasas_de_decaimiento.push_back( tasa_de_decaimiento );

	// update sources and such
	for( unsigned int i=0; i < temp_vectores_densidad1.size() ; i++ )
	{
		temp_vectores_densidad1[i].push_back( 0.0 );
		temp_vectores_densidad2[i].push_back( 0.0 );
	}

	// resize the gradient data structures
	for( unsigned int k=0 ; k < mgrilla.voxeles.size() ; k++ )
	{
		vectores_gradiente[k].resize( numero_de_densidades() );
		for( unsigned int i=0 ; i < numero_de_densidades(); i++ )
		{
			(vectores_gradiente[k])[i].resize( 3, 0.0 );
		}
	}
	vector_gradiente_calculado.resize( mgrilla.voxeles.size() , false );

	un_medio = uno;
	un_medio *= 0.5;


	un_tercio = uno;
	un_tercio /= 3.0;

	vector_valores_de_dirichlet.assign( mgrilla.voxeles.size(), uno );
	vector_activacion_dirichlet.push_back(true );
	dirichlet_vectorES_activacion.assign(mgrilla.voxeles.size(), vector_activacion_dirichlet );

	pg->vector_condicion_de_dirichlet.push_back( 1.0 );
	pg->vector_activacion_dirichlet.push_back( true );

	pg->vector_condiciones_iniciales.push_back( 1.0 );

	pg->dirichlet_todo.push_back( false );

	pg->dirichlet_xmin.push_back( false );
	pg->dirichlet_xmax.push_back( false );
	pg->dirichlet_ymin.push_back( false );
	pg->dirichlet_ymax.push_back( false );
	pg->dirichlet_zmin.push_back( false );
	pg->dirichlet_zmax.push_back( false );
	pg->dirichlet_vs.push_back( false );

	pg->dirichlet_xmin_valores.push_back( 1.0 );
	pg->dirichlet_xmax_valores.push_back( 1.0 );
	pg->dirichlet_ymin_valores.push_back( 1.0 );
	pg->dirichlet_ymax_valores.push_back( 1.0 );
	pg->dirichlet_zmin_valores.push_back( 1.0 );
	pg->dirichlet_zmax_valores.push_back( 1.0 );

	return;


}

/**
 * @brief Sets properties for an existing substrate at the specified index
 * 
 * @param indice Index of the substrate to modify
 * @param nombre Name of the substrate
 * @param unidades Units of measurement for the substrate
 * @param coeficiente_de_difusion Diffusion coefficient of the substrate
 * @param tasa_de_decaimiento Decay rate of the substrate
 * 
 * Cancer Research Context:
 * This method allows configuring specific substrates with properties that match
 * experimental measurements. For example, setting oxygen with appropriate diffusion
 * coefficients (1e5 μm²/min) and consumption rates (0.1/min) creates realistic
 * oxygen gradients observed in tumors. These correctly calibrated parameters
 * are essential for accurately modeling hypoxic regions and their effects on
 * cancer cell behavior.
 */
void Microambiente::set_densidad(int indice, std::string nombre , std::string unidades , double coeficiente_de_difusion , double tasa_de_decaimiento ){


	if( indice == 0 )
	{ pg->usar_oxigeno_como_primer_sustrato = false; }

	densidades_nombres[indice] = nombre;
	densidades_unidades[indice] = unidades;

	coeficientes_de_difusion[indice] = coeficiente_de_difusion;
	tasas_de_decaimiento[indice] = tasa_de_decaimiento;
	return;


}

void Microambiente::set_densidad( int indice , std::string nombre , std::string unidades ){

	if( indice == 0 )
	{ pg->usar_oxigeno_como_primer_sustrato = false; }

	densidades_nombres[indice] = nombre;
	densidades_unidades[indice] = unidades;
}

int Microambiente::encontrar_indice_de_densidad( std::string nombre ){

	for( unsigned int i=0; i < densidades_nombres.size() ; i++ )
	{
		if( densidades_nombres[i] == nombre )
		{ return i; }
	}
	return -1;
}

/**
 * @brief Converts 3D cartesian indices to a linear voxel index
 * 
 * @param i x-axis index
 * @param j y-axis index
 * @param k z-axis index
 * @return int Linear voxel index
 * 
 * This method maps 3D coordinates in the grid to a single linear index
 * for accessing the voxel data.
 * 
 * Cancer Research Context:
 * Efficient spatial indexing is essential for modeling how cancer cells interact
 * with their local environment. This method facilitates:
 * - Tracking cell positions relative to nutrient/oxygen sources
 * - Modeling cell-cell and cell-matrix interactions in 3D space
 * - Analyzing spatial patterns of cancer growth and invasion
 * - Simulating spatially-targeted therapeutic interventions
 */
int Microambiente::indice_de_voxel( int i, int j, int k ){

	return mgrilla.indice_de_voxel(i,j,k) ;
}

/**
 * @brief Gets the center position of a voxel
 * 
 * @param indice_del_voxel Linear index of the voxel
 * @return Vector 3D vector representing the center position
 * 
 * This method returns the spatial coordinates of a voxel's center
 * given its linear index.
 * 
 * Cancer Research Context:
 * Accurate voxel positioning is crucial for:
 * - Calculating distances between cells and resources/signals
 * - Determining cellular exposure to diffusing drugs
 * - Mapping the physical landscape of the tumor microenvironment
 * - Visualizing spatial heterogeneity within the tumor
 */
Vector Microambiente::centro_del_voxel(int indice_del_voxel){

    return mgrilla.get_centro_voxel(indice_del_voxel);

}

/**
 * @brief Converts a linear voxel index to 3D cartesian indices
 * 
 * @param n Linear voxel index
 * @return std::vector<unsigned int> Vector of cartesian indices (i,j,k)
 * 
 * This method is the inverse of indice_de_voxel, converting a linear index
 * back to 3D coordinates.
 * 
 * Cancer Research Context:
 * Converting between linear and 3D indices enables:
 * - Analysis of spatial patterns in tumor growth
 * - Identification of microenvironmental niches in specific regions
 * - Mapping of metabolic or phenotypic gradients across the tumor
 * - Correlation of molecular data with spatial location in the tumor
 */
std::vector<unsigned int> Microambiente::indices_cartesianos( int n ){

	return mgrilla.indices_cartesianos( n );

}

/**
 * @brief Finds the index of the voxel closest to a given position
 * 
 * @param posicion 3D vector representing a spatial position
 * @return int Index of the nearest voxel
 * 
 * This method finds the voxel that contains or is closest to a specified
 * spatial position.
 * 
 * Cancer Research Context:
 * Finding the nearest voxel is essential for:
 * - Determining local substrate concentrations around a cancer cell
 * - Modeling how cells interact with their immediate microenvironment
 * - Simulating cellular responses to local biochemical conditions
 * - Tracking cell migration through varying environmental gradients
 */
int Microambiente::indice_del_voxel_mas_cercano( Vector& posicion ){

	return mgrilla.indice_del_voxel_mas_cercano( posicion );

}

/**
 * @brief Gets the cartesian indices of the voxel closest to a position
 * 
 * @param posicion 3D vector representing a spatial position
 * @return Vector Vector containing the nearest voxel's (i,j,k) indices
 * 
 * This method returns the 3D coordinates of the voxel closest to a given position.
 * 
 * Cancer Research Context:
 * This method supports spatial analysis of:
 * - Cancer cell positions relative to vascular structures
 * - Cell proximity to tissue boundaries or other structures
 * - Local density of cells in specific tumor regions
 * - Spatial correlation between cell phenotypes and microenvironmental features
 */
Vector Microambiente::indices_cartesianos_mas_cercanos( Vector& posicion ){

	return mgrilla.indices_cartesianos_mas_cercanos( posicion );
}

/**
 * @brief Gets the voxel object closest to a given position
 * 
 * @param posicion 3D vector representing a spatial position
 * @return Voxel& Reference to the nearest voxel
 * 
 * This method provides direct access to the voxel object nearest
 * to a specified position.
 * 
 * Cancer Research Context:
 * Direct voxel access enables:
 * - Detailed analysis of the microenvironment at specific locations
 * - Modeling cell-microenvironment interactions at cellular scale
 * - Studying how local features affect cancer cell behavior
 * - Simulating how treatments modify the local environment around cells
 */
Voxel& Microambiente::voxel_mas_cercano( Vector& posicion ){

	return mgrilla.voxel_mas_cercano( posicion );

}

/**
 * @brief Gets a voxel by its linear index
 * 
 * @param indice_de_voxel Linear index of the voxel
 * @return Voxel& Reference to the specified voxel
 * 
 * This method provides direct access to a voxel object by its linear index.
 * 
 * Cancer Research Context:
 * Accessing voxels by index allows:
 * - Systematic analysis of the entire tumor microenvironment
 * - Mapping substrate distributions throughout the simulated domain
 * - Identifying regions with specific biochemical signatures
 * - Correlating microenvironmental features with cancer cell phenotypes
 */
Voxel& Microambiente::voxeles( int indice_de_voxel ){

	return mgrilla.voxeles[indice_de_voxel];

}

/**
 * @brief Returns the density vector of the voxel closest to the given position
 * 
 * @param posicion 3D spatial position to find nearest voxel
 * @return std::vector<double>& Reference to the density vector at the nearest voxel
 * 
 * Cancer Research Context:
 * This method enables sampling of biochemical components (oxygen, nutrients, etc.) 
 * at specific spatial positions, allowing cancer cells to sense and respond to their
 * local microenvironment. This is crucial for modeling how cancer cells adapt to
 * gradients of nutrients, oxygen, and signaling molecules within the tumor.
 */
std::vector<double>& Microambiente::vector_de_densidades_mas_cercano( Vector& posicion ){

	return (*p_vectores_densidad)[ mgrilla.indice_del_voxel_mas_cercano( posicion ) ];

}

/**
 * @brief Returns the density vector at the specified voxel index
 * 
 * @param indice_de_voxel Index of the voxel to return density vector for
 * @return std::vector<double>& Reference to the density vector at the specified voxel
 * 
 * Cancer Research Context:
 * Provides direct access to substrate concentrations at specific locations within
 * the tumor microenvironment. This allows for analyzing spatial heterogeneity in
 * biochemical components that drive cancer cell behavior, metabolism, and phenotypic
 * adaptation.
 */
std::vector<double>& Microambiente::vector_de_densidades_mas_cercano( int indice_de_voxel ){

	return (*p_vectores_densidad)[ indice_de_voxel ];

}


/**
 * @brief Operator overload for direct access to density vectors by voxel index
 * 
 * @param n Index of the voxel to return density vector for
 * @return std::vector<double>& Reference to the density vector at the specified voxel
 * 
 * Cancer Research Context:
 * Provides convenient syntax for accessing substrate concentrations throughout
 * the tumor microenvironment. This operator simplifies code for analyzing spatial
 * patterns of nutrients, oxygen, and signaling molecules that affect cancer cell
 * behavior and tumor progression.
 */
std::vector<double>& Microambiente::operator()( int n ){

	return (*p_vectores_densidad)[ n ];

}

/**
 * @brief Returns the density vector at the specified voxel index
 * 
 * @param n Index of the voxel to return density vector for
 * @return std::vector<double>& Reference to the density vector at the specified voxel
 * 
 * Cancer Research Context:
 * This method provides direct access to substrate concentrations at specific locations
 * within the tumor microenvironment. Accessing density vectors is crucial for analyzing
 * how biochemical factors like oxygen, nutrients, and signaling molecules influence
 * cancer cell behavior, metabolism, and phenotypic adaptations at different spatial
 * locations. These substrate distributions drive tumor heterogeneity and influence
 * treatment response.
 */
std::vector<double>& Microambiente::vector_de_densidades( int n ){

	return (*p_vectores_densidad)[ n ];

}

/////////////////////GRADIENTES//////////////////////////////
/**
 * @brief Returns the gradient vector at the specified 3D voxel coordinates
 * 
 * @param i X-coordinate index
 * @param j Y-coordinate index
 * @param k Z-coordinate index
 * @return std::vector<gradiente>& Reference to the gradient vector at the specified location
 * 
 * Cancer Research Context:
 * Gradients of biochemical substrates (oxygen, nutrients, signaling molecules) are 
 * critical drivers of cancer cell behavior. This method provides access to these
 * gradients for modeling chemotaxis, cellular adaptation, and migration responses.
 * Steep gradients often exist around blood vessels and necrotic regions in tumors,
 * influencing invasion patterns and treatment response.
 */
std::vector<gradiente>& Microambiente::vector_de_gradientes(int i, int j, int k){

    int n = fg_indice_de_voxel(i,j,k);
	if( vector_gradiente_calculado[n] == false )
	{
		calcular_vector_de_gradiente( n );
	}

	return vectores_gradiente[n];

}

/**
 * @brief Returns the gradient vector at the specified voxel index
 * 
 * @param n Index of the voxel to return gradient vector for
 * @return std::vector<gradiente>& Reference to the gradient vector at the specified voxel
 * 
 * Cancer Research Context:
 * This method enables modeling of directional cellular responses to biochemical
 * gradients, which is essential for simulating cancer cell migration, invasion,
 * and adaptation to the tumor microenvironment. Gradient-driven behaviors are
 * key factors in metastatic potential and treatment resistance in cancer.
 */
std::vector<gradiente>& Microambiente::vector_de_gradientes(int n ){


	if( vector_gradiente_calculado[n] == false )
	{
		calcular_vector_de_gradiente( n );
	}


	return vectores_gradiente[n];
}



/**
 * @brief Returns the gradient vector of the voxel closest to the given position
 * 
 * @param posicion 3D spatial position to find nearest voxel
 * @return std::vector<gradiente>& Reference to the gradient vector at the nearest voxel
 * 
 * Cancer Research Context:
 * This method provides convenient access to biochemical gradients at any 3D position,
 * which is essential for modeling how cancer cells sense and respond to their local
 * microenvironment. This enables simulation of:
 * - Directional migration of cancer cells along nutrient or oxygen gradients
 * - Cellular responses to cytokine or growth factor gradients
 * - Invasion patterns driven by chemotactic responses
 * - Adaptation of cancer cells to local concentration changes
 */
std::vector<gradiente>& Microambiente::vector_de_gradiente_mas_cercano( Vector& posicion ){

	int n = indice_del_voxel_mas_cercano( posicion );
	if( vector_gradiente_calculado[n] == false )
	{
		calcular_vector_de_gradiente( n );
	}

	return vectores_gradiente[n];

}

/**
 * @brief Calculates gradient vectors for all substrates at all voxels in the domain
 * 
 * Cancer Research Context:
 * This method pre-computes all biochemical gradients throughout the entire tumor
 * microenvironment, which is essential for comprehensive analysis of:
 * - Global patterns of oxygen and nutrient gradients across the tumor
 * - Spatial heterogeneity in microenvironmental factors
 * - Potential migration paths for cancer cells responding to these gradients
 * - Distribution of hypoxic regions and other microenvironmental niches
 * 
 * Pre-computing all gradients enables efficient simulation of gradient-driven
 * behaviors across the entire tumor, though it is computationally intensive.
 * The method uses OpenMP parallelization to improve performance.
 */
void Microambiente::calcular_todos_los_vectores_de_gradientes( void ){

	static double dos_dx = mgrilla.dx;
	static double dos_dy = mgrilla.dy;
	static double dos_dz = mgrilla.dz;
	static bool constantes_de_gradientes_definidas = false;
	if( constantes_de_gradientes_definidas == false )
	{
		dos_dx *= 2.0;
		dos_dy *= 2.0;
		dos_dz *= 2.0;
		constantes_de_gradientes_definidas = true;
	}

	#pragma omp parallel for
	for( unsigned int k=0; k < mgrilla.coordenadas_z.size(); k++ )
	{
		for( unsigned int j=0; j < mgrilla.coordenadas_y.size(); j++ )
		{

			for( unsigned int q=0; q < numero_de_densidades(); q++ )
			{
				int i = 0;
				int n = fg_indice_de_voxel(i,j,k);

				vectores_gradiente[n][q][0] = (*p_vectores_densidad)[n+thomas_salto_en_i][q];
				vectores_gradiente[n][q][0] -= (*p_vectores_densidad)[n][q];
				vectores_gradiente[n][q][0] /= mgrilla.dx;

				vector_gradiente_calculado[n] = true;
			}
			for( unsigned int q=0; q < numero_de_densidades() ; q++ )
			{
				int i = mgrilla.coordenadas_x.size()-1;
				int n = fg_indice_de_voxel(i,j,k);

				vectores_gradiente[n][q][0] = (*p_vectores_densidad)[n][q];
				vectores_gradiente[n][q][0] -= (*p_vectores_densidad)[n-thomas_salto_en_i][q];
				vectores_gradiente[n][q][0] /= mgrilla.dx;

				vector_gradiente_calculado[n] = true;
			}

			for( unsigned int i=1; i < mgrilla.coordenadas_x.size()-1 ; i++ )
			{
				for( unsigned int q=0; q < numero_de_densidades() ; q++ )
				{
					int n = fg_indice_de_voxel(i,j,k);

					vectores_gradiente[n][q][0] = (*p_vectores_densidad)[n+thomas_salto_en_i][q];
					vectores_gradiente[n][q][0] -= (*p_vectores_densidad)[n-thomas_salto_en_i][q];
					vectores_gradiente[n][q][0] /= dos_dx;

					vector_gradiente_calculado[n] = true;
 				}
			}

		}
	}

	#pragma omp parallel for
	for( unsigned int k=0; k < mgrilla.coordenadas_z.size() ; k++ )
	{
		for( unsigned int i=0; i < mgrilla.coordenadas_x.size() ; i++ )
		{
			// endcaps
			for( unsigned int q=0; q < numero_de_densidades() ; q++ )
			{
				int j = 0;
				int n = fg_indice_de_voxel(i,j,k);

				vectores_gradiente[n][q][1] = (*p_vectores_densidad)[n+thomas_salto_en_j][q];
				vectores_gradiente[n][q][1] -= (*p_vectores_densidad)[n][q];
				vectores_gradiente[n][q][1] /= mgrilla.dy;

				vector_gradiente_calculado[n] = true;
			}
			for( unsigned int q=0; q < numero_de_densidades() ; q++ )
			{
				int j = mgrilla.coordenadas_y.size()-1;
				int n = fg_indice_de_voxel(i,j,k);

				vectores_gradiente[n][q][1] = (*p_vectores_densidad)[n][q];
				vectores_gradiente[n][q][1] -= (*p_vectores_densidad)[n-thomas_salto_en_j][q];
				vectores_gradiente[n][q][1] /= mgrilla.dy;

				vector_gradiente_calculado[n] = true;
			}

			for( unsigned int j=1; j < mgrilla.coordenadas_y.size()-1 ; j++ )
			{
				for( unsigned int q=0; q < numero_de_densidades() ; q++ )
				{
					int n = fg_indice_de_voxel(i,j,k);

					vectores_gradiente[n][q][1] = (*p_vectores_densidad)[n+thomas_salto_en_j][q];
					vectores_gradiente[n][q][1] -= (*p_vectores_densidad)[n-thomas_salto_en_j][q];
					vectores_gradiente[n][q][1] /= dos_dy;

					vector_gradiente_calculado[n] = true;
				}
			}

		}
	}


	#pragma omp parallel for
	for( unsigned int j=0; j < mgrilla.coordenadas_y.size() ; j++ )
	{
		for( unsigned int i=0; i < mgrilla.coordenadas_x.size() ; i++ )
		{

			for( unsigned int q=0; q < numero_de_densidades() ; q++ )
			{
				int k = 0;
				int n = fg_indice_de_voxel(i,j,k);

				vectores_gradiente[n][q][2] = (*p_vectores_densidad)[n+thomas_salto_en_k][q];
				vectores_gradiente[n][q][2] -= (*p_vectores_densidad)[n][q];
				vectores_gradiente[n][q][2] /= mgrilla.dz;

				vector_gradiente_calculado[n] = true;
			}
			for( unsigned int q=0; q < numero_de_densidades() ; q++ )
			{
				int k = mgrilla.coordenadas_z.size()-1;
				int n = fg_indice_de_voxel(i,j,k);

				vectores_gradiente[n][q][2] = (*p_vectores_densidad)[n][q];
				vectores_gradiente[n][q][2] -= (*p_vectores_densidad)[n-thomas_salto_en_k][q];
				vectores_gradiente[n][q][2] /= mgrilla.dz;

				vector_gradiente_calculado[n] = true;
			}

			for( unsigned int k=1; k < mgrilla.coordenadas_z.size()-1 ; k++ )
			{
				for( unsigned int q=0; q < numero_de_densidades() ; q++ )
				{
					int n = fg_indice_de_voxel(i,j,k);

					vectores_gradiente[n][q][2] = (*p_vectores_densidad)[n+thomas_salto_en_k][q];
					vectores_gradiente[n][q][2] -= (*p_vectores_densidad)[n-thomas_salto_en_k][q];
					vectores_gradiente[n][q][2] /= dos_dz;

					vector_gradiente_calculado[n] = true;
				}
			}

		}
	}

	return;

}

/**
 * @brief Calculates the gradient vectors for all substrates at the specified voxel
 * 
 * @param n Index of the voxel to calculate gradients for
 * 
 * Cancer Research Context:
 * This method computes the spatial gradients of biochemical substrates that drive
 * cancer cell behavior. These gradients are essential for modeling:
 * - Chemotaxis of cancer cells toward nutrient sources
 * - Cellular adaptation to hypoxia gradients around necrotic regions
 * - Directional migration of cancer cells during invasion and metastasis
 * - Heterogeneous response to treatment based on microenvironmental gradients
 * 
 * The method uses central differences to calculate gradients, which provides
 * accurate spatial derivatives for modeling how cancer cells sense and respond
 * to their local environment.
 */
void Microambiente::calcular_vector_de_gradiente( int n ){

	static double dos_dx = mgrilla.dx;
	static double dos_dy = mgrilla.dy;
	static double dos_dz = mgrilla.dz;
	static bool constantes_de_gradientes_definidas = false;
	std::vector<unsigned int> indices(3,0);

	if( constantes_de_gradientes_definidas == false )
	{
		dos_dx *= 2.0;
		dos_dy *= 2.0;
		dos_dz *= 2.0;
		constantes_de_gradientes_definidas = true;
	}


	indices = indices_cartesianos( n );


	if( indices[0] > 0 && indices[0] < mgrilla.coordenadas_x.size()-1 )
	{
		for( unsigned int q=0; q < numero_de_densidades() ; q++ )
		{
			vectores_gradiente[n][q][0] = (*p_vectores_densidad)[n+thomas_salto_en_i][q];
			vectores_gradiente[n][q][0] -= (*p_vectores_densidad)[n-thomas_salto_en_i][q];
			vectores_gradiente[n][q][0] /= dos_dx;

			vector_gradiente_calculado[n] = true;
		}
	}



	if( indices[1] > 0 && indices[1] < mgrilla.coordenadas_y.size()-1 )
	{
		for( unsigned int q=0; q < numero_de_densidades() ; q++ )
		{
			vectores_gradiente[n][q][1] = (*p_vectores_densidad)[n+thomas_salto_en_j][q];
			vectores_gradiente[n][q][1] -= (*p_vectores_densidad)[n-thomas_salto_en_j][q];
			vectores_gradiente[n][q][1] /= dos_dy;

			vector_gradiente_calculado[n] = true;
		}
	}



	if( indices[2] > 0 && indices[2] < mgrilla.coordenadas_z.size()-1 )
	{
		for( unsigned int q=0; q < numero_de_densidades() ; q++ )
		{
			vectores_gradiente[n][q][2] = (*p_vectores_densidad)[n+thomas_salto_en_k][q];
			vectores_gradiente[n][q][2] -= (*p_vectores_densidad)[n-thomas_salto_en_k][q];
			vectores_gradiente[n][q][2] /= dos_dz;

			vector_gradiente_calculado[n] = true;
		}
	}

	return;


}

/**
 * @brief Resets all gradient vectors, marking them as not calculated
 * 
 * Cancer Research Context:
 * This method is important when the microenvironment has changed significantly
 * (due to diffusion, cellular consumption, or secretion), requiring gradients
 * to be recalculated. This ensures that:
 * - Cancer cell responses to gradients are based on current conditions
 * - Simulations accurately reflect dynamic changes in the tumor microenvironment
 * - Computational resources are used efficiently by only calculating gradients when needed
 * 
 * This is particularly important after diffusion steps or when modeling rapid
 * changes in the tumor microenvironment such as acute hypoxia, drug delivery,
 * or vascular disruption.
 */
void Microambiente::resetear_todos_los_vectores_de_gradientes( void ){

	for( unsigned int k=0 ; k < mgrilla.voxeles.size() ; k++ )
	{
		for( unsigned int i=0 ; i < numero_de_densidades() ; i++ )
		{
			(vectores_gradiente[k])[i].resize( 3, 0.0 );
		}
	}
	vector_gradiente_calculado.assign( mgrilla.voxeles.size() , false );


}




/**
 * @brief Simulates diffusion and decay of all substrates for one time step
 * 
 * @param dt Time step size for the simulation
 * 
 * Cancer Research Context:
 * This method is central to modeling the tumor microenvironment dynamics by simulating
 * how biochemical factors (oxygen, nutrients, signaling molecules) diffuse through
 * the tissue and decay over time. These processes create the spatial gradients that:
 * - Drive hypoxic adaptation in cancer cells
 * - Create metabolically distinct regions within tumors
 * - Influence cancer cell migration and invasion patterns
 * - Affect drug delivery and efficacy in different tumor regions
 * 
 * The method uses a locally one-dimensional (LOD) scheme for 3D diffusion-reaction
 * equations, which is computationally efficient for modeling complex tumor
 * microenvironments.
 */
void Microambiente::simular_difusion_decaimiento( double dt ){

	solver_decaimiento_de_la_difusion__coeficientes_constantes_LOD_3D( dt );

	return;

}



/**
 * @brief Adds a Dirichlet boundary condition node at the specified voxel
 * 
 * @param indice_de_voxel Index of the voxel to set as a Dirichlet node
 * @param valor Vector of substrate values to maintain at this node
 * 
 * Cancer Research Context:
 * Dirichlet boundary conditions are essential for modeling fixed concentration sources
 * in the tumor microenvironment, such as:
 * - Blood vessels that maintain constant oxygen/nutrient levels
 * - Tissue boundaries with fixed metabolite concentrations
 * - Drug delivery sources with constant concentration
 * These boundary conditions create the concentration gradients that drive diffusion
 * and influence cancer cell behavior, metabolism, and treatment response.
 */
void Microambiente::agregar_nodo_de_dirichlet( int indice_de_voxel, std::vector<double>& valor ){

	mgrilla.voxeles[indice_de_voxel].es_dirichlet=true;

	vector_valores_de_dirichlet[indice_de_voxel] = valor; // .assign( mesh.voxels.size(), one );

	return;

}

/**
 * @brief Updates the values at an existing Dirichlet boundary condition node
 * 
 * @param indice_de_voxel Index of the Dirichlet node to update
 * @param nuevo_valor New vector of substrate values for this node
 * 
 * Cancer Research Context:
 * This method allows for dynamic changes to boundary conditions, which is crucial for modeling:
 * - Temporal changes in blood vessel perfusion affecting oxygen delivery
 * - Dynamic drug delivery profiles during treatment
 * - Changes in metabolite availability at tissue boundaries
 * These dynamic changes can significantly impact tumor growth patterns, hypoxic regions,
 * and treatment efficacy in cancer simulations.
 */
void Microambiente::actualizar_nodo_de_dirichlet( int indice_de_voxel , std::vector<double>& nuevo_valor ){

	mgrilla.voxeles[indice_de_voxel].es_dirichlet = true;
	vector_valores_de_dirichlet[indice_de_voxel] = nuevo_valor;

	return;
}

/**
 * @brief Updates a single substrate value at an existing Dirichlet boundary condition node
 * 
 * @param indice_de_voxel Index of the Dirichlet node to update
 * @param indice_de_sustrato Index of the substrate to update
 * @param nuevo_valor New value for the specified substrate
 * 
 * Cancer Research Context:
 * This method enables selective modification of individual substrate concentrations
 * at boundary points, which is important for modeling:
 * - Selective changes in specific signaling molecules
 * - Targeted drug delivery affecting only certain substrates
 * - Differential regulation of nutrients or metabolites
 * This granular control is essential for studying how specific microenvironmental
 * factors influence cancer cell behavior and treatment response.
 */
void Microambiente::actualizar_nodo_de_dirichlet( int indice_de_voxel , int indice_de_sustrato, double nuevo_valor){

	mgrilla.voxeles[indice_de_voxel].es_dirichlet = true;
	vector_valores_de_dirichlet[indice_de_voxel][indice_de_sustrato] = nuevo_valor;

	dirichlet_vectorES_activacion[indice_de_voxel][indice_de_sustrato] = true;

	return;
}

/**
 * @brief Applies all Dirichlet boundary conditions to enforce fixed concentrations
 * 
 * Cancer Research Context:
 * This method enforces fixed concentration boundary conditions throughout the
 * tumor microenvironment, which is essential for modeling:
 * - Blood vessels that maintain constant oxygen/nutrient supply
 * - Interfaces with surrounding normal tissue
 * - Experimental conditions with controlled nutrient or drug delivery
 * 
 * By applying these boundary conditions, the method creates realistic concentration
 * gradients within the tumor that drive critical biological processes including
 * hypoxic adaptation, metabolic zonation, and directional cell migration.
 * This enforces physiologically relevant constraints on the diffusion process.
 */
void Microambiente::aplicar_condiciones_de_dirichlet( void ){

        #pragma omp parallel for
        for( unsigned int i=0 ; i < mgrilla.voxeles.size() ;i++ )
        {
                if( mgrilla.voxeles[i].es_dirichlet == true )
                {
                        for( unsigned int j=0; j < vector_valores_de_dirichlet[i].size(); j++ )
                        {
                                if( dirichlet_vectorES_activacion[i][j] == true )
                                {
                                        vector_de_densidades(i)[j] = vector_valores_de_dirichlet[i][j]; //38;
                                }
                        }

                }
        }
        return;
}

/**
 * @brief Sets the activation state of a substrate at all Dirichlet nodes
 * 
 * @param indice_de_sustrato Index of the substrate to set activation for
 * @param nuevo_valor New activation state (true = active, false = inactive)
 * 
 * Cancer Research Context:
 * This method enables global control over whether specific substrates are
 * maintained at fixed concentrations at boundary points. This is important for modeling:
 * - Selective permeability of tumor vasculature to different molecules
 * - Changes in availability of specific nutrients or drugs
 * - Experimental conditions where certain factors are controlled
 * This control is essential for studying how specific microenvironmental factors
 * influence tumor growth, invasion, and response to therapy.
 */
void Microambiente::set_activacion_de_sustrato_de_dirichlet( int indice_de_sustrato , bool nuevo_valor ){

	vector_activacion_dirichlet[indice_de_sustrato] = nuevo_valor;

	for( unsigned int n = 0 ; n < mgrilla.voxeles.size() ; n++ )
	{ dirichlet_vectorES_activacion[n][indice_de_sustrato] = nuevo_valor; }

	return;
}

/**
 * @brief Sets the activation state of a substrate at a specific Dirichlet node
 * 
 * @param indice_de_sustrato Index of the substrate to set activation for
 * @param indice Index of the Dirichlet node
 * @param nuevo_valor New activation state (true = active, false = inactive)
 * 
 * Cancer Research Context:
 * This method provides localized control over boundary conditions for specific
 * substrates at specific locations. This is crucial for modeling:
 * - Heterogeneous vascular permeability within tumors
 * - Localized drug delivery or metabolite sources
 * - Spatial variations in tissue-tumor interfaces
 * This spatial control allows for more realistic modeling of the heterogeneous
 * tumor microenvironment that influences cancer progression and treatment response.
 */
void Microambiente::set_activacion_de_sustrato_de_dirichlet( int indice_de_sustrato , int indice, bool nuevo_valor ){

	dirichlet_vectorES_activacion[indice][indice_de_sustrato] = nuevo_valor;
	return;

}

/**
 * @brief Sets the activation state of all substrates at a specific Dirichlet node
 * 
 * @param indice Index of the Dirichlet node
 * @param nuevo_valor Vector of activation states for all substrates
 * 
 * Cancer Research Context:
 * This method allows for comprehensive control of all substrates at a specific
 * boundary location, which is important for modeling:
 * - Blood vessels with specific permeability profiles
 * - Tissue interfaces with selective transport properties
 * - Experimental boundary conditions with controlled substrate availability
 * This comprehensive control enables realistic simulation of complex boundary
 * conditions that affect multiple aspects of the tumor microenvironment simultaneously.
 */
void Microambiente::set_activacion_de_sustrato_de_dirichlet( int indice, std::vector<bool>& nuevo_valor ){

	dirichlet_vectorES_activacion[indice] = nuevo_valor;
	return;

}

/**
 * @brief Gets the activation state of a substrate at a specific Dirichlet node
 * 
 * @param indice_de_sustrato Index of the substrate to check
 * @param indice Index of the Dirichlet node
 * @return bool Current activation state (true = active, false = inactive)
 * 
 * Cancer Research Context:
 * This method enables querying the current state of boundary conditions, which is
 * important for:
 * - Adaptive modeling based on current microenvironmental state
 * - Analysis of boundary effects on tumor growth and behavior
 * - Verification of simulation conditions during complex scenarios
 * This information helps track how boundary conditions influence cancer cell
 * behavior and microenvironmental evolution during simulation.
 */
bool Microambiente::get_activacion_de_sustrato_de_dirichlet( int indice_de_sustrato, int indice ){

	return vector_activacion_dirichlet[indice_de_sustrato];

}

/**
 * @brief Gets a reference to the Dirichlet node status of the specified voxel
 * 
 * @param indice_de_voxel Index of the voxel to check
 * @return bool& Reference to the Dirichlet status (true = is a Dirichlet node, false = not a Dirichlet node)
 * 
 * Cancer Research Context:
 * This method provides access to check or modify the Dirichlet boundary status of specific
 * locations in the tumor microenvironment. This is important for:
 * - Identifying fixed concentration sources like blood vessels in the tumor
 * - Determining which regions maintain constant biochemical conditions
 * - Analyzing how boundary conditions affect local cancer cell behavior
 * - Tracking the spatial distribution of nutrient/oxygen sources in the simulation
 */
bool& Microambiente::es_nodo_de_dirichlet( int indice_de_voxel ){

        return mgrilla.voxeles[indice_de_voxel].es_dirichlet;
}

/**
 * @brief Implements a locally one-dimensional (LOD) 3D diffusion-reaction solver
 * 
 * @param dt Time step size for the simulation
 * 
 * Cancer Research Context:
 * This method is the core numerical implementation for simulating how biochemical
 * factors diffuse and decay within the tumor microenvironment. It is critical for modeling:
 * - Oxygen gradients that drive cancer cell hypoxic adaptation
 * - Nutrient distribution that affects cancer cell metabolism and proliferation
 * - Drug penetration through tumor tissue for treatment efficacy studies
 * - Signaling molecule gradients that influence cancer cell phenotypes
 * 
 * The method implements a computationally efficient LOD (Alternating Direction Implicit)
 * approach that solves the diffusion equation sequentially along each spatial dimension.
 * This enables realistic 3D modeling of complex biochemical gradients that are known to
 * significantly influence tumor progression, invasion, and treatment response.
 */
void Microambiente::solver_decaimiento_de_la_difusion__coeficientes_constantes_LOD_3D( double dt ){



	if( !setup_del_solver_de_difusion_hecho )
	{


		thomas_denomx.resize( mgrilla.coordenadas_x.size() , cero );
		thomas_cx.resize( mgrilla.coordenadas_x.size() , cero );

		thomas_denomy.resize( mgrilla.coordenadas_y.size() , cero );
		thomas_cy.resize( mgrilla.coordenadas_y.size() , cero );

		thomas_denomz.resize( mgrilla.coordenadas_z.size() , cero );
		thomas_cz.resize( mgrilla.coordenadas_z.size() , cero );

		thomas_salto_en_i = 1;
		thomas_salto_en_j = mgrilla.coordenadas_x.size();
		thomas_salto_en_k = thomas_salto_en_j * mgrilla.coordenadas_y.size();

		thomas_constante1 =  coeficientes_de_difusion;
		thomas_constante1a = cero;
		thomas_constante2 =  tasas_de_decaimiento;
		thomas_constante3 = uno;
		thomas_constante3a = uno;

		thomas_constante1 *= dt;
		thomas_constante1 /= mgrilla.dx;
		thomas_constante1 /= mgrilla.dx;

		thomas_constante1a = thomas_constante1;
		thomas_constante1a *= -1.0;

		thomas_constante2 *= dt;
		thomas_constante2 /= 3.0;

		thomas_constante3 += thomas_constante1;
		thomas_constante3 += thomas_constante1;
		thomas_constante3 += thomas_constante2;

		thomas_constante3a += thomas_constante1;
		thomas_constante3a += thomas_constante2;

		// Thomas solver coefficients

		thomas_cx.assign( mgrilla.coordenadas_x.size() , thomas_constante1a );
		thomas_denomx.assign( mgrilla.coordenadas_x.size()  , thomas_constante3 );
		thomas_denomx[0] = thomas_constante3a;
		thomas_denomx[ mgrilla.coordenadas_x.size()-1 ] = thomas_constante3a;
		if( mgrilla.coordenadas_x.size() == 1 )
		{ thomas_denomx[0] = uno; thomas_denomx[0] += thomas_constante2; }

		thomas_cx[0] /= thomas_denomx[0];
		for( unsigned int i=1 ; i <= mgrilla.coordenadas_x.size()-1 ; i++ )
		{
			axpy( &thomas_denomx[i] , thomas_constante1 , thomas_cx[i-1] );
			thomas_cx[i] /= thomas_denomx[i];
		}

		thomas_cy.assign( mgrilla.coordenadas_y.size() , thomas_constante1a );
		thomas_denomy.assign( mgrilla.coordenadas_y.size()  , thomas_constante3 );
		thomas_denomy[0] = thomas_constante3a;
		thomas_denomy[ mgrilla.coordenadas_y.size()-1 ] = thomas_constante3a;
		if( mgrilla.coordenadas_y.size() == 1 )
		{thomas_denomy[0] = uno; thomas_denomy[0] += thomas_constante2; }

		thomas_cy[0] /= thomas_denomy[0];
		for( unsigned int i=1 ; i <= mgrilla.coordenadas_y.size()-1 ; i++ )
		{
			axpy( &thomas_denomy[i] , thomas_constante1 , thomas_cy[i-1] );
			thomas_cy[i] /= thomas_denomy[i];
		}

		thomas_cz.assign( mgrilla.coordenadas_z.size() , thomas_constante1a );
		thomas_denomz.assign( mgrilla.coordenadas_z.size()  , thomas_constante3 );
		thomas_denomz[0] = thomas_constante3a;
		thomas_denomz[ mgrilla.coordenadas_z.size()-1 ] = thomas_constante3a;
		if( mgrilla.coordenadas_z.size() == 1 )
		{ thomas_denomz[0] = uno; thomas_denomz[0] += thomas_constante2; }

		thomas_cz[0] /= thomas_denomz[0];
		for( unsigned int i=1 ; i <= mgrilla.coordenadas_z.size()-1 ; i++ )
		{
			axpy( &thomas_denomz[i] , thomas_constante1 , thomas_cz[i-1] );
			thomas_cz[i] /= thomas_denomz[i];
		}

		setup_del_solver_de_difusion_hecho = true;

	}





	aplicar_condiciones_de_dirichlet();

	#pragma omp parallel for
	for( unsigned int k=0; k < mgrilla.coordenadas_z.size() ; k++ )
	{
		for( unsigned int j=0; j < mgrilla.coordenadas_y.size() ; j++ )
		{



			int n = fg_indice_de_voxel(0,j,k);
			(*p_vectores_densidad)[n] /= thomas_denomx[0];

			for( unsigned int i=1; i < mgrilla.coordenadas_x.size() ; i++ )
			{
				n = fg_indice_de_voxel(i,j,k);

                for( unsigned int m=0; m < (*p_vectores_densidad)[n].size() ; m++ ){
                    (*p_vectores_densidad)[n][m] += thomas_constante1[m] * (*p_vectores_densidad)[n-thomas_salto_en_i][m];
                }
				(*p_vectores_densidad)[n] /= thomas_denomx[i];
			}

			for( int i = mgrilla.coordenadas_x.size()-2 ; i >= 0 ; i-- )
			{
				n = fg_indice_de_voxel(i,j,k);

                for( unsigned int m=0; m < (*p_vectores_densidad)[n].size() ; m++ ){
                    (*p_vectores_densidad)[n][m] -= thomas_cx[i][m] * (*p_vectores_densidad)[n+thomas_salto_en_i][m];
                }                
			}

		}
	}



	aplicar_condiciones_de_dirichlet();

	#pragma omp parallel for
	for( unsigned int k=0; k < mgrilla.coordenadas_z.size() ; k++ )
	{
		for( unsigned int i=0; i < mgrilla.coordenadas_x.size() ; i++ )
		{




	int n = fg_indice_de_voxel(i,0,k);
	(*p_vectores_densidad)[n] /= thomas_denomy[0];

	for( unsigned int j=1; j < mgrilla.coordenadas_y.size() ; j++ )
	{
		n = fg_indice_de_voxel(i,j,k);

        for( unsigned int m=0; m < (*p_vectores_densidad)[n].size() ; m++ ){
            (*p_vectores_densidad)[n][m] += thomas_constante1[m] * (*p_vectores_densidad)[n-thomas_salto_en_j][m];
        }        
		(*p_vectores_densidad)[n] /= thomas_denomy[j];
	}




	for( int j = mgrilla.coordenadas_y.size()-2 ; j >= 0 ; j-- )
	{
		n = fg_indice_de_voxel(i,j,k);

        for( unsigned int m=0; m < (*p_vectores_densidad)[n].size() ; m++ ){
            (*p_vectores_densidad)[n][m] -= thomas_cy[j][m] * (*p_vectores_densidad)[n+thomas_salto_en_j][m];
        }        
	}

  }
 }



	aplicar_condiciones_de_dirichlet();

 #pragma omp parallel for
 for( unsigned int j=0; j < mgrilla.coordenadas_y.size() ; j++ )
 {

  for( unsigned int i=0; i < mgrilla.coordenadas_x.size() ; i++ )
  {




	int n = fg_indice_de_voxel(i,j,0);
	(*p_vectores_densidad)[n] /= thomas_denomz[0];


	for( unsigned int k=1; k < mgrilla.coordenadas_z.size() ; k++ )
	{
		n = fg_indice_de_voxel(i,j,k);

        for( unsigned int m=0; m < (*p_vectores_densidad)[n].size() ; m++ ){
            (*p_vectores_densidad)[n][m] += thomas_constante1[m] * (*p_vectores_densidad)[n-thomas_salto_en_k][m];
        }           
		(*p_vectores_densidad)[n] /= thomas_denomz[k];
	}




	for( int k = mgrilla.coordenadas_z.size()-2 ; k >= 0 ; k-- )
	{
		n = fg_indice_de_voxel(i,j,k);

        for( unsigned int m=0; m < (*p_vectores_densidad)[n].size() ; m++ ){
            (*p_vectores_densidad)[n][m] -= thomas_cz[k][m] * (*p_vectores_densidad)[n+thomas_salto_en_k][m];
        }        

	}
  }
 }

	aplicar_condiciones_de_dirichlet();




	return;
}

/**
 * @brief Displays detailed information about the microenvironment configuration
 * 
 * @param os Output stream where information will be written
 * 
 * Cancer Research Context:
 * Provides comprehensive information about the microenvironment setup,
 * including diffusion lengths, which are critical for understanding
 * how far oxygen and nutrients can penetrate into tumor tissue. This
 * helps researchers interpret simulation results in the context of
 * known tumor biology, where oxygen diffusion typically limits to
 * 100-200 microns from blood vessels.
 */
void Microambiente::mostrar_informacion( std::ostream& os )
{
	os << std::endl << "Resumen del Microambiente: " << nombre << ": " << std::endl;
	mgrilla.mostrar_informacion_cartesiano( os );
	os << "Densidades: (" << numero_de_densidades() << " en total)" << std::endl;
	for( unsigned int i = 0 ; i < densidades_nombres.size() ; i++ )
	{
		os << "   " << densidades_nombres[i] << ":" << std::endl
		<< "     unidades: " << densidades_unidades[i] << std::endl
		<< "     coeficiente de difusion: " << coeficientes_de_difusion[i]
			<< " " << unidades_espaciales << "^2 / " << unidades_temporales << std::endl
		<< "     tasa de decaimiento: " << tasas_de_decaimiento[i]
			<< " " << unidades_temporales << "^-1" << std::endl
		<< "     longitud de escala de la difusion: " << sqrt( coeficientes_de_difusion[i] / ( 1e-12 + tasas_de_decaimiento[i] ) )
			<< " " << unidades_espaciales << std::endl
		<< "     condicion inicial: " << pg->vector_condiciones_iniciales[i]
			<< " " << densidades_unidades[i] << std::endl
		<< "     condiciones de borde: " << pg->vector_condicion_de_dirichlet[i]
			<< " " << densidades_unidades[i] << " (activo: ";
		if( vector_activacion_dirichlet[i] == true )
		{ os << "si"; }
		else
		{ os << "no"; }
		os << ")" << std::endl;
	}
	os << std::endl;

	return;
}



//Microambiente microambiente;
//Microambiente_Parametros p;

/**
 * @brief Initializes the microenvironment with parameters from global settings
 * 
 * Sets up the microenvironment based on global parameter settings, including
 * substrate types, domain size, boundary conditions, and initial concentrations.
 * By default, configures oxygen as the first substrate if specified.
 * 
 * Cancer Research Context:
 * Proper initialization of the microenvironment is essential for cancer simulation,
 * as it establishes the biochemical landscape that influences tumor growth.
 * Oxygen is typically set as the primary substrate due to its critical role in
 * tumor metabolism, with values calibrated to physiological levels (38 mmHg).
 * The option to include immune factors enables modeling of cancer-immune interactions.
 */
void Microambiente::inicializar_microambiente(){


	nombre = pg->m_nombre;


	if( pg->usar_oxigeno_como_primer_sustrato == true )
	{
		set_densidad(0, "oxigeno" , "mmHg");
		coeficientes_de_difusion[0] = 1e5;
		tasas_de_decaimiento[0] = 0.1;
	}

	if( pg->activar_respuesta_inmune == true )
	{
		agregar_densidad("immunostimulatory factor" , "dimensionless", 1000, .016);
		vector_activacion_dirichlet[1] = false;
		pg->vector_condiciones_iniciales[1] = 0.0;
	}


	redimensionar_espacio( pg->rango_en_X[0], pg->rango_en_X[1] ,
		pg->rango_en_Y[0], pg->rango_en_Y[1],
		pg->rango_en_Z[0], pg->rango_en_Z[1],
		pg->m_dx,pg->m_dy,pg->m_dz );


	unidades_espaciales = pg->unidades_espaciales;
	unidades_temporales = pg->unidades_temporales;
	mgrilla.unidades = pg->unidades_espaciales;




	if( pg->vector_condiciones_iniciales.size() !=
		numero_de_densidades() )
	{
		pg->vector_condiciones_iniciales = pg->vector_condicion_de_dirichlet;
	}


	for( unsigned int n=0; n < numero_de_voxeles() ; n++ )
	{ vector_de_densidades(n) = pg->vector_condiciones_iniciales; }




	pg->dirichlet_xmin_valores = pg->vector_condicion_de_dirichlet;
	pg->dirichlet_xmax_valores = pg->vector_condicion_de_dirichlet;
	pg->dirichlet_ymin_valores = pg->vector_condicion_de_dirichlet;
	pg->dirichlet_ymax_valores = pg->vector_condicion_de_dirichlet;
	pg->dirichlet_zmin_valores = pg->vector_condicion_de_dirichlet;
	pg->dirichlet_zmax_valores = pg->vector_condicion_de_dirichlet;



	bool xmin = false;
	bool xmax = false;
	bool ymin = false;
	bool ymax = false;
	bool zmin = false;
	bool zmax = false;

	if( pg->condiciones_de_Dirichlet_externas == true )
	{
		for( unsigned int n=0 ; n < numero_de_densidades() ; n++ )
		{
			if( pg->dirichlet_todo[n] ||
				pg->dirichlet_xmin[n] )
				{ xmin = true; }

			if( pg->dirichlet_todo[n] ||
				pg->dirichlet_xmax[n] )
				{ xmax = true; }

			if( pg->dirichlet_todo[n] ||
				pg->dirichlet_ymin[n] )
				{ ymin = true; }

			if( pg->dirichlet_todo[n] ||
				pg->dirichlet_ymax[n] )
				{ ymax = true; }

			if( pg->dirichlet_todo[n] ||
				pg->dirichlet_zmin[n] )
				{ zmin = true; }

			if( pg->dirichlet_todo[n] ||
				pg->dirichlet_zmax[n] )
				{ zmax = true; }
		}



	}




	if( pg->condiciones_de_Dirichlet_externas == true )
	{

		if( xmin == true )
		{
			for( unsigned int k=0 ; k < mgrilla.coordenadas_z.size() ; k++ )
			{
				int I = 0;

				for( unsigned int j=0 ; j < mgrilla.coordenadas_y.size() ; j++ )
				{

					agregar_nodo_de_dirichlet( fg_indice_de_voxel(I,j,k) , pg->dirichlet_xmin_valores );


					set_activacion_de_sustrato_de_dirichlet( fg_indice_de_voxel(I,j,k) ,
					pg->dirichlet_xmin );

				}
			}
		}


		if( xmax == true )
		{
			for( unsigned int k=0 ; k < mgrilla.coordenadas_z.size() ; k++ )
			{
				int I = mgrilla.coordenadas_x.size()-1;;

				for( unsigned int j=0 ; j < mgrilla.coordenadas_y.size() ; j++ )
				{

					agregar_nodo_de_dirichlet( fg_indice_de_voxel(I,j,k) , pg->dirichlet_xmax_valores );


					set_activacion_de_sustrato_de_dirichlet( fg_indice_de_voxel(I,j,k) ,
					pg->dirichlet_xmax );
				}
			}
		}


		if( ymin == true )
		{
			for( unsigned int k=0 ; k < mgrilla.coordenadas_z.size() ; k++ )
			{
				int J = 0;

				for( unsigned int i=0 ; i < mgrilla.coordenadas_x.size() ; i++ )
				{

					agregar_nodo_de_dirichlet( fg_indice_de_voxel(i,J,k) , pg->dirichlet_ymin_valores );


					set_activacion_de_sustrato_de_dirichlet( fg_indice_de_voxel(i,J,k) ,
					pg->dirichlet_ymin );
				}
			}
		}


		if( ymax == true )
		{
			for( unsigned int k=0 ; k < mgrilla.coordenadas_z.size() ; k++ )
			{
				int J = mgrilla.coordenadas_y.size()-1;;

				for( unsigned int i=0 ; i < mgrilla.coordenadas_x.size() ; i++ )
				{

					agregar_nodo_de_dirichlet( fg_indice_de_voxel(i,J,k) , pg->dirichlet_ymax_valores );


					set_activacion_de_sustrato_de_dirichlet(fg_indice_de_voxel(i,J,k) ,
					pg->dirichlet_ymax );
				}
			}
		}


			if( zmin == true )
			{
				for( unsigned int j=0 ; j < mgrilla.coordenadas_y.size() ; j++ )
				{
					int K = 0;

					for( unsigned int i=0 ; i < mgrilla.coordenadas_x.size() ; i++ )
					{

						agregar_nodo_de_dirichlet( fg_indice_de_voxel(i,j,K) , pg->dirichlet_zmin_valores );


						set_activacion_de_sustrato_de_dirichlet( fg_indice_de_voxel(i,j,K),
						pg->dirichlet_zmin );
					}
				}
			}


			if( zmax == true )
			{
				for( unsigned int j=0 ; j < mgrilla.coordenadas_y.size() ; j++ )
				{
					int K = mgrilla.coordenadas_z.size()-1;;

					for( unsigned int i=0 ; i < mgrilla.coordenadas_x.size() ; i++ )
					{

						agregar_nodo_de_dirichlet( fg_indice_de_voxel(i,j,K) , pg->dirichlet_zmax_valores );


						set_activacion_de_sustrato_de_dirichlet( fg_indice_de_voxel(i,j,K) ,
						pg->dirichlet_zmax );
					}
				}
			}
		}



	for( unsigned int i=0 ; i < pg->vector_activacion_dirichlet.size(); i++ )
	{
		set_activacion_de_sustrato_de_dirichlet( i , pg->vector_activacion_dirichlet[i] );
	}

	mostrar_informacion( std::cout );
	return;
}



/**
 * @brief Creates a blood vessel within the specified 3D region
 * 
 * @param xmin Minimum x-coordinate of the vessel region
 * @param ymin Minimum y-coordinate of the vessel region
 * @param zmin Minimum z-coordinate of the vessel region
 * @param xmax Maximum x-coordinate of the vessel region
 * @param ymax Maximum y-coordinate of the vessel region
 * @param zmax Maximum z-coordinate of the vessel region
 * 
 * Cancer Research Context:
 * Blood vessels are critical components of the tumor microenvironment that:
 * - Supply oxygen and nutrients to cancer cells
 * - Create oxygen and nutrient gradients that drive tumor heterogeneity
 * - Serve as routes for metastatic dissemination
 * - Provide delivery pathways for therapeutic agents
 * 
 * This method creates a blood vessel by setting Dirichlet boundary conditions
 * within the specified region, maintaining fixed substrate concentrations that
 * model the continuous supply from the vasculature. This is essential for
 * realistic modeling of tumor growth, hypoxia development, and treatment response.
 */
void Microambiente::crear_vaso_sanguineo( int xmin, int ymin, int zmin, int xmax, int ymax, int zmax ){


    int x1 = xmax, y1 = ymax, z1 = zmax, x0 = xmin, y0 = ymin, z0 = zmin;

    int dx = std::max(abs(x1 - x0), 1); //abs(x1 - x0);
    int dy = std::max(abs(y1 - y0), 1); //abs(y1 - y0);
    int dz = std::max(abs(z1 - z0), 1); //abs(z1 - z0);


//    int stepX;
//    if (x0 < x1){
//        stepX = 1;
//    }else if (x0 > x1){
//        stepX = -1;
//    }else{
//        stepX = 0;
//    }
    int stepX = x0 < x1 ? 1 : -1;

//    int stepY;
//    if (y0 < y1){
//        stepY = 1;
//    }else if (y0 > y1){
//        stepY = -1;
//    }else{
//        stepY = 0;
//    }
    int stepY = y0 < y1 ? 1 : -1;

//    int stepZ;
//    if (z0 < z1){
//        stepZ = 1;
//    }else if (z0 > z1){
//        stepZ = -1;
//    }else{
//        stepZ = 0;
//    }
    int stepZ = z0 < z1 ? 1 : -1;



    double hipotenusa = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));


    double tMaxX = hipotenusa*0.5 / dx;
    double tMaxY = hipotenusa*0.5 / dy;
    double tMaxZ = hipotenusa*0.5 / dz;


    double tDeltaX = hipotenusa / dx;
    double tDeltaY = hipotenusa / dy;
    double tDeltaZ = hipotenusa / dz;


    int voxel_actual = -99;

    while (x0 != x1 || y0 != y1 || z0 != z1){
        if(tMaxX < tMaxY){
            if(tMaxX < tMaxZ){
                x0 = x0 + stepX;
                tMaxX = tMaxX + tDeltaX;
            }else if(tMaxX < tMaxZ){
                z0 = z0 + stepZ;
                tMaxZ = tMaxZ + tDeltaZ;
            }else{
                x0 = x0 + stepX;
                tMaxX = tMaxX + tDeltaX;
                z0 = z0 + stepZ;
                tMaxZ = tMaxZ + tDeltaZ;
            }
        }else if(tMaxX > tMaxY){
            if(tMaxY < tMaxZ){
                y0 = y0 + stepY;
                tMaxY = tMaxY + tDeltaY;
            }else if(tMaxY > tMaxZ){
                z0 = z0 + stepZ;
                tMaxZ = tMaxZ + tDeltaZ;
            }else{
                y0 = y0 + stepY;
                tMaxY = tMaxY + tDeltaY;
                z0 = z0 + stepZ;
                tMaxZ = tMaxZ + tDeltaZ;
            }
        }else{
            if(tMaxY < tMaxZ){
                y0 = y0 + stepY;
                tMaxY = tMaxY + tDeltaY;
                x0 = x0 + stepX;
                tMaxX = tMaxX + tDeltaX;
            }else if(tMaxY > tMaxZ){
                z0 = z0 + stepZ;
                tMaxZ = tMaxZ + tDeltaZ;
            }else{
                x0 = x0 + stepX;
                tMaxX = tMaxX + tDeltaX;
                y0 = y0 + stepY;
                tMaxY = tMaxY + tDeltaY;
                z0 = z0 + stepZ;
                tMaxZ = tMaxZ + tDeltaZ;
            }
        }

        Vector posicion;
        posicion.x=x0;
        posicion.y=y0;
        posicion.z=z0;

        int voxel_siguiente = indice_del_voxel_mas_cercano(posicion);



        if(voxel_actual != voxel_siguiente){

            voxeles_del_vaso_sanguineo.push_back(indice_del_voxel_mas_cercano(posicion));




            agregar_nodo_de_dirichlet(indice_del_voxel_mas_cercano(posicion), pg->vector_condicion_de_dirichlet);



            set_activacion_de_sustrato_de_dirichlet(indice_del_voxel_mas_cercano(posicion), pg->dirichlet_vs);

            voxel_actual = voxel_siguiente;
        }
    }

//    for(long unsigned int j=0; j<voxeles_del_vaso_sanguineo.size(); j++){
//
//        std::cout << "Voxel: " << voxeles_del_vaso_sanguineo[j] << "\n";
//
//    };
//	std::cout << "pulse enter para continuar \n";
//
//    std::cin.get();
    //voxeles_del_vaso_sanguineo.clear();
////////////////////////////////////////////////////////////////
//	for( int k=zmin ; k < zmax ; k++ )
//	{
//		 //set Dirichlet conditions along the xmin outer edges
//		for( int j=ymin ; j < ymax ; j++ )
//		{
//            for( int i=xmin ; i < xmax ; i++ )
//            {
//
//                 //set the value
//                m->agregar_nodo_de_dirichlet( m->indice_de_voxel(i,j,k) , pg->vector_condicion_de_dirichlet );
//
//                 s//et the activation
//                m->set_activacion_de_sustrato_de_dirichlet( m->indice_de_voxel(i,j,k), true );
//
//            }
//        }
//    }

    return;
}
