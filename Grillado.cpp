/**
 * @file Grillado.cpp
 *
 * @author Luciana Melina Luque
 *
 * @details Implementation of spatial grid system for cancer simulation microenvironment.
 * 
 * This file implements the Cartesian grid structure that provides spatial organization for the
 * cancer simulation microenvironment. The grid system enables precise modeling of spatial
 * heterogeneity within tumors, which is critical for understanding cancer progression.
 * 
 * The implementation includes:
 * - General grid functionality (Grillado_General)
 * - Cartesian grid specialization (Grillado_Cartesiano)
 * - Methods for spatial indexing and voxel connectivity
 * - Grid resizing and domain configuration
 * 
 * In cancer research, spatial discretization is essential for accurately modeling:
 * - Oxygen and nutrient gradients that drive tumor growth patterns
 * - Heterogeneous microenvironments that influence cell phenotypes
 * - Spatial organization of cells with different proliferation capacities
 * - Invasion patterns and tumor morphology
 * 
 * Inbound Dependencies: Vector.h, Voxel.h, Parametros_globales.h, Macros.h
 *
 * Outbound Dependencies: Microambiente.h, Contenedor_de_Celulas.h
 *
 * Usage: Used by the microenvironment to implement the finite volume method for diffusion and
 *        by the cell container for spatial organization of cells.
 *
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#include "Grillado.h"

const int indice_de_x_min=0;
const int indice_de_y_min=1;
const int indice_de_z_min=2;
const int indice_de_x_max=3;
const int indice_de_y_max=4;
const int indice_de_z_max=5;

/**
 * @brief Constructor for the Grillado_General class.
 * 
 * Initializes a default grid with a small test domain (-0.5 to 0.5 in each dimension).
 * Sets up the basic spatial structure that will contain voxels for the simulation,
 * which is critical for modeling the tumor microenvironment in cancer research.
 * 
 * In cancer simulations, this grid provides the spatial framework for:
 * - Positioning cells within the tumor and surrounding tissue
 * - Representing biochemical gradients within and around tumors
 * - Modeling cell-cell and cell-microenvironment interactions
 * - Supporting finite volume methods for diffusion calculations
 * 
 * @return Newly constructed Grillado_General object with default parameters
 */
Grillado_General::Grillado_General()
{
	// x1, x2, y1, y2, z1, z2
	caja.assign(6,0.0);
	caja[indice_de_x_min] = -0.5;
	caja[indice_de_y_min] = -0.5;
	caja[indice_de_z_min] = -0.5;
	caja[indice_de_x_max] = 0.5;
	caja[indice_de_y_max] = 0.5;
	caja[indice_de_z_max] = 0.5;

	voxeles.resize( 1 );

	indices_de_voxeles_conectados.resize( 1 );
	indices_de_voxeles_conectados[0].clear();



	unidades = "micrometros";
}

/**
 * @brief Stream output operator for Grillado_General class.
 * 
 * Enables formatted output of a grid to ostreams like cout, file streams, or string streams.
 * Displays voxel properties and connectivity information, which is useful for debugging
 * and visualizing the spatial organization within the cancer microenvironment.
 * 
 * In cancer research, visualizing grid properties helps:
 * - Verify grid resolution for accurate tumor microenvironment modeling
 * - Understand voxel connectivity patterns that affect tumor growth dynamics
 * - Debug spatial components of the simulation
 * - Generate connectivity information for analysis of spatial tumor organization
 * 
 * @param os The output stream where the grid information will be written
 * @param grillado The grid object to display
 * @return Reference to the output stream for chaining
 */
std::ostream& operator<<(std::ostream& os, const Grillado_General& grillado)
{
	std::boolalpha( os );
	static std::string tabbing = "\t";
	static std::string tabbing2 = "\t\t";
	static std::string tabbing3 = "\t\t\t";

	for( unsigned int i=0; i < grillado.voxeles.size() ; i++ )
	{ os << grillado.voxeles[i] << std::endl; }
	os 	<< tabbing2 << "</voxeles>" << std::endl;

	os	<< tabbing2 << "<conexiones>" << std::endl;
	for( unsigned int i=0 ; i < grillado.indices_de_voxeles_conectados.size() ; i++ )
	{
		os << tabbing3 << "<ID de los voxeles conectados=\"" << i << "\">" << std::endl;
		for( unsigned int j=0; j < grillado.indices_de_voxeles_conectados[i].size() ; j++ )
		{
			os 	<< tabbing3 << "\t<indice>" << (grillado.indices_de_voxeles_conectados[i])[j] << "</indice>" << std::endl;
		}
		os << tabbing3 << "</ID de los voxeles conectados>" << std::endl;
	}

	os << tabbing2 << "</conexiones>" << std::endl;

 return os;
}

/**
 * @brief Validates whether a position is within the allowed spatial domain.
 * 
 * Checks if a given 3D coordinate is within the bounding box of the grid.
 * In cancer simulations, this function helps ensure cells remain within
 * the valid computational domain during movement and interactions.
 * 
 * Note: The current implementation returns true for all positions, allowing
 * cells to extend beyond the formal grid boundaries. This can be useful for
 * modeling cancer invasion beyond the initial tumor region.
 * 
 * @param x X-coordinate to validate
 * @param y Y-coordinate to validate
 * @param z Z-coordinate to validate
 * @return Boolean indicating whether the position is valid (currently always true)
 */
bool Grillado_General::es_valida_la_posicion(double x, double y, double z)
{
//	if(x< caja[indice_de_x_min] || x>caja[indice_de_x_max])
//		return false;
//	if(y< caja[indice_de_y_min] || y>caja[indice_de_y_max])
//		return false;
//	if(z< caja[indice_de_z_min] || z>caja[indice_de_z_max])
//		return false;
	return true;
}

/**
 * @brief Creates a bidirectional connection between two voxels.
 * 
 * Establishes neighbor relationships between voxels, which is essential for:
 * - Calculating biochemical diffusion in the tumor microenvironment
 * - Determining cell-cell interactions across voxel boundaries
 * - Simulating gradient-driven processes like chemotaxis and haptotaxis
 * - Modeling cell migration pathways in cancer invasion
 * 
 * In cancer research, voxel connectivity affects how biochemical signals propagate
 * through the tumor microenvironment, influencing cancer cell behavior and organization.
 * 
 * @param i Index of the first voxel to connect
 * @param j Index of the second voxel to connect
 */
void Grillado_General::conectar_voxeles(int i,int j) // done
{

	indices_de_voxeles_conectados[i].push_back( j );
	indices_de_voxeles_conectados[j].push_back( i );

	return;
}

/**
 * @brief Displays general information about the grid to an output stream.
 * 
 * Outputs key grid properties including:
 * - Grid type
 * - Spatial domain dimensions with units
 * - Number of voxels
 * - Total volume of the grid
 * 
 * This function is valuable in cancer research for:
 * - Verifying the domain size matches the intended tumor microenvironment scale
 * - Confirming voxel count and resolution for computational analyses
 * - Documenting grid properties in research reports and publications
 * - Debugging spatial components of cancer simulations
 * 
 * @param os The output stream where grid information will be written
 */
void Grillado_General::mostrar_informacion_general( std::ostream& os){


	os << std::endl << "Informacion del Grillado: " << std::endl
	<< "tipo: grillado general" << std::endl
	<< "Dominio: "
	<< "[" << caja[0] << "," << caja[3] << "] " <<  unidades << " x "
	<< "[" << caja[1] << "," << caja[4] << "] " <<  unidades << " x "
	<< "[" << caja[2] << "," << caja[5] << "] " <<  unidades << std::endl
	<< "   voxeles: " << voxeles.size() << std::endl
	<< "   volumen: ";

	double volumen_total = 0.0;
	for( unsigned int i=0; i < voxeles.size(); i++ )
	{ volumen_total += voxeles[i].volumen; }
	os << volumen_total << " " << unidades << " cubicos " << std::endl;

	return;

}


///////////////////////// CARTESIANO ///////////////////////////////////

/**
 * @brief Default constructor for the Grillado_Cartesiano class.
 * 
 * Initializes a minimal Cartesian grid with a single voxel.
 * This creates a basic spatial structure for testing or as a starting point
 * for more complex grid configurations in cancer simulations.
 * 
 * In cancer research, properly configured Cartesian grids are essential for:
 * - Accurately resolving tumor microenvironment heterogeneity
 * - Modeling diffusion-limited hypoxic regions within tumors
 * - Representing spatial variations in cell phenotypes
 * - Supporting agent-based modeling of individual cancer cells
 * 
 * @return Newly constructed Grillado_Cartesiano object with minimal configuration
 */
Grillado_Cartesiano::Grillado_Cartesiano()
{

	coordenadas_x.assign( 1 , 0.0 );
	coordenadas_y.assign( 1 , 0.0 );
	coordenadas_z.assign( 1 , 0.0 );

	dx = caja[3] - caja[0];
	dy = caja[4] - caja[1];
	dz = caja[5] - caja[2];

	dV = dx*dy*dz;
	dS = dx*dy;
	dS_xy = dx*dy;
	dS_yz = dy*dz;
	dS_xz = dx*dz;

	Voxel template_voxel;
	template_voxel.volumen = dV;

	voxeles.assign( coordenadas_x.size() * coordenadas_y.size() * coordenadas_z.size() , template_voxel );
	voxeles[0].centro.x = coordenadas_x[0];
	voxeles[0].centro.y = coordenadas_y[0];
	voxeles[0].centro.z = coordenadas_z[0];
}

/**
 * @brief Parameterized constructor for the Grillado_Cartesiano class.
 * 
 * Creates a Cartesian grid with the specified number of nodes in each dimension.
 * This constructor establishes a uniform grid with unit spacing, which serves as
 * the foundation for spatial discretization in cancer simulations.
 * 
 * The resolution of the grid directly affects:
 * - Accuracy of biochemical gradient calculations in tumor models
 * - Spatial resolution of oxygen and nutrient distributions
 * - Computational requirements for the simulation
 * - Ability to represent heterogeneous regions within tumors
 * 
 * @param xnodos Number of nodes in the x-dimension
 * @param ynodos Number of nodes in the y-dimension
 * @param znodos Number of nodes in the z-dimension
 * @return Newly constructed Grillado_Cartesiano with specified dimensions
 */
Grillado_Cartesiano::Grillado_Cartesiano( int xnodos, int ynodos, int znodos )
{
	coordenadas_x.assign( xnodos , 0.0 );
	coordenadas_y.assign( ynodos , 0.0 );
	coordenadas_z.assign( znodos , 0.0 );

	dx = 1;
	dy = 1;
	dz = 1;

	dV = dx*dx*dz;
	dS = dx*dy;

	dS_xy = dS;
	dS_yz = dS;
	dS_xz = dS;

	for( unsigned int i=0; i < coordenadas_x.size() ; i++ )
	{ coordenadas_x[i] = i*dx; }
	for( unsigned int i=0; i < coordenadas_y.size() ; i++ )
	{ coordenadas_y[i] = i*dy; }
	for( unsigned int i=0; i < coordenadas_z.size() ; i++ )
	{ coordenadas_z[i] = i*dz; }

	caja[0] = coordenadas_x[0]-dx/2.0;
	caja[3] = coordenadas_x[coordenadas_x.size()-1]+dx/2.0;
	caja[1] = coordenadas_y[0]-dy/2.0;
	caja[4] = coordenadas_y[coordenadas_y.size()-1]+dy/2.0;
	caja[2] = coordenadas_z[0]-dz/2.0;
	caja[5] = coordenadas_z[coordenadas_z.size()-1]+dz/2.0;

	Voxel template_voxel;
	template_voxel.volumen = dV;

	unidades = "micrometros";

	voxeles.assign( coordenadas_x.size() * coordenadas_y.size() * coordenadas_z.size() , template_voxel );



	int n=0;
	for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
	{
		for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
		{
			for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
			{
				voxeles[n].centro.x = coordenadas_x[i];
				voxeles[n].centro.y = coordenadas_y[j];
				voxeles[n].centro.z = coordenadas_z[k];
				voxeles[n].indice = n;
				voxeles[n].volumen = dV;

				n++;
			}
		}
	}



	indices_de_voxeles_conectados.resize( voxeles.size() );

	int i_salto = 1;
	int j_salto = coordenadas_x.size();
	int k_salto = coordenadas_x.size() * coordenadas_y.size();

	// x-aligned connections
	for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
	{
		for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
		{
			for( unsigned int i=0 ; i < coordenadas_x.size()-1 ; i++ )
			{
				int n = g_indice_de_voxel(i,j,k);
				conectar_voxeles(n,n+i_salto);
			}
		}
	}

	for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
	{
		for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
		{
			for( unsigned int j=0 ; j < coordenadas_y.size()-1 ; j++ )
			{
				int n = g_indice_de_voxel(i,j,k);
				conectar_voxeles(n,n+j_salto);
			}
		}
	}

	for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
	{
		for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
		{
			for( unsigned int k=0 ; k < coordenadas_z.size()-1 ; k++ )
			{
				int n = g_indice_de_voxel(i,j,k);
				conectar_voxeles(n,n+k_salto);
			}
		}
	}

}

/**
 * @brief Creates Moore neighborhood connections between voxels.
 * 
 * Establishes connections between each voxel and its 26 adjacent neighbors (Moore neighborhood)
 * within the valid grid boundaries. These connections define the diffusion pathways and
 * interaction networks within the cancer microenvironment simulation.
 * 
 * In cancer research, Moore neighborhoods are critical for:
 * - Creating realistic 3D diffusion patterns for oxygen and nutrients
 * - Modeling biochemical signal propagation through tumor tissue
 * - Enabling diagonal cell migration paths important for cancer invasion
 * - Simulating interactions between cancer cells across voxel boundaries
 * 
 * Note: This function respects domain boundaries and does not establish periodic connections.
 */
void Grillado_Cartesiano::crear_vecindario_moore()
{
	indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
	for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
	{
		for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
		{
			for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
			{
				int indice_centro = g_indice_de_voxel(i,j,k);
				for(int ii=-1;ii<=1;ii++)
					for(int jj=-1;jj<=1;jj++)
						for(int kk=-1;kk<=1;kk++)
							if(i+ii>=0 && i+ii<coordenadas_x.size() &&
								j+jj>=0 && j+jj<coordenadas_y.size() &&
								k+kk>=0 && k+kk<coordenadas_z.size() &&
								!(ii==0 && jj==0 && kk==0))
								{
									int indice_vecino= g_indice_de_voxel(i+ii,j+jj,k+kk);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
								}
			}
		}
	}
}

/**
 * @brief Creates optimized Moore neighborhood connections between voxels.
 * 
 * Establishes an optimized implementation of the 26-neighbor Moore connectivity pattern.
 * This optimized version prioritizes performance for large-scale cancer simulations
 * while maintaining the same connectivity patterns as the standard version.
 * 
 * The optimization is particularly valuable in cancer research for:
 * - Improving computational efficiency in large tumor simulations
 * - Enabling higher resolution tumor microenvironment models
 * - Supporting more complex diffusion-reaction systems with lower overhead
 * - Allowing for more extensive spatial heterogeneity in cancer models
 * 
 * Note: This function respects domain boundaries and does not establish periodic connections.
 */
void Grillado_Cartesiano::crear_vecindario_moore_optimizado()
{
	indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
	for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
	{
		for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
		{
			for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
			{
				int indice_centro = g_indice_de_voxel(i,j,k);
                for(int kk=0;kk<=1;kk++)
                    for(int ii=-1;ii<=1;ii++)
                        if(kk==0){
                            for(int jj=0;jj<=1;jj++){
                                if(i+ii>=0 && i+ii<coordenadas_x.size() &&
								j+jj>=0 && j+jj<coordenadas_y.size() &&
								k+kk>=0 && k+kk<coordenadas_z.size() &&
								!(ii==0 && jj==0 && kk==0))
								{
                                    if(jj==0 && ii==-1){ continue;}else{
									int indice_vecino= g_indice_de_voxel(i+ii,j+jj,k+kk);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }
								}
                                }
                            }else{
                                for(int jj=-1;jj<=1;jj++){
                                    if(i+ii>=0 && i+ii<coordenadas_x.size() &&
                                        j+jj>=0 && j+jj<coordenadas_y.size() &&
                                        k+kk>=0 && k+kk<coordenadas_z.size() &&
                                        !(ii==0 && jj==0 && kk==0)){
                                        
                                        int indice_vecino= g_indice_de_voxel(i+ii,j+jj,k+kk);
                                        indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    
                                        }
                                    }
                                }
                
 /*               
				for(int ii=-1;ii<=1;ii++)
					for(int jj=0;jj<=1;jj++)
						for(int kk=0;kk<=1;kk++)
							if(i+ii>=0 && i+ii<coordenadas_x.size() &&
								j+jj>=0 && j+jj<coordenadas_y.size() &&
								k+kk>=0 && k+kk<coordenadas_z.size() &&
								!(ii==0 && jj==0 && kk==0))
								{
                                    if(jj==0 && ii==-1){ continue;}else{
									int indice_vecino= g_indice_de_voxel(i+ii,j+jj,k+kk);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }
								}
*/								
			}
		}
	}
	
}

void Grillado_Cartesiano::crear_vecindario_moore_periodico()
{
    unsigned int menosuno = -1;
    if(pg->condiciones_de_periodicidad_x && pg->condiciones_de_periodicidad_y && pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=-1;jj<=1;jj++)
                            for(int kk=-1;kk<=1;kk++)
                                if(!(ii==0 && jj==0 && kk==0))
								{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

                                    if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
									if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}
									if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

									if(i+ii == coordenadas_x.size()){c_x=0;}
									if(j+jj == coordenadas_y.size()){c_y=0;}
									if(k+kk == coordenadas_z.size()){c_z=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);

								}
                }
            }
        }
    }


        if(pg->condiciones_de_periodicidad_x && pg->condiciones_de_periodicidad_y && !pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=-1;jj<=1;jj++)
                            for(int kk=-1;kk<=1;kk++)
                                if(k+kk>=0 && k+kk<coordenadas_z.size() &&
                                    !(ii==0 && jj==0 && kk==0))
								{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

                                    if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
									if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}

									if(i+ii == coordenadas_x.size()){c_x=0;}
									if(j+jj == coordenadas_y.size()){c_y=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);

								}
                }
            }
        }
    }

        if(pg->condiciones_de_periodicidad_x && !pg->condiciones_de_periodicidad_y && !pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=-1;jj<=1;jj++)
                            for(int kk=-1;kk<=1;kk++)
                                if(j+jj>=0 && j+jj<coordenadas_y.size() &&
                                   k+kk>=0 && k+kk<coordenadas_z.size() &&
                                   !(ii==0 && jj==0 && kk==0))
								{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

                                    if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
									if(i+ii == coordenadas_x.size()){c_x=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);

								}
                }
            }
        }
    }

        if(!pg->condiciones_de_periodicidad_x && pg->condiciones_de_periodicidad_y && !pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=-1;jj<=1;jj++)
                            for(int kk=-1;kk<=1;kk++)
                                if(i+ii>=0 && i+ii<coordenadas_x.size() &&
                                   k+kk>=0 && k+kk<coordenadas_z.size() &&
                                   !(ii==0 && jj==0 && kk==0))
								{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

									if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}

									if(j+jj == coordenadas_y.size()){c_y=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);

								}
                }
            }
        }
    }

        if(!pg->condiciones_de_periodicidad_x && !pg->condiciones_de_periodicidad_y && pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=-1;jj<=1;jj++)
                            for(int kk=-1;kk<=1;kk++)
                                if(i+ii>=0 && i+ii<coordenadas_x.size() &&
                                   j+jj>=0 && j+jj<coordenadas_y.size() &&
                                   !(ii==0 && jj==0 && kk==0))
								{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

									if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

									if(k+kk == coordenadas_z.size()){c_z=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);

								}
                }
            }
        }
    }

        if(!pg->condiciones_de_periodicidad_x && pg->condiciones_de_periodicidad_y && pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=-1;jj<=1;jj++)
                            for(int kk=-1;kk<=1;kk++)
                                if(i+ii>=0 && i+ii<coordenadas_x.size() &&
                                   !(ii==0 && jj==0 && kk==0))
								{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

									if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}
									if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

									if(j+jj == coordenadas_y.size()){c_y=0;}
									if(k+kk == coordenadas_z.size()){c_z=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);

								}
                }
            }
        }
    }

        if(pg->condiciones_de_periodicidad_x && !pg->condiciones_de_periodicidad_y && pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=-1;jj<=1;jj++)
                            for(int kk=-1;kk<=1;kk++)
                                if(j+jj>=0 && j+jj<coordenadas_y.size() &&
                                  !(ii==0 && jj==0 && kk==0))
								{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

                                    if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
									if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

									if(i+ii == coordenadas_x.size()){c_x=0;}
									if(k+kk == coordenadas_z.size()){c_z=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);

								}
                }
            }
        }
    }



}

void Grillado_Cartesiano::crear_vecindario_moore_periodico_optimizado()
{
    unsigned int menosuno = -1;
    if(pg->condiciones_de_periodicidad_x && pg->condiciones_de_periodicidad_y && pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int kk=0;kk<=1;kk++)
                        for(int ii=-1;ii<=1;ii++)
                            if(kk==0){
                                for(int jj=0;jj<=1;jj++)
                                    if(!(ii==0 && jj==0 && kk==0)){
                                        if(jj==0 && ii==-1){ continue;}else{
                                            int c_x = i+ii;
                                            int c_y = j+jj;
                                            int c_z = k+kk;

                                            if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
                                            if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}
                                            if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

                                            if(i+ii == coordenadas_x.size()){c_x=0;}
                                            if(j+jj == coordenadas_y.size()){c_y=0;}
                                            if(k+kk == coordenadas_z.size()){c_z=0;}

                                            int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                            indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                        }

                                    }
                            }else{
                                for(int jj=-1;jj<=1;jj++)
                                    if(!(ii==0 && jj==0 && kk==0)){
                                        int c_x = i+ii;
                                        int c_y = j+jj;
                                        int c_z = k+kk;

                                        if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
                                        if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}
                                        if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

                                        if(i+ii == coordenadas_x.size()){c_x=0;}
                                        if(j+jj == coordenadas_y.size()){c_y=0;}
                                        if(k+kk == coordenadas_z.size()){c_z=0;}

                                        int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                        indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }
                            }
                                
                        
/*                            
                    
                    
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=0;jj<=1;jj++)
                            for(int kk=0;kk<=1;kk++)
                                if(!(ii==0 && jj==0 && kk==0))
								{
                                    if(jj==0 && ii==-1){ continue;}else{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

                                    if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
									if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}
									if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

									if(i+ii == coordenadas_x.size()){c_x=0;}
									if(j+jj == coordenadas_y.size()){c_y=0;}
									if(k+kk == coordenadas_z.size()){c_z=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }

								}
								
*/
                }
            }
        }
    }
/////////////////////////COND 2 ///////////////////////////////////

        if(pg->condiciones_de_periodicidad_x && pg->condiciones_de_periodicidad_y && !pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int kk=0;kk<=1;kk++)
                        for(int ii=-1;ii<=1;ii++)
                            if(kk==0){
                                for(int jj=0;jj<=1;jj++)                    
                                    if(k+kk>=0 && k+kk<coordenadas_z.size() &&
                                        !(ii==0 && jj==0 && kk==0))
                                    {
                                        if(jj==0 && ii==-1){ continue;}else{
                                            int c_x = i+ii;
                                            int c_y = j+jj;
                                            int c_z = k+kk;

                                            if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
                                            if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}
    
                                            if(i+ii == coordenadas_x.size()){c_x=0;}
                                            if(j+jj == coordenadas_y.size()){c_y=0;}

                                            int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                            indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                        }

                                    }
                            }else{
                                for(int jj=-1;jj<=1;jj++)                    
                                    if(k+kk>=0 && k+kk<coordenadas_z.size() &&
                                        !(ii==0 && jj==0 && kk==0))
                                    {
                                        if(jj==0 && ii==-1){ continue;}else{
                                            int c_x = i+ii;
                                            int c_y = j+jj;
                                            int c_z = k+kk;

                                            if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
                                            if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}
    
                                            if(i+ii == coordenadas_x.size()){c_x=0;}
                                            if(j+jj == coordenadas_y.size()){c_y=0;}

                                            int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                            indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                        }

                                    }
                                
                            }
                    
                    
/*                    
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=0;jj<=1;jj++)
                            for(int kk=0;kk<=1;kk++)
                                if(k+kk>=0 && k+kk<coordenadas_z.size() &&
                                    !(ii==0 && jj==0 && kk==0))
								{
                                    if(jj==0 && ii==-1){ continue;}else{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

                                    if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
									if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}

									if(i+ii == coordenadas_x.size()){c_x=0;}
									if(j+jj == coordenadas_y.size()){c_y=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }

								}
								
*/								
                }
            }
        }
    }
    
////////////////////COND 3 ///////////////////

        if(pg->condiciones_de_periodicidad_x && !pg->condiciones_de_periodicidad_y && !pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int kk=0;kk<=1;kk++)
                        for(int ii=-1;ii<=1;ii++)
                            if(kk==0){
                                for(int jj=0;jj<=1;jj++)  
                                    if(j+jj>=0 && j+jj<coordenadas_y.size() &&
                                    k+kk>=0 && k+kk<coordenadas_z.size() &&
                                    !(ii==0 && jj==0 && kk==0)){
                                        if(jj==0 && ii==-1){ continue;}else{
                                            int c_x = i+ii;
                                            int c_y = j+jj;
                                            int c_z = k+kk;

                                            if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
                                            if(i+ii == coordenadas_x.size()){c_x=0;}

                                            int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                            indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                        }
                                        
                                    }
                            }else{
                                for(int jj=-1;jj<=1;jj++)  
                                    if(j+jj>=0 && j+jj<coordenadas_y.size() &&
                                    k+kk>=0 && k+kk<coordenadas_z.size() &&
                                    !(ii==0 && jj==0 && kk==0)){
                                        if(jj==0 && ii==-1){ continue;}else{
                                            int c_x = i+ii;
                                            int c_y = j+jj;
                                            int c_z = k+kk;

                                            if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
                                            if(i+ii == coordenadas_x.size()){c_x=0;}

                                            int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                            indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                        }
                                        
                                    }                                
                            }
/*                    
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=0;jj<=1;jj++)
                            for(int kk=0;kk<=1;kk++)
                                if(j+jj>=0 && j+jj<coordenadas_y.size() &&
                                   k+kk>=0 && k+kk<coordenadas_z.size() &&
                                   !(ii==0 && jj==0 && kk==0))
								{
                                    if(jj==0 && ii==-1){ continue;}else{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

                                    if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
									if(i+ii == coordenadas_x.size()){c_x=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }

								}
*/
                }
            }
        }
    }

 ///////////// COND 4 ///////////////////   
    
        if(!pg->condiciones_de_periodicidad_x && pg->condiciones_de_periodicidad_y && !pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int kk=0;kk<=1;kk++)
                        for(int ii=-1;ii<=1;ii++)
                            if(kk==0){
                                for(int jj=0;jj<=1;jj++)  
                                    if(i+ii>=0 && i+ii<coordenadas_x.size() &&
                                    k+kk>=0 && k+kk<coordenadas_z.size() &&
                                    !(ii==0 && jj==0 && kk==0)){
                                        if(jj==0 && ii==-1){ continue;}else{
                                            int c_x = i+ii;
                                            int c_y = j+jj;
                                            int c_z = k+kk;

                                            if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}

                                            if(j+jj == coordenadas_y.size()){c_y=0;}

                                            int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                            indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }

								}
                            }else{
                                for(int jj=-1;jj<=1;jj++)  
                                    if(i+ii>=0 && i+ii<coordenadas_x.size() &&
                                        k+kk>=0 && k+kk<coordenadas_z.size() &&
                                        !(ii==0 && jj==0 && kk==0)){
                                            if(jj==0 && ii==-1){ continue;}else{
                                                int c_x = i+ii;
                                                int c_y = j+jj;
                                                int c_z = k+kk;

                                                if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}

                                                if(j+jj == coordenadas_y.size()){c_y=0;}

                                                int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                                indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }

								}                 
                            }
                                    
/*                                    
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=0;jj<=1;jj++)
                            for(int kk=0;kk<=1;kk++)
                                if(i+ii>=0 && i+ii<coordenadas_x.size() &&
                                   k+kk>=0 && k+kk<coordenadas_z.size() &&
                                   !(ii==0 && jj==0 && kk==0))
								{
                                    if(jj==0 && ii==-1){ continue;}else{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

									if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}

									if(j+jj == coordenadas_y.size()){c_y=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }

								}
*/
                }
            }
        }
    }
    
 ///////////// COND 5 ///////////////////   

        if(!pg->condiciones_de_periodicidad_x && !pg->condiciones_de_periodicidad_y && pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int kk=0;kk<=1;kk++)
                        for(int ii=-1;ii<=1;ii++)
                            if(kk==0){
                                for(int jj=0;jj<=1;jj++)
                                    if(i+ii>=0 && i+ii<coordenadas_x.size() &&
                                    j+jj>=0 && j+jj<coordenadas_y.size() &&
                                    !(ii==0 && jj==0 && kk==0)){
                                        if(jj==0 && ii==-1){ continue;}else{
                                            int c_x = i+ii;
                                            int c_y = j+jj;
                                            int c_z = k+kk;

                                            if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

                                            if(k+kk == coordenadas_z.size()){c_z=0;}

                                            int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                            indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }

								}
                            }else{
                                for(int jj=-1;jj<=1;jj++)
                                    if(i+ii>=0 && i+ii<coordenadas_x.size() &&
                                    j+jj>=0 && j+jj<coordenadas_y.size() &&
                                    !(ii==0 && jj==0 && kk==0)){
                                        if(jj==0 && ii==-1){ continue;}else{
                                            int c_x = i+ii;
                                            int c_y = j+jj;
                                            int c_z = k+kk;

                                            if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

                                            if(k+kk == coordenadas_z.size()){c_z=0;}

                                            int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                            indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }

								}
                            }
                    
/*                    
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=0;jj<=1;jj++)
                            for(int kk=0;kk<=1;kk++)
                                if(i+ii>=0 && i+ii<coordenadas_x.size() &&
                                   j+jj>=0 && j+jj<coordenadas_y.size() &&
                                   !(ii==0 && jj==0 && kk==0))
								{
                                    if(jj==0 && ii==-1){ continue;}else{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

									if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

									if(k+kk == coordenadas_z.size()){c_z=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }

								}
*/
                }
            }
        }
    }
    
///////////// COND 6 ///////////////////    

        if(!pg->condiciones_de_periodicidad_x && pg->condiciones_de_periodicidad_y && pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int kk=0;kk<=1;kk++)
                        for(int ii=-1;ii<=1;ii++)
                            if(kk==0){
                                for(int jj=0;jj<=1;jj++)
                                    if(i+ii>=0 && i+ii<coordenadas_x.size() &&
                                    !(ii==0 && jj==0 && kk==0)){
                                        if(jj==0 && ii==-1){ continue;}else{
                                            int c_x = i+ii;
                                            int c_y = j+jj;
                                            int c_z = k+kk;

                                            if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}
                                            if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

                                            if(j+jj == coordenadas_y.size()){c_y=0;}
                                            if(k+kk == coordenadas_z.size()){c_z=0;}

                                            int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                            indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                        }

                                    }
                            }else{
                                for(int jj=-1;jj<=1;jj++)
                                    if(i+ii>=0 && i+ii<coordenadas_x.size() &&
                                    !(ii==0 && jj==0 && kk==0)){
                                        if(jj==0 && ii==-1){ continue;}else{
                                            int c_x = i+ii;
                                            int c_y = j+jj;
                                            int c_z = k+kk;

                                            if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}
                                            if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

                                            if(j+jj == coordenadas_y.size()){c_y=0;}
                                            if(k+kk == coordenadas_z.size()){c_z=0;}

                                            int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                            indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                        }

                                    }                                
                            }
                            
/*
                    
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=0;jj<=1;jj++)
                            for(int kk=0;kk<=1;kk++)
                                if(i+ii>=0 && i+ii<coordenadas_x.size() &&
                                   !(ii==0 && jj==0 && kk==0))
								{
                                    if(jj==0 && ii==-1){ continue;}else{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

									if(j+jj==menosuno){c_y=coordenadas_y.size()-1;}
									if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

									if(j+jj == coordenadas_y.size()){c_y=0;}
									if(k+kk == coordenadas_z.size()){c_z=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }

								}
*/
                }
            }
        }
    }

///////////// COND 7 ///////////////////    
    
        if(pg->condiciones_de_periodicidad_x && !pg->condiciones_de_periodicidad_y && pg->condiciones_de_periodicidad_z){

        indices_de_voxeles_conectados_tipo_moore.resize( voxeles.size() );
        for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
        {
            for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
            {
                for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
                {
                    int indice_centro = g_indice_de_voxel(i,j,k);
                    for(int kk=0;kk<=1;kk++)
                        for(int ii=-1;ii<=1;ii++)
                            if(kk==0){
                                for(int jj=0;jj<=1;jj++)
                                    if(j+jj>=0 && j+jj<coordenadas_y.size() &&
                                    !(ii==0 && jj==0 && kk==0)){
                                        if(jj==0 && ii==-1){ continue;}else{
                                        int c_x = i+ii;
                                        int c_y = j+jj;
                                        int c_z = k+kk;

                                        if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
                                        if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

                                        if(i+ii == coordenadas_x.size()){c_x=0;}
                                        if(k+kk == coordenadas_z.size()){c_z=0;}

                                        int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                        indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                        }
                                    }
                            }else{
                                for(int jj=-1;jj<=1;jj++)
                                    if(j+jj>=0 && j+jj<coordenadas_y.size() &&
                                    !(ii==0 && jj==0 && kk==0)){
                                        if(jj==0 && ii==-1){ continue;}else{
                                        int c_x = i+ii;
                                        int c_y = j+jj;
                                        int c_z = k+kk;

                                        if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
                                        if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

                                        if(i+ii == coordenadas_x.size()){c_x=0;}
                                        if(k+kk == coordenadas_z.size()){c_z=0;}

                                        int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
                                        indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                        }
                                    }                                
                            }

								
/*								
                    
                    for(int ii=-1;ii<=1;ii++)
                        for(int jj=0;jj<=1;jj++)
                            for(int kk=0;kk<=1;kk++)
                                if(j+jj>=0 && j+jj<coordenadas_y.size() &&
                                  !(ii==0 && jj==0 && kk==0))
								{
                                    if(jj==0 && ii==-1){ continue;}else{
                                    int c_x = i+ii;
                                    int c_y = j+jj;
                                    int c_z = k+kk;

                                    if(i+ii==menosuno){c_x=coordenadas_x.size()-1;}
									if(k+kk==menosuno){c_z=coordenadas_z.size()-1;}

									if(i+ii == coordenadas_x.size()){c_x=0;}
									if(k+kk == coordenadas_z.size()){c_z=0;}

                                    int indice_vecino= g_indice_de_voxel(c_x,c_y,c_z);
									indices_de_voxeles_conectados_tipo_moore[indice_centro].push_back( indice_vecino);
                                    }

								}
*/								
                }
            }
        }
    }



}


void Grillado_Cartesiano::mostrar_vecindarios_moore(std::ostream& os){

	std::boolalpha( os );
	static std::string tabbing = "\t";
	static std::string tabbing2 = "\t\t";
	static std::string tabbing3 = "\t\t\t";

	for( unsigned int i=0; i < voxeles.size() ; i++ )
	{ os << voxeles[i] << std::endl; }
	os 	<< tabbing2 << "</voxeles>" << std::endl;

	os	<< tabbing2 << "<conexiones>" << std::endl;
	for( unsigned int i=0 ; i < indices_de_voxeles_conectados_tipo_moore.size() ; i++ )
	{
		os << tabbing3 << "<ID de los voxeles conectados=\"" << i << "\">" << std::endl;
		for( unsigned int j=0; j < indices_de_voxeles_conectados_tipo_moore[i].size() ; j++ )
		{
			os 	<< tabbing3 << "\t<indice>" << (indices_de_voxeles_conectados_tipo_moore[i])[j] << "</indice>" << std::endl;
		}
		os << tabbing3 << "</ID de los voxeles conectados>" << std::endl;
	}

	os << tabbing2 << "</conexiones>" << std::endl;

 return;
}

/**
 * @brief Converts 3D Cartesian indices to a linear voxel index.
 * 
 * Maps three-dimensional coordinates (i,j,k) to a one-dimensional array index.
 * This function is fundamental for efficient spatial operations in the cancer simulation,
 * enabling quick access to voxel-specific data like biochemical concentrations and cell populations.
 * 
 * In cancer research, this mapping is critical for:
 * - Accessing specific locations within the tumor microenvironment
 * - Implementing efficient diffusion algorithms for oxygen and nutrients
 * - Tracking spatial patterns of cancer cell distribution
 * - Correlating cell behaviors with their positions in the tumor
 * 
 * @param i Index in the x-dimension
 * @param j Index in the y-dimension
 * @param k Index in the z-dimension
 * @return Linear index corresponding to the specified 3D coordinates
 */
unsigned int Grillado_Cartesiano::indice_de_voxel( unsigned int i, unsigned int j, unsigned int k )
{
 return ( k*coordenadas_y.size() + j )*coordenadas_x.size() + i;
}

/**
 * @brief Converts a linear voxel index to 3D Cartesian indices.
 * 
 * The inverse operation of indice_de_voxel, mapping a one-dimensional array index
 * back to three-dimensional coordinates (i,j,k). This function enables spatial analysis
 * and visualization of tumor patterns in cancer research.
 * 
 * This conversion is essential for:
 * - Visualizing spatial patterns in tumor growth
 * - Analyzing regional variations in cancer cell phenotypes
 * - Identifying spatial correlations between microenvironmental factors and cell behaviors
 * - Implementing boundary conditions at specific faces of the simulation domain
 * 
 * @param n Linear index to convert
 * @return Vector of 3 unsigned integers representing the x, y, and z indices
 */
std::vector<unsigned int> Grillado_Cartesiano::indices_cartesianos( unsigned int n )
{
	std::vector<unsigned int> out(3, -1 );


	unsigned int XY = coordenadas_x.size() * coordenadas_y.size();
	out[0] = (unsigned int) floor( n/XY );


	out[1] = (unsigned int) floor(   (n - out[0]*XY) / coordenadas_x.size() );


	out[2] = n - coordenadas_x.size()*(   out[1] + coordenadas_y.size()*out[0] );

	return out;
}

/**
 * @brief Resizes the grid with specified domain boundaries and node counts.
 * 
 * Reconfigures the Cartesian grid with new spatial dimensions and resolution.
 * This method is critical for setting up the appropriate spatial scale and resolution
 * for cancer simulations, allowing researchers to focus on specific tumor regions
 * or adjust computational requirements.
 * 
 * In cancer research, grid dimensioning affects:
 * - The physical scale of the simulated tumor microenvironment
 * - Spatial resolution for capturing tumor heterogeneity
 * - Computational efficiency of diffusion and cell interaction calculations
 * - Ability to represent fine-scale features like blood vessels or tissue boundaries
 * 
 * @param x_ini Minimum x-coordinate of the domain
 * @param x_fin Maximum x-coordinate of the domain
 * @param y_ini Minimum y-coordinate of the domain
 * @param y_fin Maximum y-coordinate of the domain
 * @param z_ini Minimum z-coordinate of the domain
 * @param z_fin Maximum z-coordinate of the domain
 * @param x_nodos Number of nodes in the x-dimension
 * @param y_nodos Number of nodes in the y-dimension
 * @param z_nodos Number of nodes in the z-dimension
 */
void Grillado_Cartesiano::redimensionar( double x_ini, double x_fin, double y_ini, double y_fin, double z_ini, double z_fin , int x_nodos, int y_nodos, int z_nodos )
{
	coordenadas_x.assign( x_nodos , 0.0 );
	coordenadas_y.assign( y_nodos , 0.0 );
	coordenadas_z.assign( z_nodos , 0.0 );

	dx = ( x_fin - x_ini )/( (double) x_nodos );
	if( x_nodos < 2 )
	{ dx = 1; }
	dy = ( y_fin - y_ini )/( (double) y_nodos );
	if( y_nodos < 2 )
	{ dy = 1; }
	dz = ( z_fin - z_ini )/( (double) z_nodos  );
	if( z_nodos < 2 )
	{ dz = 1; }


	for( unsigned int i=0; i < coordenadas_x.size() ; i++ )
	{ coordenadas_x[i] = x_ini + (i+0.5)*dx; }
	for( unsigned int i=0; i < coordenadas_y.size() ; i++ )
	{ coordenadas_y[i] = y_ini + (i+0.5)*dy; }
	for( unsigned int i=0; i < coordenadas_z.size() ; i++ )
	{ coordenadas_z[i] = z_ini + (i+0.5)*dz; }

	caja[0] = x_ini;
	caja[3] = x_fin;
	caja[1] = y_ini;
	caja[4] = y_fin;
	caja[2] = z_ini;
	caja[5] = z_fin;

	dV = dx*dy*dz;
	dS = dx*dy;

	dS_xy = dx*dy;
	dS_yz = dy*dz;
	dS_xz = dx*dz;

	Voxel template_voxel;
	template_voxel.volumen = dV;

	voxeles.assign( coordenadas_x.size() * coordenadas_y.size() * coordenadas_z.size() , template_voxel );

	int n=0;
	for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
	{
		for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
		{
			for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
			{
				voxeles[n].centro.x = coordenadas_x[i];
				voxeles[n].centro.y = coordenadas_y[j];
				voxeles[n].centro.z = coordenadas_z[k];
				voxeles[n].indice = n;
				voxeles[n].volumen = dV;

				n++;
			}
		}
	}



	indices_de_voxeles_conectados.resize( voxeles.size() );

	for( unsigned int i=0; i < indices_de_voxeles_conectados.size() ; i++ )
	{ indices_de_voxeles_conectados[i].clear(); }
	int i_salto = 1;
	int j_salto = coordenadas_x.size();
	int k_salto = coordenadas_x.size() * coordenadas_y.size();


	int cuenta = 0;
	for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
	{
		for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
		{
			for( unsigned int i=0 ; i < coordenadas_x.size()-1 ; i++ )
			{
				int n = g_indice_de_voxel(i,j,k);
				conectar_voxeles(n,n+i_salto);
				cuenta++;
			}
		}
	}

	for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
	{
		for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
		{
			for( unsigned int j=0 ; j < coordenadas_y.size()-1 ; j++ )
			{
				int n = g_indice_de_voxel(i,j,k);
				conectar_voxeles(n,n+j_salto);
			}
		}
	}

	for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
	{
		for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
		{
			for( unsigned int k=0 ; k < coordenadas_z.size()-1 ; k++ )
			{
				int n = g_indice_de_voxel(i,j,k);
				conectar_voxeles(n,n+k_salto);
			}
		}
	}

	crear_vecindario_moore();
	return;
}


void Grillado_Cartesiano::redimensionar( double x_ini, double x_fin, double y_ini, double y_fin, double z_ini, double z_fin , double dx_nuevo, double dy_nuevo, double dz_nuevo)
{
	dx = dx_nuevo;
	dy = dy_nuevo;
	dz = dz_nuevo;

	double eps = 1e-16;
	int x_nodos = (int) ceil( eps + (x_fin-x_ini)/dx );
	int y_nodos = (int) ceil( eps + (y_fin-y_ini)/dy );
	int z_nodos = (int) ceil( eps + (z_fin-z_ini)/dz );

	coordenadas_x.assign( x_nodos , 0.0 );
	coordenadas_y.assign( y_nodos , 0.0 );
	coordenadas_z.assign( z_nodos , 0.0 );

	for( unsigned int i=0; i < coordenadas_x.size() ; i++ )
	{ coordenadas_x[i] = x_ini + (i+0.5)*dx; }
	for( unsigned int i=0; i < coordenadas_y.size() ; i++ )
	{ coordenadas_y[i] = y_ini + (i+0.5)*dy; }
	for( unsigned int i=0; i < coordenadas_z.size() ; i++ )
	{ coordenadas_z[i] = z_ini + (i+0.5)*dz; }

	caja[0] = x_ini;
	caja[3] = x_fin;
	caja[1] = y_ini;
	caja[4] = y_fin;
	caja[2] = z_ini;
	caja[5] = z_fin;

	dV = dx*dy*dz;
	dS = dx*dy;

	dS_xy = dx*dy;
	dS_yz = dy*dz;
	dS_xz = dx*dz;

	Voxel template_voxel;
	template_voxel.volumen = dV;

	voxeles.assign( coordenadas_x.size() * coordenadas_y.size() * coordenadas_z.size() , template_voxel );

	int n=0;
	for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
	{
		for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
		{
			for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
			{
				voxeles[n].centro.x = coordenadas_x[i];
				voxeles[n].centro.y = coordenadas_y[j];
				voxeles[n].centro.z = coordenadas_z[k];
				voxeles[n].indice = n;
				voxeles[n].volumen = dV;

				n++;
			}
		}
	}



	indices_de_voxeles_conectados.resize( voxeles.size() );

	for( unsigned int i=0; i < indices_de_voxeles_conectados.size() ; i++ )
	{ indices_de_voxeles_conectados[i].clear(); }

	int i_salto = 1;
	int j_salto = coordenadas_x.size();
	int k_salto = coordenadas_x.size() * coordenadas_y.size();


	int cuenta = 0;
	for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
	{
		for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
		{
			for( unsigned int i=0 ; i < coordenadas_x.size()-1 ; i++ )
			{
				int n = g_indice_de_voxel(i,j,k);
				conectar_voxeles(n,n+i_salto);
				cuenta++;
			}
		}
	}


	for( unsigned int k=0 ; k < coordenadas_z.size() ; k++ )
	{
		for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
		{
			for( unsigned int j=0 ; j < coordenadas_y.size()-1 ; j++ )
			{
				int n = g_indice_de_voxel(i,j,k);
				conectar_voxeles(n,n+j_salto);
			}
		}
	}


	for( unsigned int j=0 ; j < coordenadas_y.size() ; j++ )
	{
		for( unsigned int i=0 ; i < coordenadas_x.size() ; i++ )
		{
			for( unsigned int k=0 ; k < coordenadas_z.size()-1 ; k++ )
			{
				int n = g_indice_de_voxel(i,j,k);
				conectar_voxeles(n,n+k_salto);
			}
		}
	}

	if (pg->condiciones_de_periodicidad) {
		crear_vecindario_moore_periodico_optimizado();
	} else {
        crear_vecindario_moore_optimizado();
		//crear_vecindario_moore();
	}
	return;
}

void Grillado_Cartesiano::redimensionar( int x_nodos, int y_nodos, int z_nodos )
{ return redimensionar( 0-0.5, x_nodos-1+0.5 , 0-0.5 , y_nodos-1+0.5 , 0-0.5 , z_nodos - 1+0.5 , x_nodos, y_nodos, z_nodos ); }


int Grillado_Cartesiano::indice_del_voxel_mas_cercano( Vector& posicion )
{
	unsigned int i = (unsigned int) floor( (posicion.x-caja[0])/dx );
	unsigned int j = (unsigned int) floor( (posicion.y-caja[1])/dy );
	unsigned int k = (unsigned int) floor( (posicion.z-caja[2])/dz );



	if( i >= coordenadas_x.size() ){ i = coordenadas_x.size()-1; }
	if( i < 0 ){ i = 0; }

	if( j >= coordenadas_y.size() ){ j = coordenadas_y.size()-1; }
	if( j < 0 ){ j = 0; }

	if( k >= coordenadas_z.size() ){ k = coordenadas_z.size()-1; }
	if( k < 0 ){ k = 0; }

	return ( k*coordenadas_y.size() + j )*coordenadas_x.size() + i;
}


Vector Grillado_Cartesiano::indices_cartesianos_mas_cercanos( Vector& posicion )
{
	Vector out;
	out.x = (unsigned int) floor( (posicion.x-caja[0])/dx );
	out.y = (unsigned int) floor( (posicion.y-caja[1])/dy );
	out.z = (unsigned int) floor( (posicion.z-caja[2])/dz );



	if( out.x >= coordenadas_x.size() ){ out.x = coordenadas_x.size()-1; }
	if( out.x < 0 ){ out.x = 0; }

	if( out.y >= coordenadas_y.size() ){ out.y = coordenadas_y.size()-1; }
	if( out.y < 0 ){ out.y = 0; }

	if( out.z >= coordenadas_z.size() ){ out.z = coordenadas_z.size()-1; }
	if( out.z < 0 ){ out.z = 0; }

	return out;
}


Voxel& Grillado_Cartesiano::voxel_mas_cercano( Vector& posicion ){
 return voxeles[ indice_del_voxel_mas_cercano( posicion ) ];
}

/**
 * @brief Gets the center coordinates of a voxel by its index.
 * 
 * Returns the 3D spatial coordinates of the center of a specified voxel.
 * This function is essential for mapping between discrete voxel indices and
 * continuous spatial coordinates in the cancer simulation.
 * 
 * In cancer research, this mapping is important for:
 * - Positioning cells accurately within the tumor microenvironment
 * - Calculating distances between cells and biochemical sources
 * - Visualizing spatial patterns in tumor growth
 * - Implementing boundary conditions at specific locations
 * 
 * @param indice_de_voxel Linear index of the voxel
 * @return Vector containing the x, y, z coordinates of the voxel center
 */
Vector Grillado_Cartesiano::get_centro_voxel(int indice_de_voxel){

    Vector centro;
    centro = voxeles[ indice_de_voxel ].centro;
    return centro;


}

/**
 * @brief Displays detailed information about the Cartesian grid to an output stream.
 * 
 * Outputs comprehensive grid properties including:
 * - Grid type and dimensions
 * - Spatial domain boundaries with units
 * - Voxel counts in each dimension
 * - Voxel sizes (dx, dy, dz)
 * - Total volume and surface areas
 * 
 * This function is valuable in cancer research for:
 * - Verifying grid configuration for tumor microenvironment modeling
 * - Documenting simulation parameters for reproducibility
 * - Debugging spatial components of cancer simulations
 * - Reporting grid properties in research publications
 * 
 * @param os The output stream where grid information will be written
 */
void Grillado_Cartesiano::mostrar_informacion_cartesiano( std::ostream& os )
{
	os << std::endl << "Informacion del Grillado: " << std::endl;
	os << "tipo: Cartesiano" << std::endl;
	os << "Dominio: "
	<< "[" << caja[0] << "," << caja[3] << "] " <<  unidades << " x "
	<< "[" << caja[1] << "," << caja[4] << "] " <<  unidades << " x "
	<< "[" << caja[2] << "," << caja[5] << "] " <<  unidades << std::endl
	<< "   resolucion: dx = " << dx << " " << unidades;
	os << std::endl
	<< "   voxeles: " << voxeles.size() << std::endl
	<< "   volumen: " << ( caja[3]-caja[0] )*( caja[4]-caja[1] )*( caja[5]-caja[2] )
		<< " " << unidades << " cubicos " << std::endl;

	return;
}


