/**
 * @file: Main.cpp
 *
 * @author: Luciana Melina Luque
 *
 * @details Entry point and main simulation control for the cell-based cancer model.
 * This file implements the primary simulation loop and orchestrates the three main processes:
 * diffusion, mechanics, and cell biology.
 * 
 * Inbound Dependencies:
 * - Windows.h/unistd.h - Platform-specific timing functions
 * - iostream, fstream - I/O operations
 * - cstdio, cstdlib, ctime, cmath - Standard C++ utilities
 * - iomanip, chrono - Formatting and timing utilities
 * - Macros.h - Constants and utility macros
 * - Tejido.h - Tissue implementation (which includes all other components)
 * - Ciclos_estandares.h - Standard cell cycle models
 * 
 * Outbound Dependencies:
 * - None (this is the top-level control file)
 * 
 * Usage:
 * The program is executed with a parameter file path:
 * ./program parameter_file.txt
 * 
 * @license Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif
#include <fstream>
#include<iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include<iomanip>
#include <chrono>
#include <ctime>
#include <sys/stat.h> // For mkdir
#include <errno.h>    // For error codes
#include <string.h>   // For strerror

#include "Macros.h"
#include "Tejido.h"


/**
 * Global variables used throughout the simulation
 * 
 * @var rng Random number generator for stochastic processes
 * @var v Vector utility object
 * @var c Constants object containing simulation parameters
 * @var Ki67, vida, necrosis, apoptosis Cell cycle model objects
 * @var necrosis_parametros, apoptosis_parametros Death process parameters
 * @var parametros_globales Global parameters object
 * @var pg Pointer to the global parameters
 */
RNG *rng;
Vector *v;
Constantes *c;
Ciclo_Modelo Ki67, vida, necrosis, apoptosis;
Muerte_parametros necrosis_parametros, apoptosis_parametros;

Parametros_globales parametros_globales;
Parametros_globales *pg = &parametros_globales;



using namespace std;

///*
//Para debuggear el código
#define INFO(msg) \
    fprintf(stderr, "info: %s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, "%s\n", msg);
//*/

/**
 * Cell collection vectors used to track different cell states
 * 
 * @var todas_las_celulas Vector containing all cells in the simulation
 * @var celulas_listas_para_dividirse Vector of cells ready to divide
 * @var celulas_listas_para_remover Vector of cells to be removed
 * @var celulas_para_registrar_en_voxeles Vector of cells to register in voxels
 * @var ciclo_celular_estandar_inicializado Flag indicating if standard cell cycle is initialized
 * @var ciclo_celular_de_muerte_inicializado Flag indicating if death cycle is initialized
 * @var datos_finales_T, datos_finales_V Vectors for storing final time and volume data
 */
std::vector<Celula*> todas_las_celulas;
std::vector<Celula*> celulas_listas_para_dividirse;
std::vector<Celula*> celulas_listas_para_remover;
std::vector<Celula*> celulas_para_registrar_en_voxeles;
bool ciclo_celular_estandar_inicializado;
bool ciclo_celular_de_muerte_inicializado;
std::vector<double> datos_finales_T;
std::vector<double> datos_finales_V;

/**
 * @brief Creates directory if it doesn't exist.
 * 
 * Creates a directory at the specified path if it doesn't already exist.
 * 
 * @param path The directory path to create
 * @return true if directory exists or was created successfully, false otherwise
 */
bool create_directory_if_not_exists(const std::string& path) {
    struct stat info;
    
    // Check if directory already exists
    if (stat(path.c_str(), &info) == 0) {
        return true; // Directory exists
    }
    
    // Directory doesn't exist, create it
    #ifdef _WIN32
    int result = mkdir(path.c_str());
    #else
    int result = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    #endif
    
    if (result == 0) {
        std::cout << "Created directory: " << path << std::endl;
        return true;
    }
    
    // Check if directory was created by another process in the meantime
    if (errno == EEXIST) {
        return true;
    }
    
    std::cerr << "Error creating directory " << path << ": " << strerror(errno) << std::endl;
    return false;
}

/**
 * @brief Main entry point for the cell-based cancer simulation.
 * 
 * Controls the overall simulation flow:
 * 1. Loads parameters from input file
 * 2. Initializes the tissue and cell structures
 * 3. Executes the main simulation loop (diffusion, mechanics, biology)
 * 4. Outputs results at specified intervals
 * 
 * @param argc Number of command line arguments
 * @param argv Array of command line arguments, argv[1] should be the parameter file path
 * @return 0 on successful completion
 * 
 * Cancer Research Context:
 * This simulation enables studying cancer cell behavior, tumor growth dynamics,
 * and potential therapeutic interventions in a controlled in silico environment.
 */
int main(int argc, char *argv[]){
    // Create results directory
    const std::string results_dir = "results";
    create_directory_if_not_exists(results_dir);
    
    // Read parameter file
    std::string inputFile(argv[1]);
    std::string parametro;
	std::string valor;



	ifstream param(inputFile);
	if (!param.is_open()) {
		cout << "ERROR: No se puede abrir el archivo " << inputFile << endl;
		exit(1);
	}
	while (!param.eof()) {
		param >> parametro;
		param >> valor;
		cout << parametro << " = " << valor << endl;
		pg->set_parametros(parametro, valor);
	}

    // Initialize random number generator and time settings
	rng = new RNG(1.0, pg->seed, pg->imm_mean, pg->imm_sd);
	pg->tiempo_total=0.0;
	// pg->tiempo_final = 60*24*12+0.01;//Debug
	pg->tiempo_final = 43200.01;//Debug Uncomment
    double tiempo_de_escritura = 720.0;


    // Initialize tissue and cell cycle models
	Tejido tejido;
	crear_ciclo_celular_estandar();
	crear_ciclo_de_muerte_estandar();


	pg->ciclo = Ki67;
	pg->numero_id = 0;


    // Initialize tissue with cells
	tejido.inicializar_tejido();
    
//    for(int j=0; j<20; j++){
//        std::cout << rng->RandomNumber(pg->rango_en_X[0],pg->rango_en_X[1]) << "\n";
//        std::cin.get();
//    }
    
//    cout << "Cantidad de celulas: " << todas_las_celulas.size() << "\n";

//	tejido.microambiente.mostrar_informacion(cout);
//	cout << "pulse enter para continuar \n";
//	tejido.cdc.grillado.mostrar_informacion_cartesiano(cout);
//	cout << "pulse enter para continuar \n";
//	tejido.cdc.celula->mostrar_informacion_de_la_celula(cout);
//	tejido.cdc.celula->fenotipo.ciclo.pCiclo_Modelo->mostrar_ciclo(cout);
//	cout << "pulse enter para continuar \n";


//	cout << "TASAS DE TRANSICION \n";
//	cout << tejido.cdc.celula->fenotipo.ciclo.tasas_de_transicion[0][0] << "\n";
//	cout << tejido.cdc.celula->fenotipo.ciclo.tasas_de_transicion[1][0] << "\n";
//	cout << tejido.cdc.celula->fenotipo.ciclo.tasas_de_transicion[2][0] << "\n";
//	cout << "pulse enter para continuar \n";
//	cin.get();

    /////////////////CHEQUEAR VASO SANGUINEO CON CELULAS////////////////////
/*
	tejido.microambiente.crear_vaso_sanguineo(pg->rango_en_X[0],pg->rango_en_Y[0],pg->rango_en_Z[0],pg->rango_en_X[0]+1,pg->rango_en_Y[0]+1,pg->rango_en_Z[0]+1);

	std::cout << "cantidad de voxeles: " << tejido.microambiente.voxeles_del_vaso_sanguineo.size() << "\n";
    //std::cin.get();
*/

////////////////////////////////////////////////////////////////////////
    tejido.geometria_del_tumor();
    std::ofstream outfile(results_dir + "/DatosFinales.dat", ios::app);
    outfile << "#tiempo, volumen, volumen2, radio , celulas tumorales, dead_cancer_cells, todas las celulas, cart_alive \n";
    outfile.flush();




    // Variables for tracking timing of mechanics, output, and other operations
    double n = 0;  // Last mechanics time
    double m = 0;  // Last output time
    double m2 = 0; // Additional timing variable
    int itr = 0;
    static double tolerancia_de_la_mecanica = 0.001 * c->dt_mecanica;   
    static double tolerancia_de_la_escritura = 0.001 * c->dt_ciclo;
    
    // Main simulation loop
    for (pg->tiempo_total =0; pg->tiempo_total <= pg->tiempo_final; pg->tiempo_total+=c->dt_difusion){
        // Conditional introduction of immune cells at specific time points
        static bool introduje_celulas_inmunes = false;
        if( pg->activar_respuesta_inmune == true && pg->tiempo_total > pg->tiempo_de_imm - 0.01*c->dt_difusion && introduje_celulas_inmunes == false ){
            tejido.introducir_linfocitos_aleatorios(pg->cantidad_de_linfocitos);
            introduje_celulas_inmunes = true;
        }

        static bool introduje_celulas_inmunes2 = false;
        if( pg->activar_respuesta_inmune == true && pg->tiempo_total > pg->tiempo_de_imm_2 - 0.01*c->dt_difusion && introduje_celulas_inmunes2 == false ){
            tejido.introducir_linfocitos_aleatorios(pg->cantidad_de_linfocitos);
            introduje_celulas_inmunes2 = true;
        }

        // Diffusion simulation step
        tejido.microambiente.simular_difusion_decaimiento( c->dt_difusion );

        // Gradient calculation at mechanical time steps
        double T_mec = pg->tiempo_total - n;
        if( pg->calcular_gradientes && fabs(T_mec - c->dt_mecanica) < tolerancia_de_la_mecanica){
                tejido.microambiente.calcular_todos_los_vectores_de_gradientes();
                n = pg->tiempo_total;
        }

        // Cell updates (mechanics and biology)
        tejido.cdc.actualizar_todas_las_celulas(pg->tiempo_total, c->dt_difusion, c->dt_mecanica, c->dt_ciclo);

        // Output data at specific intervals
        double T_esc = pg->tiempo_total - m;
        double T_esc2 = pg->tiempo_total - m2;
        //if(fabs(T_esc2 - 60) < tolerancia_de_la_escritura || pg->tiempo_total < 0.01){
            //std::cout << pg->tiempo_total << "\n";
            //m2 = pg->tiempo_total;
        //}
        if( fabs(T_esc - tiempo_de_escritura) < tolerancia_de_la_escritura || pg->tiempo_total < 0.01){
                // Calculate tumor geometry metrics
                tejido.geometria_del_tumor();
                
                // Count dead cancer cells
                int cancer_muerto = 0;
                // Count alive CAR-T cells
                int alive_cart = 0;
                for( unsigned int i=0; i < todas_las_celulas.size(); i++ ){
                    if(todas_las_celulas[i]->tipo == 0 && todas_las_celulas[i]->fenotipo.muerte.muerta){
                        cancer_muerto = cancer_muerto + 1;
                    }
                    
                    if(todas_las_celulas[i]->tipo == 2 && // Linfocites
                    !todas_las_celulas[i]->fenotipo.muerte.muerta) { // Alive
                        alive_cart++;
                    }

                }


                // Write summary data to output file
                outfile << std::setprecision(12) << pg->tiempo_total << " " << tejido.volumen_del_tumor << " " << tejido.volumen_del_tumor2 << " " << tejido.radio_del_tumor << " " << tejido.celulas_tumorales << " " << cancer_muerto << " " << todas_las_celulas.size() << " " << alive_cart << "\n";
                outfile.flush();
                m = pg->tiempo_total;
                
                // Write detailed cell state data to XYZ file
                std::ofstream outfile2 (results_dir + "/Datos_" + std::to_string( (int) pg->tiempo_total ) + ".xyz" );
                outfile2 << todas_las_celulas.size() << "\n";
                outfile2 << "\n";
                int fallecida;
                string adh;
                int adhint;
                string nombreciclo;

                // Output each cell's state and properties
                for( unsigned int i=0; i < todas_las_celulas.size(); i++ ){
                    int onco = 5;
                    fallecida = todas_las_celulas[i]->tipo;
                    
                    // Determine cell state based on death status and oncoprotein levels
                    if(todas_las_celulas[i]->fenotipo.muerte.muerta){
                        fallecida = 9;
                        onco = 0;
                    }else if(todas_las_celulas[i]->fenotipo.secrecion.oncoproteina >=1.5 && todas_las_celulas[i]->tipo != 2){
                        onco = 1;
                    }else if(todas_las_celulas[i]->fenotipo.secrecion.oncoproteina >=1.0 && todas_las_celulas[i]->fenotipo.secrecion.oncoproteina < 1.5 && todas_las_celulas[i]->tipo != 2){
                        onco = 2;
                    }else if(todas_las_celulas[i]->fenotipo.secrecion.oncoproteina >=0.5 && todas_las_celulas[i]->fenotipo.secrecion.oncoproteina < 1.0 && todas_las_celulas[i]->tipo != 2){
                        onco = 3;
                    }else if (todas_las_celulas[i]->fenotipo.secrecion.oncoproteina >=0.0 && todas_las_celulas[i]->fenotipo.secrecion.oncoproteina < 0.5 && todas_las_celulas[i]->tipo != 2){
                        onco = 4;
                    }
                    
                    nombreciclo = todas_las_celulas[i]->fenotipo.ciclo.pCiclo_Modelo->nombre;
                    
                    // Check for cell adhesion
                    if(todas_las_celulas[i]->adherida){
                        adh = "adherida";
                        adhint = todas_las_celulas[i]->celula_adherida->id;
                    }else{
                        adh = "no";
                        adhint = 0;
                    }
                    
                    // Write cell data to file
                    outfile2 << todas_las_celulas[i]->id << " " << todas_las_celulas[i]->voxel << " " << todas_las_celulas[i]->posicion.x << " " << todas_las_celulas[i]->posicion.y << " " << todas_las_celulas[i]->posicion.z << " " << todas_las_celulas[i]->fenotipo.geometria.radio << " " << todas_las_celulas[i]->voxel_del_microambiente << " " << todas_las_celulas[i]->vector_de_densidades_mas_cercano()[0] << " " << todas_las_celulas[i]->fenotipo.secrecion.oncoproteina << " " << todas_las_celulas[i]->madre << " " << todas_las_celulas[i]->tipo << " " << fallecida << " " << onco << " " << todas_las_celulas[i]->vector_de_densidades_mas_cercano()[1]<< " " << nombreciclo << " " << adh << " " << adhint << "\n";
                }
                
                // The VTK file generation code is commented out but would write visualization files
               
                     std::ofstream outfile3 (results_dir + "/PM1_" + std::to_string( (int) pg->tiempo_total ) + ".vtu" );
                    outfile3 << "<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'> \n";
                    outfile3 << "\t<UnstructuredGrid>\n";
                    outfile3 << "\t\t<Piece NumberOfPoints='" << todas_las_celulas.size() << "' NumberOfCells='0'>\n";
                    outfile3 << "\t\t\t<Points>\n";
                    outfile3 << "\t\t\t\t<DataArray name='Position' type='Float32' NumberOfComponents='3' format='ascii'>\n";
                    for(long unsigned int j=0; j<todas_las_celulas.size(); j++){
                        outfile3 << "\t\t\t\t" << todas_las_celulas[j]->posicion.x << " " << todas_las_celulas[j]->posicion.y << " " << todas_las_celulas[j]->posicion.z << "\n";
                    }
                    outfile3 << "\t\t\t\t</DataArray>\n";
                    outfile3 << "\t\t\t</Points>\n";
                    outfile3 << "\t\t\t<PointData  Vectors='vector'>\n";
                    outfile3 << "\t\t\t\t<DataArray type='Float32' Name='Radio' format='ascii'>\n";
                    for(long unsigned int j=0; j<todas_las_celulas.size(); j++){
                        outfile3 << "\t\t\t\t" << todas_las_celulas[j]->fenotipo.geometria.radio << "\n";
                    }
                    outfile3 << "\t\t\t\t</DataArray>\n";
                    outfile3 << "\t\t\t\t<DataArray type='Float32' Name='Oncoproteina' format='ascii'>\n";
                    for(long unsigned int j=0; j<todas_las_celulas.size(); j++){
                        outfile3 << "\t\t\t\t" << todas_las_celulas[j]->fenotipo.secrecion.oncoproteina << "\n";
                    }
                    outfile3 << "\t\t\t\t</DataArray>\n";
                    outfile3 << "\t\t\t\t<DataArray type='Float32' Name='Tipo' format='ascii'>\n";

                    for(long unsigned int j=0; j<todas_las_celulas.size(); j++){
                    int tipo2 = todas_las_celulas[j]->tipo;
                    if(todas_las_celulas[j]->fenotipo.muerte.muerta == true){
                        tipo2=1;
                    }
                        outfile3 << "\t\t\t\t" << tipo2 << "\n";
                    }
                    outfile3 << "\t\t\t\t</DataArray>\n";
                    outfile3 << "\t\t\t</PointData>\n";
                    outfile3 << "\t\t\t<Cells>\n";
                    outfile3 << "\t\t\t\t<DataArray type='Int32' Name='connectivity' format='ascii'>\n";
                    outfile3 << "\t\t\t\t</DataArray>\n";
                    outfile3 << "\t\t\t\t<DataArray type='Int32' Name='offsets' format='ascii'>\n";
                    outfile3 << "\t\t\t\t</DataArray>\n";
                    outfile3 << "\t\t\t\t<DataArray type='UInt8' Name='types' format='ascii'>\n";
                    outfile3 << "\t\t\t\t</DataArray>\n";
                    outfile3 << "\t\t\t</Cells>\n";
                    outfile3 << "\t\t</Piece>\n";
                    outfile3 << "\t</UnstructuredGrid>\n";
                    outfile3 << "</VTKFile>";



                    std::ofstream outfile4(results_dir + "/HM1_" + std::to_string( (int) pg->tiempo_total ) + ".vtu" );
                    outfile4 << "<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'> \n";
                    outfile4 << "\t<UnstructuredGrid>\n";
                    outfile4 << "\t\t<Piece NumberOfPoints='" << tejido.microambiente.mgrilla.voxeles.size() << "' NumberOfCells='0'>\n";
                    outfile4 << "\t\t\t<Points>\n";
                    outfile4 << "\t\t\t\t<DataArray name='Position' type='Float32' NumberOfComponents='3' format='ascii'>\n";
                    for(long unsigned int j=0; j<tejido.microambiente.mgrilla.voxeles.size(); j++){
                        outfile4 << "\t\t\t\t" << tejido.microambiente.mgrilla.voxeles[j].centro.x << " " << tejido.microambiente.mgrilla.voxeles[j].centro.y << " " << tejido.microambiente.mgrilla.voxeles[j].centro.z << "\n";
                    }
                    outfile4 << "\t\t\t\t</DataArray>\n";
                    outfile4 << "\t\t\t</Points>\n";
                    outfile4 << "\t\t\t<PointData  Vectors='vector'>\n";
                    outfile4 << "\t\t\t\t<DataArray type='Float32' Name='Radio' format='ascii'>\n";
                    for(long unsigned int j=0; j<tejido.microambiente.mgrilla.voxeles.size(); j++){
                        outfile4 << "\t\t\t\t" << pg->m_dx << "\n";
                    }
                    outfile4 << "\t\t\t\t</DataArray>\n";
                    outfile4 << "\t\t\t\t<DataArray type='Float32' Name='Oxigeno' format='ascii'>\n";
                    for(long unsigned int j=0; j<tejido.microambiente.mgrilla.voxeles.size(); j++){
                        outfile4 << "\t\t\t\t" << tejido.microambiente.vector_de_densidades(j)[0] << "\n";
                    }
                    outfile4 << "\t\t\t\t</DataArray>\n";
                    outfile4 << "\t\t\t\t<DataArray type='Float32' Name='Oncoproteina' format='ascii'>\n";
                    for(long unsigned int j=0; j<tejido.microambiente.mgrilla.voxeles.size(); j++){
                        outfile4 << "\t\t\t\t" << tejido.microambiente.vector_de_densidades(j)[1] << "\n";
                    }
                    outfile4 << "\t\t\t\t</DataArray>\n";
                    outfile4 << "\t\t\t</PointData>\n";
                    outfile4 << "\t\t\t<Cells>\n";
                    outfile4 << "\t\t\t\t<DataArray type='Int32' Name='connectivity' format='ascii'>\n";
                    outfile4 << "\t\t\t\t</DataArray>\n";
                    outfile4 << "\t\t\t\t<DataArray type='Int32' Name='offsets' format='ascii'>\n";
                    outfile4 << "\t\t\t\t</DataArray>\n";
                    outfile4 << "\t\t\t\t<DataArray type='UInt8' Name='types' format='ascii'>\n";
                    outfile4 << "\t\t\t\t</DataArray>\n";
                    outfile4 << "\t\t\t</Cells>\n";
                    outfile4 << "\t\t</Piece>\n";
                    outfile4 << "\t</UnstructuredGrid>\n";
                    outfile4 << "</VTKFile>";

  
  
                
/////////////////////////////////////////////////////////////////////////                
                
                
        }

        // Periodically update cell voxel assignments for spatial organization
        if(fmod(pg->tiempo_total, c->dt_mecanica*20)<0.1){
            tejido.cdc.actualizar_voxeles_de_celulas();
        }
	}


      //Debug
  std::cout << "Final accumulated probability: " 
          << std::setprecision(30) << Linfocito::acumulator_probabilities 
          << std::endl;
  // Write the final value to a file
  std::ofstream file("out/final_accumulated_probability.txt");
  if (file.is_open()) {
    file << std::setprecision(30) << Linfocito::acumulator_probabilities;
    file.close();
  }

	return 0;

}




/*
Para debuggear el código
#define INFO(msg) \
    fprintf(stderr, "info: %s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, "%s\n", msg);
*/
