/**
 * @page implementation Implementation Details
 * 
 * @section implementation_overview Overview
 * 
 * <div class="implementation-overview">
 * <p>
 * This page documents the implementation details of the Cancer Cell Agent Simulation, including class hierarchies,
 * data structures, algorithms, and the technical approach to modeling complex biological processes such as
 * CAR T-cell interactions with tumor organoids.
 * </p>
 * </div>
 * 
 * @section class_hierarchy Class Hierarchy
 * 
 * <div class="class-hierarchy">
 * <p>
 * The simulation is structured around the following class hierarchy:
 * </p>
 * 
 * <h3>Core Classes</h3>
 * <ul>
 * <li><b>Celula</b>: Base class for all cell agents
 *   <ul>
 *     <li><b>Linfocito</b>: Specialized class for immune cells, including CAR T-cells</li>
 *     <li><b>Celula_Cancerosa</b>: Specialized class for cancer cells with heterogeneous properties</li>
 *   </ul>
 * </li>
 * <li><b>Microambiente</b>: Manages the 3D diffusion of substrates and chemical signals</li>
 * <li><b>Tejido</b>: Coordinates tissue-level organization and cell-microenvironment interactions</li>
 * <li><b>Contenedor_de_Celulas</b>: Handles efficient storage and retrieval of cells</li>
 * <li><b>Ciclo_Modelo</b>: Implements the cell cycle model with phase transitions</li>
 * </ul>
 * 
 * <h3>Supporting Classes</h3>
 * <ul>
 * <li><b>Vector</b>: Utility class for 3D vector operations</li>
 * <li><b>Voxel</b>: Represents a discrete volume element in the microenvironment</li>
 * <li><b>Muerte</b>: Handles cell death processes (apoptosis, necrosis)</li>
 * <li><b>Volumen</b>: Manages cell volume changes</li>
 * <li><b>Secrecion</b>: Controls cellular secretion of substrates</li>
 * <li><b>Random</b>: Provides random number generation functionality</li>
 * <li><b>Parametros_globales</b>: Stores global simulation parameters</li>
 * </ul>
 * </div>
 * 
 * @section celula_impl Cell Implementation
 * 
 * <div class="celula-implementation">
 * <p>
 * The Celula class is the foundation of the agent-based approach:
 * </p>
 * 
 * <h3>Key Properties</h3>
 * <ul>
 * <li><b>Spatial Properties</b>: Position, velocity, and orientation in 3D space</li>
 * <li><b>Physical Properties</b>: Volume, radius, surface area</li>
 * <li><b>Mechanical Properties</b>: Adhesion coefficients, repulsion parameters</li>
 * <li><b>Biological Properties</b>: Cell cycle state, phenotype, age</li>
 * <li><b>Microenvironment Interaction</b>: Consumption rates, secretion rates</li>
 * </ul>
 * 
 * <h3>Key Methods</h3>
 * <ul>
 * <li><b>Mechanics Functions</b>: Calculate and apply forces between cells and boundaries</li>
 * <li><b>Migration Functions</b>: Handle cell movement in response to cues</li>
 * <li><b>Division Functions</b>: Manage cell division process and daughter cell creation</li>
 * <li><b>Death Functions</b>: Handle different modes of cell death</li>
 * <li><b>Sensing Functions</b>: Sample and respond to the local microenvironment</li>
 * </ul>
 * 
 * <h3>Cell Specialization</h3>
 * <p>
 * The implementation includes specialized cell types:
 * </p>
 * <ul>
 * <li><b>Linfocito Class</b>: Extends Celula with immune cell behaviors
 *   <ul>
 *     <li><code>chequear_vecinos_para_adherirse</code>: Scans for potential target cells</li>
 *     <li><code>adherirse_a_objetivo</code>: Forms immune synapse with target cell</li>
 *     <li><code>desencadenar_apoptosis</code>: Initiates apoptosis in target cell</li>
 *     <li><code>actualizar_estado_agotamiento</code>: Updates T-cell exhaustion level</li>
 *   </ul>
 * </li>
 * <li><b>Celula_Cancerosa Class</b>: Extends Celula with cancer-specific behaviors
 *   <ul>
 *     <li><code>set_nivel_expresion_antigeno</code>: Sets antigen expression level</li>
 *     <li><code>regular_antigeno</code>: Adjusts antigen expression (e.g., in response to stress)</li>
 *     <li><code>responder_a_senales</code>: Modifies behavior based on microenvironment signals</li>
 *   </ul>
 * </li>
 * </ul>
 * </div>
 * 
 * @section microambiente_impl Microenvironment Implementation
 * 
 * <div class="microambiente-implementation">
 * <p>
 * The Microambiente class implements a finite volume approach to substrate diffusion:
 * </p>
 * 
 * <h3>Implementation Details</h3>
 * <ul>
 * <li><b>Spatial Discretization</b>: The domain is divided into voxels (Voxel objects)</li>
 * <li><b>Diffusion Solver</b>: Implements numerical solutions to the diffusion equation
 *   <ul>
 *     <li>Forward Euler method for time discretization</li>
 *     <li>Central differences for spatial discretization</li>
 *     <li>Stability-enforced time step constraints</li>
 *   </ul>
 * </li>
 * <li><b>Boundary Conditions</b>: Supports Dirichlet, Neumann, and periodic boundaries</li>
 * <li><b>Source/Sink Terms</b>: Accounts for cell-based consumption and secretion</li>
 * </ul>
 * 
 * <h3>Key Methods</h3>
 * <ul>
 * <li><code>inicializar_densidades</code>: Sets up initial substrate distributions</li>
 * <li><code>difundir_densidades</code>: Solves the diffusion equations for one time step</li>
 * <li><code>actualizar_nodo_de_dirichlet</code>: Applies Dirichlet boundary conditions</li>
 * <li><code>obtener_gradiente</code>: Calculates concentration gradients for cell sensing</li>
 * <li><code>registrar_secretor</code>: Registers a cell as a source/sink for substrates</li>
 * </ul>
 * </div>
 * 
 * @section tejido_impl Tissue Implementation
 * 
 * <div class="tejido-implementation">
 * <p>
 * The Tejido class coordinates the overall tissue organization:
 * </p>
 * 
 * <h3>Implementation Details</h3>
 * <ul>
 * <li><b>Initialization</b>: Sets up initial cell populations and microenvironment</li>
 * <li><b>Integration</b>: Connects cells with the microenvironment</li>
 * <li><b>Organization</b>: Manages tissue structure and geometry</li>
 * <li><b>Measurement</b>: Computes tissue-level metrics</li>
 * </ul>
 * 
 * <h3>Key Methods</h3>
 * <ul>
 * <li><code>inicializar_tejido</code>: Creates the initial tissue configuration</li>
 * <li><code>crear_esfera_celulas</code>: Creates a spherical arrangement of cells</li>
 * <li><code>introducir_linfocitos</code>: Adds T-cells at specified times and locations</li>
 * <li><code>actualizar_tejido</code>: Updates the entire tissue for one time step</li>
 * <li><code>calcular_metricas</code>: Computes various tissue metrics (volume, cell counts, etc.)</li>
 * </ul>
 * </div>
 * 
 * @section car_tcell_impl CAR T-cell Interaction Implementation
 * 
 * <div class="car-tcell-implementation">
 * <p>
 * The CAR T-cell interaction with tumor cells is implemented through several coordinated mechanisms:
 * </p>
 * 
 * <h3>CAR T-cell Simulation Elements</h3>
 * <ul>
 * <li><b>Migration Algorithm</b>: T-cells move through the tissue using a biased random walk algorithm
 *   <ul>
 *     <li>Direction is influenced by chemokine gradients</li>
 *     <li>Persistence term maintains directional movement</li>
 *     <li>Physical barriers affect movement speed and direction</li>
 *   </ul>
 * </li>
 * <li><b>Target Recognition</b>: Implemented through neighbor scanning and probability-based binding
 *   <ul>
 *     <li>Binding probability proportional to target antigen expression</li>
 *     <li>Stochastic process for initial recognition</li>
 *     <li>Adhesion forces established upon successful binding</li>
 *   </ul>
 * </li>
 * <li><b>Killing Mechanism</b>: Multi-step process for target cell elimination
 *   <ul>
 *     <li>Engagement time tracking for minimal contact duration</li>
 *     <li>Probability-based killing trigger</li>
 *     <li>Apoptosis pathway activation in target cell</li>
 *     <li>T-cell detachment after successful killing</li>
 *   </ul>
 * </li>
 * <li><b>Exhaustion Model</b>: Represents diminishing T-cell functionality
 *   <ul>
 *     <li>Exhaustion state variable incremented with each killing</li>
 *     <li>Threshold-based function for efficiency reduction</li>
 *     <li>Eventually leads to T-cell apoptosis when fully exhausted</li>
 *   </ul>
 * </li>
 * </ul>
 * 
 * <h3>Implementation of Cancer Cell Heterogeneity</h3>
 * <ul>
 * <li><b>Antigen Expression Variation</b>: Implemented as cell-specific probability distribution</li>
 * <li><b>Susceptibility Differences</b>: Varied killing probability based on cell properties</li>
 * <li><b>Resistance Development</b>: Stochastic transitions to less susceptible states</li>
 * <li><b>Phenotype Adaptation</b>: Dynamic changes in antigen expression based on microenvironment</li>
 * </ul>
 * </div>
 * 
 * @section algorithms_impl Key Algorithms
 * 
 * <div class="algorithms-implementation">
 * <p>
 * The simulation employs several specialized algorithms:
 * </p>
 * 
 * <h3>Mechanical Interactions</h3>
 * <ul>
 * <li><b>Cell-Cell Force Calculation</b>: Modified Hertz contact model for repulsion, potential-based adhesion</li>
 * <li><b>Cell Movement</b>: Overdamped equation of motion (neglecting inertial terms)</li>
 * <li><b>Collision Detection</b>: Grid-based spatial binning for efficient neighbor finding</li>
 * </ul>
 * 
 * <h3>Division Algorithm</h3>
 * <ul>
 * <li><b>Division Axis Selection</b>: Based on local mechanical constraints and cell polarity</li>
 * <li><b>Volume Distribution</b>: Slightly asymmetric volume allocation to daughter cells</li>
 * <li><b>Property Inheritance</b>: Stochastic variation in inherited properties</li>
 * </ul>
 * 
 * <h3>Substrate Diffusion</h3>
 * <ul>
 * <li><b>Finite Volume Method</b>: Conservation-preserving discretization</li>
 * <li><b>Thomas Algorithm</b>: For tridiagonal system solving in each dimension</li>
 * <li><b>Operator Splitting</b>: For 3D diffusion with reaction terms</li>
 * </ul>
 * </div>
 * 
 * @section performance_impl Performance Optimizations
 * 
 * <div class="performance-implementation">
 * <p>
 * Several optimizations ensure computational efficiency:
 * </p>
 * 
 * <ul>
 * <li><b>Spatial Partitioning</b>: Grid-based partitioning for O(1) neighbor lookup</li>
 * <li><b>Memory Management</b>: Custom memory pool for frequent cell creation/destruction</li>
 * <li><b>Computation Scheduling</b>: Adaptive computation focusing on active regions</li>
 * <li><b>Lazy Evaluation</b>: Deferred computation of gradients and other expensive operations</li>
 * <li><b>Multi-threading</b>: Parallel processing of independent cell updates and diffusion steps</li>
 * </ul>
 * </div>
 */ 