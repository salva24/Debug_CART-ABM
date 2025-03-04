/**
 * @page usage_guide Usage Guide
 * 
 * @section usage_overview Overview
 * 
 * <div class="usage-overview">
 * <p>
 * This guide provides detailed instructions for setting up, configuring, running, and analyzing 
 * simulations using the Cancer Cell Agent Simulation framework, with a particular focus on CAR T-cell 
 * therapy simulations as described in the associated <i>Scientific Reports</i> publication.
 * </p>
 * </div>
 * 
 * @section system_requirements System Requirements
 * 
 * <div class="system-requirements">
 * <p>
 * To run the simulation, your system should meet the following requirements:
 * </p>
 * 
 * <h3>Hardware Requirements</h3>
 * <ul>
 * <li><b>Processor</b>: Multi-core CPU recommended (e.g., Intel i5/i7 or equivalent)</li>
 * <li><b>Memory</b>: Minimum 8GB RAM, 16GB+ recommended for larger simulations</li>
 * <li><b>Storage</b>: At least 1GB free space for the software and simulation outputs</li>
 * </ul>
 * 
 * <h3>Software Requirements</h3>
 * <ul>
 * <li><b>Operating System</b>: Linux (preferred), macOS, or Windows</li>
 * <li><b>Compiler</b>: GCC 7.0+ or compatible C++ compiler</li>
 * <li><b>Libraries</b>: Standard C++ libraries</li>
 * <li><b>Visualization Tools</b>: ParaView or VisIt for viewing VTK output files (optional)</li>
 * </ul>
 * </div>
 * 
 * @section installation Installation
 * 
 * <div class="installation">
 * <p>
 * To install the simulation framework:
 * </p>
 * 
 * <ol>
 * <li><b>Clone the Repository</b>
 * <pre class="fragment">
 * git clone [repository-url]
 * cd [repository-directory]
 * </pre>
 * </li>
 * 
 * <li><b>Build the Project</b>
 * <pre class="fragment">
 * make -f MakefileRo
 * </pre>
 * The executable will be created in the <code>build/</code> directory as <code>ciclo</code>.
 * </li>
 * 
 * <li><b>Verify Installation</b>
 * <pre class="fragment">
 * ./build/ciclo --version
 * </pre>
 * This should display the version information.
 * </li>
 * </ol>
 * </div>
 * 
 * @section configuration Configuration
 * 
 * <div class="configuration">
 * <p>
 * The simulation is configured through parameter files. These files specify all aspects of the simulation,
 * from environment setup to cell properties and simulation duration.
 * </p>
 * 
 * <h3>Parameter File Structure</h3>
 * <p>
 * Parameter files are structured in sections with key-value pairs:
 * </p>
 * 
 * <pre class="fragment">
 * [Section1]
 * parameter1 = value1
 * parameter2 = value2
 * 
 * [Section2]
 * parameter3 = value3
 * ...
 * </pre>
 * 
 * <h3>Key Configuration Sections</h3>
 * <p>
 * The main configuration sections include:
 * </p>
 * 
 * <ul>
 * <li><b>Simulation</b>: General simulation parameters (duration, time step, output frequency)</li>
 * <li><b>Microenvironment</b>: Domain size, substrates, boundary conditions</li>
 * <li><b>Cell Properties</b>: Cell types, phenotypes, behaviors</li>
 * <li><b>Immune System</b>: CAR T-cell properties, introduction timing, counts</li>
 * <li><b>Tumor</b>: Initial tumor configuration, heterogeneity parameters</li>
 * <li><b>Output</b>: Specification of output files and formats</li>
 * </ul>
 * 
 * <h3>Example Parameter File</h3>
 * <p>
 * A minimal example parameter file for a CAR T-cell therapy simulation:
 * </p>
 * 
 * <pre class="fragment">
 * [Simulation]
 * duration = 48.0     # Simulation duration in hours
 * time_step = 0.1     # Time step in hours
 * output_interval = 1.0  # Output frequency in hours
 * 
 * [Microenvironment]
 * domain_size = 1000.0, 1000.0, 1000.0  # Domain dimensions in µm
 * voxel_size = 20.0    # Voxel size in µm
 * substrates = oxygen, glucose, waste   # Diffusing substrates
 * 
 * [Tumor]
 * initial_radius = 200.0  # Initial tumor radius in µm
 * cell_count = 1000       # Initial tumor cell count
 * heterogeneity = true    # Enable tumor heterogeneity
 * antigen_expression = 0.8, 0.5, 0.2  # Antigen expression levels (fraction of cells)
 * 
 * [CAR_T_Cells]
 * introduction_time = 24.0  # When to introduce CAR T-cells (hours)
 * count = 5000              # Number of CAR T-cells to introduce
 * binding_affinity = 0.85   # CAR-antigen binding strength (0-1)
 * killing_efficiency = 0.75 # Probability of inducing apoptosis (0-1)
 * exhaustion_rate = 0.1     # Rate of T-cell exhaustion per engagement
 * </pre>
 * </div>
 * 
 * @section running Running Simulations
 * 
 * <div class="running">
 * <p>
 * To run a simulation:
 * </p>
 * 
 * <ol>
 * <li><b>Basic Execution</b>
 * <pre class="fragment">
 * ./build/ciclo params.ini
 * </pre>
 * This runs the simulation using parameters from <code>params.ini</code>.
 * </li>
 * 
 * <li><b>Output Directory Specification</b>
 * <pre class="fragment">
 * ./build/ciclo params.ini --output=my_results
 * </pre>
 * This directs output to the <code>my_results</code> directory.
 * </li>
 * 
 * <li><b>Resume Simulation</b>
 * <pre class="fragment">
 * ./build/ciclo params.ini --resume=checkpoint_file.dat
 * </pre>
 * This resumes a simulation from a checkpoint file.
 * </li>
 * </ol>
 * 
 * <h3>Runtime Information</h3>
 * <p>
 * During execution, the simulation provides progress updates:
 * </p>
 * <ul>
 * <li>Current simulation time</li>
 * <li>Cell counts (tumor cells, CAR T-cells)</li>
 * <li>Estimated completion time</li>
 * <li>Output file generation notifications</li>
 * </ul>
 * </div>
 * 
 * @section output Output and Analysis
 * 
 * <div class="output-analysis">
 * <p>
 * The simulation produces several types of output files for analysis:
 * </p>
 * 
 * <h3>Summary Data</h3>
 * <p>
 * <code>DatosFinales.dat</code> contains time series data including:
 * </p>
 * <ul>
 * <li>Time points</li>
 * <li>Total cell counts by type</li>
 * <li>Tumor volume</li>
 * <li>Tumor radius</li>
 * <li>Treatment efficacy metrics</li>
 * </ul>
 * 
 * <h3>Detailed Cell Data</h3>
 * <p>
 * <code>Datos_[time].xyz</code> files contain detailed cell-level information at specific time points:
 * </p>
 * <ul>
 * <li>Cell positions</li>
 * <li>Cell types and states</li>
 * <li>Cell-specific properties (volume, age, etc.)</li>
 * <li>Antigen expression levels</li>
 * <li>T-cell activation status</li>
 * </ul>
 * 
 * <h3>Visualization Files</h3>
 * <p>
 * VTK files for visualization in ParaView or VisIt:
 * </p>
 * <ul>
 * <li><code>PM1_[time].xml</code>: Cell-based visualization data</li>
 * <li><code>HM1_[time].xml</code>: Microenvironment visualization data</li>
 * </ul>
 * 
 * <h3>Analysis Tools</h3>
 * <p>
 * Common analysis approaches include:
 * </p>
 * <ul>
 * <li><b>Time Series Analysis</b>: Plotting tumor volume and cell counts over time</li>
 * <li><b>Spatial Analysis</b>: Examining cell distributions and clustering</li>
 * <li><b>Efficacy Assessment</b>: Calculating tumor reduction and treatment response metrics</li>
 * <li><b>Heterogeneity Analysis</b>: Quantifying changes in subpopulation distributions</li>
 * <li><b>3D Visualization</b>: Creating renders and animations of the simulation</li>
 * </ul>
 * </div>
 * 
 * @section case_studies Case Studies
 * 
 * <div class="case-studies">
 * <p>
 * These case studies demonstrate how to set up and analyze specific scenarios:
 * </p>
 * 
 * <h3>Case Study 1: Homogeneous vs. Heterogeneous Tumors</h3>
 * <p>
 * This case study compares CAR T-cell efficacy against:
 * </p>
 * <ul>
 * <li>A homogeneous tumor with uniform antigen expression</li>
 * <li>A heterogeneous tumor with varied antigen expression</li>
 * </ul>
 * <p>
 * Parameter variations focus on:
 * </p>
 * <pre class="fragment">
 * [Tumor]
 * # Homogeneous case
 * heterogeneity = false
 * antigen_expression = 0.9
 * 
 * # Heterogeneous case
 * heterogeneity = true
 * antigen_expression = 0.9, 0.5, 0.1
 * subpopulation_fractions = 0.3, 0.4, 0.3
 * </pre>
 * 
 * <h3>Case Study 2: T-cell Dosing Optimization</h3>
 * <p>
 * This study explores different CAR T-cell dosing strategies:
 * </p>
 * <ul>
 * <li>Single high dose</li>
 * <li>Multiple smaller doses</li>
 * <li>Escalating dose schedule</li>
 * </ul>
 * <p>
 * Parameter variations include:
 * </p>
 * <pre class="fragment">
 * [CAR_T_Cells]
 * # Single dose
 * introduction_time = 24.0
 * count = 10000
 * 
 * # Multiple doses
 * introduction_time = 24.0, 48.0, 72.0
 * count = 3000, 3000, 3000
 * </pre>
 * </div>
 * 
 * @section troubleshooting Troubleshooting
 * 
 * <div class="troubleshooting">
 * <p>
 * Common issues and their solutions:
 * </p>
 * 
 * <h3>Performance Issues</h3>
 * <ul>
 * <li><b>Problem</b>: Simulation running too slowly
 * <br><b>Solution</b>: Reduce domain size or increase voxel size to reduce computational load</li>
 * <li><b>Problem</b>: Memory usage too high
 * <br><b>Solution</b>: Reduce cell count or limit output frequency</li>
 * </ul>
 * 
 * <h3>Numerical Stability</h3>
 * <ul>
 * <li><b>Problem</b>: Diffusion solver instability
 * <br><b>Solution</b>: Decrease time step or increase voxel size</li>
 * <li><b>Problem</b>: Cell movement instability
 * <br><b>Solution</b>: Reduce cell movement speed or adhesion strength</li>
 * </ul>
 * 
 * <h3>Output Issues</h3>
 * <ul>
 * <li><b>Problem</b>: Missing output files
 * <br><b>Solution</b>: Check write permissions and disk space</li>
 * <li><b>Problem</b>: Corrupted visualization files
 * <br><b>Solution</b>: Verify VTK file format compatibility with visualization software</li>
 * </ul>
 * </div>
 */ 