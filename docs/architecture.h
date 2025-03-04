/**
 * @page architecture Architecture
 * 
 * @section arch_overview Architecture Overview
 * 
 * <div class="architecture-overview">
 * <p>
 * The Cancer Cell Agent Simulation is built on a modular, object-oriented architecture that enables
 * detailed modeling of cancer cell behavior in a 3D microenvironment. The architecture follows these
 * key design principles:
 * </p>
 * 
 * <ul>
 * <li><b>Separation of Concerns</b>: Different aspects of the simulation (cells, microenvironment, mechanics, etc.) are handled by specialized components.</li>
 * <li><b>Extensibility</b>: The system is designed to be extended with new cell types, behaviors, and environmental factors.</li>
 * <li><b>Performance Optimization</b>: Critical operations are optimized for computational efficiency in large-scale simulations.</li>
 * <li><b>Research Applicability</b>: The architecture prioritizes features relevant to cancer research.</li>
 * </ul>
 * </div>
 * 
 * @section arch_components Core Architectural Components
 * 
 * <div class="arch-components">
 * <p>
 * The simulation architecture consists of several interconnected components:
 * </p>
 * 
 * <h3>Agent System</h3>
 * <p>
 * The agent system is responsible for:
 * </p>
 * <ul>
 * <li>Representing individual cells as autonomous agents</li>
 * <li>Managing agent properties, states, and behaviors</li>
 * <li>Handling agent-agent and agent-environment interactions</li>
 * <li>Implementing specialized cell types (tumor cells, immune cells, etc.)</li>
 * </ul>
 * 
 * <h3>Microenvironment System</h3>
 * <p>
 * The microenvironment system handles:
 * </p>
 * <ul>
 * <li>3D spatial discretization using the finite volume method</li>
 * <li>Diffusion of multiple chemical substrates</li>
 * <li>Boundary conditions and source/sink terms</li>
 * <li>Gradient calculations for directional cell behaviors</li>
 * </ul>
 * 
 * <h3>Mechanics System</h3>
 * <p>
 * The mechanics system computes:
 * </p>
 * <ul>
 * <li>Cell-cell adhesion and repulsion forces</li>
 * <li>Cell-ECM (extracellular matrix) interactions</li>
 * <li>Cell motility and migration</li>
 * <li>Physical constraints and mechanical feedback</li>
 * </ul>
 * 
 * <h3>Cell Biology System</h3>
 * <p>
 * The cell biology system models:
 * </p>
 * <ul>
 * <li>Cell cycle progression and division</li>
 * <li>Cell death (apoptosis, necrosis)</li>
 * <li>Metabolic processes (consumption, secretion)</li>
 * <li>Phenotype transitions based on microenvironmental cues</li>
 * </ul>
 * 
 * <h3>Tissue Organization System</h3>
 * <p>
 * The tissue system coordinates:
 * </p>
 * <ul>
 * <li>Overall tissue structure and geometry</li>
 * <li>Spatial organization of cells</li>
 * <li>Tissue-level properties and measurements</li>
 * <li>Multi-scale integration between cells and tissue</li>
 * </ul>
 * 
 * <h3>Simulation Core</h3>
 * <p>
 * The simulation core manages:
 * </p>
 * <ul>
 * <li>Simulation initialization and configuration</li>
 * <li>The main simulation loop and time stepping</li>
 * <li>Process scheduling and synchronization</li>
 * <li>Data collection, output, and visualization</li>
 * </ul>
 * </div>
 * 
 * @section arch_workflow Architectural Workflow
 * 
 * <div class="arch-workflow">
 * <p>
 * The architectural components interact in a defined workflow:
 * </p>
 * 
 * <ol>
 * <li><b>Initialization Phase</b>
 *   <ul>
 *     <li>Load simulation parameters</li>
 *     <li>Initialize microenvironment grid</li>
 *     <li>Create initial cell population</li>
 *     <li>Setup data structures and output files</li>
 *   </ul>
 * </li>
 * 
 * <li><b>Simulation Loop</b> - For each time step:
 *   <ul>
 *     <li><b>Microenvironment Update</b>: Solve diffusion equations for all substrates</li>
 *     <li><b>Cell Sensing</b>: Cells probe their local microenvironment</li>
 *     <li><b>Cell Decisions</b>: Cells update their internal state and make decisions</li>
 *     <li><b>Mechanics</b>: Calculate forces and update cell positions</li>
 *     <li><b>Cell Processes</b>: Execute cell division, death, and other processes</li>
 *     <li><b>Data Collection</b>: Record system state at specified intervals</li>
 *   </ul>
 * </li>
 * 
 * <li><b>Finalization Phase</b>
 *   <ul>
 *     <li>Generate final outputs and statistics</li>
 *     <li>Close files and clean up resources</li>
 *     <li>Prepare visualization data</li>
 *   </ul>
 * </li>
 * </ol>
 * </div>
 * 
 * @section arch_design Architectural Design Decisions
 * 
 * <div class="arch-design">
 * <p>
 * Key architectural design decisions include:
 * </p>
 * 
 * <ul>
 * <li><b>Off-Lattice Approach</b>: Cells have continuous positions rather than being confined to grid points, allowing for more realistic physical interactions.</li>
 * <li><b>Hybrid Discrete-Continuous Model</b>: Cells are discrete entities while the microenvironment is modeled as continuous fields, combining the advantages of both approaches.</li>
 * <li><b>Modular Component Design</b>: Components can be modified or replaced independently, facilitating model experimentation and extension.</li>
 * <li><b>Process-Based Simulation Loop</b>: The simulation alternates between different processes (diffusion, mechanics, biology) to maintain physical and biological realism.</li>
 * <li><b>Hierarchical Organization</b>: The architecture spans multiple scales, from subcellular processes to tissue-level organization.</li>
 * </ul>
 * </div>
 * 
 * @section arch_extensions Architectural Extensions for CAR T-cell Therapy
 * 
 * <div class="arch-extensions">
 * <p>
 * To support CAR T-cell therapy simulations, the architecture includes these specialized extensions:
 * </p>
 * 
 * <ul>
 * <li><b>Immune Cell Agents</b>: Enhanced agent models for CAR T-cells with specialized behaviors</li>
 * <li><b>Recognition Dynamics</b>: Mechanisms for CAR T-cells to identify and target cancer cells expressing specific antigens</li>
 * <li><b>Cytotoxic Interactions</b>: Models for immune-mediated killing of cancer cells</li>
 * <li><b>T-cell Exhaustion</b>: Representation of diminishing T-cell efficacy over multiple cancer cell engagements</li>
 * <li><b>Heterogeneous Response</b>: Framework for modeling varied responses among cancer cell subpopulations</li>
 * </ul>
 * </div>
 */ 