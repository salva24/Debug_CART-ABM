/**
 * @mainpage Cancer Cell Agent Simulation
 * 
 * @section intro_sec Introduction
 * 
 * <div class="introduction">
 * <p>
 * This simulation is an off-lattice agent-based model designed to simulate the behavior of individual 
 * cells within a 3D tissue-like microenvironment. Each cell is an agent with its own properties 
 * (position, size, state, etc.) that interacts with a biochemical microenvironment containing 
 * various substrates and factors.
 * </p>
 * 
 * <p>
 * The model enables researchers to study tumor growth, cell interactions, immune responses, 
 * and the effects of microenvironmental factors on cancer progression, with a specific focus 
 * on studying heterogeneous tumor-derived organoid response to CAR T-cell therapy.
 * </p>
 * </div>
 * 
 * @section publication_sec Publication
 * 
 * <div class="publication">
 * <p>
 * This simulation framework has been published in <b>Scientific Reports</b> (Nature Portfolio):
 * </p>
 * 
 * <blockquote>
 * <p>
 * <b>In silico study of heterogeneous tumour-derived organoid response to CAR T-cell therapy</b><br>
 * <i>Scientific Reports</i>, volume 14, Article number: 11999 (2024)<br>
 * Published: May 29, 2024<br>
 * DOI: <a href="https://www.nature.com/articles/s41598-024-63125-5" target="_blank">https://www.nature.com/articles/s41598-024-63125-5</a>
 * </p>
 * </blockquote>
 * 
 * <p>
 * This publication describes the application of this agent-based model to study how heterogeneous 
 * tumor-derived organoids respond to CAR T-cell therapy, providing insights into the complex 
 * interactions between immune cells and tumor microenvironments.
 * </p>
 * 
 * <p>
 * The model builds upon previous research published in <b>PLOS Computational Biology</b>:
 * </p>
 * 
 * <blockquote>
 * <p>
 * <b>Physics-based tissue simulator to model multicellular systems: A study of liver regeneration and hepatocellular carcinoma recurrence</b><br>
 * <i>PLOS Computational Biology</i>, 19(3): e1010920 (2023)<br>
 * Published: March 6, 2023<br>
 * DOI: <a href="https://doi.org/10.1371/journal.pcbi.1010920" target="_blank">https://doi.org/10.1371/journal.pcbi.1010920</a>
 * </p>
 * </blockquote>
 * 
 * <p>
 * This earlier work introduced the core concepts of the agent-based model, focusing on tissue regeneration 
 * and tumor development in the context of liver regeneration after surgical hepatectomy.
 * </p>
 * </div>
 * 
 * @section arch_sec System Architecture
 * 
 * <div class="architecture">
 * <p>
 * The simulation is structured as an object-oriented system with the following key architectural elements:
 * </p>
 * 
 * <ul>
 * <li><b>Agent-Based Approach</b>: Individual cells are modeled as autonomous agents with their own properties, states, and behaviors.</li>
 * <li><b>Off-Lattice Design</b>: Cells can move continuously in 3D space rather than being confined to fixed grid positions.</li>
 * <li><b>Microenvironment</b>: A 3D grid system (finite volume method) models biochemical substrates like oxygen, nutrients, and signaling molecules that diffuse through the tissue.</li>
 * <li><b>Process-Based Simulation Loop</b>: The simulation alternates between three main processes:
 *   <ul>
 *     <li><b>Diffusion</b>: Updating biochemical concentrations throughout the microenvironment</li>
 *     <li><b>Mechanics</b>: Calculating physical forces and cell movements</li>
 *     <li><b>Biology</b>: Processing cell behaviors like division, death, and phenotype changes</li>
 *   </ul>
 * </li>
 * </ul>
 * </div>
 * 
 * @section car_tcell_sec CAR T-cell Therapy Simulation
 * 
 * <div class="car-tcell">
 * <p>
 * This model specifically simulates the interaction between CAR T-cells (Chimeric Antigen Receptor T-cells) 
 * and tumor-derived organoids. Key aspects of this simulation include:
 * </p>
 * 
 * <ul>
 * <li><b>T-cell Infiltration</b>: Modeling how CAR T-cells enter and navigate the tumor microenvironment</li>
 * <li><b>Target Recognition</b>: Simulating how CAR T-cells identify and bind to tumor cells</li>
 * <li><b>Cytotoxic Response</b>: Replicating the mechanisms by which CAR T-cells induce apoptosis in tumor cells</li>
 * <li><b>T-cell Exhaustion</b>: Accounting for the limited capacity of T-cells to repeatedly engage and kill tumor cells</li>
 * <li><b>Tumor Heterogeneity</b>: Modeling how varied tumor cell populations respond differently to CAR T-cell therapy</li>
 * </ul>
 * 
 * <p>
 * This sophisticated modeling approach provides researchers with a powerful tool to explore the complex 
 * dynamics of immunotherapy interventions in heterogeneous tumor environments.
 * </p>
 * </div>
 * 
 * @section core_sec Core Components
 * 
 * <div class="core-components">
 * <h3>Cell System (Celula)</h3>
 * <p>
 * The Celula class represents individual cells with:
 * </p>
 * <ul>
 * <li>Physical properties (position, velocity, size)</li>
 * <li>Phenotypic state (cycle phase, volume, death status)</li>
 * <li>Interaction with the microenvironment (consumption, secretion)</li>
 * <li>Mechanics (movement, adhesion)</li>
 * </ul>
 * <p>
 * Specialized cell types like lymphocytes (Linfocito) extend the base cell class with immune-specific behaviors.
 * </p>
 * <p>
 * See Celula.h and Celula.cpp for implementation details.
 * </p>
 * 
 * <h3>Microenvironment System (Microambiente)</h3>
 * <p>
 * The Microambiente class manages:
 * </p>
 * <ul>
 * <li>Diffusion of biochemical substrates through the tissue</li>
 * <li>Concentration gradients that influence cell behavior</li>
 * <li>Boundary conditions and substrate sources/sinks</li>
 * <li>Voxel-based spatial discretization for the finite volume method</li>
 * </ul>
 * <p>
 * See Microambiente.h and Microambiente.cpp for implementation details.
 * </p>
 * 
 * <h3>Tissue Organization (Tejido)</h3>
 * <p>
 * The Tejido class coordinates:
 * </p>
 * <ul>
 * <li>Overall tissue structure and initialization</li>
 * <li>Integration of cells with the microenvironment</li>
 * <li>Geometric arrangements of cells</li>
 * <li>Tissue-level properties and measurements</li>
 * </ul>
 * <p>
 * See Tejido.h and Tejido.cpp for implementation details.
 * </p>
 * 
 * <h3>Cell Container (Contenedor_de_Celulas)</h3>
 * <p>
 * The Contenedor_de_Celulas class handles:
 * </p>
 * <ul>
 * <li>Efficient cell organization and retrieval</li>
 * <li>Spatial indexing for neighbor interactions</li>
 * <li>Cell registration and management</li>
 * <li>Cell updates during simulation cycles</li>
 * </ul>
 * <p>
 * See Contenedor_de_Celulas.h and Contenedor_de_celulas.cpp for implementation details.
 * </p>
 * 
 * <h3>Cell Cycle System (Ciclo_Modelo)</h3>
 * <p>
 * The Ciclo_Modelo class implements:
 * </p>
 * <ul>
 * <li>Phases of the cell cycle</li>
 * <li>Transitions between phases</li>
 * <li>Rules for cell division and death</li>
 * <li>Phase-specific behaviors</li>
 * </ul>
 * <p>
 * See Ciclo_Modelo.h and Ciclo_Modelo.cpp for implementation details.
 * </p>
 * </div>
 * 
 * @section workflow_sec Simulation Workflow
 * 
 * <div class="workflow">
 * <p>
 * The simulation follows this general workflow:
 * </p>
 * 
 * <h3>1. Initialization</h3>
 * <ul>
 * <li>Load parameters from configuration file</li>
 * <li>Set up microenvironment with substrates</li>
 * <li>Create initial cell population</li>
 * <li>Initialize data structures for simulation</li>
 * </ul>
 * 
 * <h3>2. Main Simulation Loop</h3>
 * <ul>
 * <li>For each time step until end time:
 *   <ul>
 *     <li>Update diffusion of substrates</li>
 *     <li>Calculate and apply cell mechanics (forces, movement)</li>
 *     <li>Process cell biological behaviors (cycle progression, division, death)</li>
 *     <li>Update cell states and interactions</li>
 *     <li>Record data at specified intervals</li>
 *   </ul>
 * </li>
 * </ul>
 * 
 * <h3>3. Output</h3>
 * <ul>
 * <li>Generate cell state data</li>
 * <li>Create visualizations (VTK files)</li>
 * <li>Summarize simulation statistics</li>
 * </ul>
 * </div>
 * 
 * @section build_sec Build Instructions
 * 
 * <div class="build">
 * <p>
 * To build the project:
 * </p>
 * 
 * <ol>
 * <li>Clone the repository</li>
 * <li>From the project root directory, run:
 *    <pre class="fragment">make -f MakefileRo</pre>
 * </li>
 * <li>The executable will be generated in the <code>build/</code> directory as <code>ciclo</code></li>
 * </ol>
 * 
 * <h3>Build System</h3>
 * <ul>
 * <li>All compiled objects (.o files) are placed in the <code>build/</code> directory</li>
 * <li>The executable is generated in the <code>build/</code> directory</li>
 * <li>The build directory is automatically created if it doesn't exist</li>
 * <li>To clean the build files, run <pre class="fragment">make -f MakefileRo clean</pre></li>
 * </ul>
 * </div>
 * 
 * @section run_sec Running a Simulation
 * 
 * <div class="run-simulation">
 * <p>
 * To run a simulation:
 * </p>
 * 
 * <ol>
 * <li>Build the project as described above</li>
 * <li>Run the executable with a parameter file:
 *    <pre class="fragment">./build/ciclo params.ini</pre>
 * </li>
 * <li>The simulation will generate output data based on the parameters specified</li>
 * </ol>
 * </div>
 * 
 * @section output_sec Output Files
 * 
 * <div class="output">
 * <p>
 * The simulation creates output files in a "results" directory (which is created automatically if it doesn't exist):
 * </p>
 * 
 * <ul>
 * <li><code>DatosFinales.dat</code> - Contains summary data for the simulation (time, tumor volume, radius, cell counts)</li>
 * <li><code>Datos_[time].xyz</code> - Contains detailed cell state data at specific time points</li>
 * <li><code>PM1_[time].xml</code> - VTK visualization data for cells (optional)</li>
 * <li><code>HM1_[time].xml</code> - VTK visualization data for microenvironment (optional)</li>
 * </ul>
 * 
 * <p>
 * Additionally, the project includes an <code>/outcomes</code> folder containing example output files and documentation 
 * for visualization and analysis. This folder provides researchers with sample data to understand the expected outputs 
 * and demonstrates how to visualize and interpret simulation results graphically.
 * </p>
 * </div>
 * 
 * @section research_sec Research Applications
 * 
 * <div class="research">
 * <p>
 * This model is particularly suited for cancer research in several key areas:
 * </p>
 * 
 * <ol>
 * <li><b>Tumor Growth Dynamics</b>: Understanding how tumors develop and expand in tissue environments</li>
 * <li><b>Microenvironmental Influences</b>: Studying how factors like hypoxia affect cancer cell behavior</li>
 * <li><b>Heterogeneity</b>: Exploring how diverse cell populations within a tumor contribute to progression and treatment resistance</li>
 * <li><b>Immune Interactions</b>: Investigating how immune cells interact with and potentially eliminate cancer cells</li>
 * <li><b>Therapeutic Response</b>: Modeling how cancer cells respond to treatments like CAR T-cell therapy</li>
 * </ol>
 * 
 * <p>
 * The model has been specifically applied to study the response of heterogeneous tumor-derived organoids to CAR T-cell therapy,
 * as published in Scientific Reports (Nature Portfolio).
 * </p>
 * </div>
 * 
 * @section docs_sec Documentation
 * 
 * <div class="documentation">
 * <p>
 * The source code is thoroughly documented with:
 * </p>
 * 
 * <ul>
 * <li><b>File Headers</b>: Each file includes a detailed description of its purpose and functionality</li>
 * <li><b>Class Documentation</b>: Classes have comprehensive documentation explaining their role in the simulation</li>
 * <li><b>Method Documentation</b>: Methods are documented with descriptions, parameters, and return values</li>
 * <li><b>Research Context</b>: Documentation includes relevant cancer research context</li>
 * </ul>
 * </div>
 * 
 * @section contributors_sec Contributors
 * 
 * <div class="contributors">
 * <p>
 * This software was developed as part of research published in Scientific Reports (2024):
 * </p>
 * 
 * <ul>
 * <li>Luciana Melina Luque - Primary Developer</li>
 * <li>Additional contributors as listed in the publication: <a href="https://www.nature.com/articles/s41598-024-63125-5" target="_blank">https://www.nature.com/articles/s41598-024-63125-5</a></li>
 * </ul>
 * </div>
 * 
 * @section license_sec License
 * 
 * <div class="license">
 * <p>
 * Private - For research purposes only.
 * </p>
 * <p>
 * Copyright © 2024. All rights reserved. This software is provided for research purposes only. 
 * Any commercial use is subject to licensing. For licensing inquiries, please contact the authors.
 * </p>
 * </div>
 * 
 * @section cite_sec How to Cite
 * 
 * <div class="citation">
 * <p>
 * If you use this software in your research, please cite:
 * </p>
 * 
 * <blockquote>
 * <p>
 * <i>In silico study of heterogeneous tumour-derived organoid response to CAR T-cell therapy</i><br>
 * Scientific Reports, volume 14, Article number: 11999 (2024)<br>
 * DOI: <a href="https://www.nature.com/articles/s41598-024-63125-5" target="_blank">https://www.nature.com/articles/s41598-024-63125-5</a>
 * </p>
 * </blockquote>
 * 
 * <p>
 * For work building on the original model design, you may also cite:
 * </p>
 * 
 * <blockquote>
 * <p>
 * <i>Physics-based tissue simulator to model multicellular systems: A study of liver regeneration and hepatocellular carcinoma recurrence</i><br>
 * PLOS Computational Biology, 19(3): e1010920 (2023)<br>
 * DOI: <a href="https://doi.org/10.1371/journal.pcbi.1010920" target="_blank">https://doi.org/10.1371/journal.pcbi.1010920</a>
 * </p>
 * </blockquote>
 * </div>
 */ 