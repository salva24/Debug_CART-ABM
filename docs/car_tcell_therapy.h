/**
 * @page car_tcell_therapy CAR T-cell Therapy Simulation
 * 
 * @section car_tcell_overview Overview
 * 
 * <div class="car-tcell-overview">
 * <p>
 * This simulation framework is specifically designed to model the complex interactions between
 * CAR T-cells (Chimeric Antigen Receptor T-cells) and tumor-derived organoids. The model incorporates
 * key biological processes and mechanisms involved in CAR T-cell therapy, providing a comprehensive platform
 * for studying treatment efficacy, tumor response, and factors affecting therapeutic outcomes.
 * </p>
 * 
 * <p>
 * As published in <i>Scientific Reports</i> (Nature Portfolio), this model has been applied to study
 * the response of heterogeneous tumor-derived organoids to CAR T-cell therapy, yielding valuable insights
 * into the factors that influence treatment efficacy and potential strategies for optimization.
 * </p>
 * 
 * <blockquote>
 * <p>
 * <b>In silico study of heterogeneous tumour-derived organoid response to CAR T-cell therapy</b><br>
 * <i>Scientific Reports</i>, volume 14, Article number: 11999 (2024)<br>
 * DOI: <a href="https://www.nature.com/articles/s41598-024-63125-5" target="_blank">https://www.nature.com/articles/s41598-024-63125-5</a>
 * </p>
 * </blockquote>
 * </div>
 * 
 * @section car_tcell_biology CAR T-cell Biology
 * 
 * <div class="car-tcell-biology">
 * <p>
 * The simulation implements the following key aspects of CAR T-cell biology:
 * </p>
 * 
 * <h3>CAR T-cell Generation and Properties</h3>
 * <ul>
 * <li><b>Engineered Receptors</b>: T-cells modified with chimeric antigen receptors targeting specific tumor antigens</li>
 * <li><b>Activation Dynamics</b>: Signal transduction processes when CAR T-cells encounter target antigens</li>
 * <li><b>Proliferation Capacity</b>: Self-replication in response to antigen recognition</li>
 * <li><b>Cytokine Production</b>: Release of immune-stimulating factors</li>
 * <li><b>Migration Behavior</b>: Movement through tissues towards chemical gradients</li>
 * </ul>
 * 
 * <h3>T-cell Interactions</h3>
 * <ul>
 * <li><b>Antigen Recognition</b>: Binding to specific surface markers on tumor cells</li>
 * <li><b>Immune Synapse Formation</b>: Creation of specialized contact zones between T-cells and targets</li>
 * <li><b>Cytotoxic Response</b>: Release of granzymes and perforins to induce target cell death</li>
 * <li><b>Serial Killing</b>: Ability to engage and eliminate multiple tumor cells sequentially</li>
 * <li><b>Exhaustion Mechanisms</b>: Reduced functionality after multiple engagements or in inhibitory environments</li>
 * </ul>
 * </div>
 * 
 * @section tumor_organoid Tumor-Derived Organoid Modeling
 * 
 * <div class="tumor-organoid">
 * <p>
 * The simulation represents tumor-derived organoids with these key features:
 * </p>
 * 
 * <h3>Organoid Structure</h3>
 * <ul>
 * <li><b>3D Architecture</b>: Spatially organized cellular arrangements mimicking in vivo tumors</li>
 * <li><b>Heterogeneous Composition</b>: Multiple cell types and states within the organoid</li>
 * <li><b>Growth Dynamics</b>: Realistic expansion patterns based on cell proliferation and death</li>
 * <li><b>Extracellular Matrix</b>: Representation of structural components influencing cell behavior</li>
 * </ul>
 * 
 * <h3>Heterogeneity Aspects</h3>
 * <ul>
 * <li><b>Antigen Expression Variation</b>: Different levels of target antigen expression across tumor cells</li>
 * <li><b>Metabolic Diversity</b>: Varied metabolic states affecting cell survival and response</li>
 * <li><b>Proliferation Rates</b>: Differences in cell cycle progression and division frequency</li>
 * <li><b>Resistance Mechanisms</b>: Various pathways that can reduce susceptibility to T-cell killing</li>
 * </ul>
 * </div>
 * 
 * @section therapy_dynamics Therapy Dynamics and Outcomes
 * 
 * <div class="therapy-dynamics">
 * <p>
 * The simulation enables analysis of various therapy-related processes and outcomes:
 * </p>
 * 
 * <h3>Treatment Efficacy Measures</h3>
 * <ul>
 * <li><b>Tumor Volume Reduction</b>: Changes in overall tumor size over time</li>
 * <li><b>Cell Killing Rates</b>: Quantification of tumor cell elimination by CAR T-cells</li>
 * <li><b>Survival Analysis</b>: Estimation of long-term treatment outcomes</li>
 * <li><b>Antigen Escape</b>: Emergence of antigen-negative populations</li>
 * </ul>
 * 
 * <h3>Resistance Mechanisms</h3>
 * <ul>
 * <li><b>Target Antigen Downregulation</b>: Reduced expression of CAR-targeted antigens</li>
 * <li><b>Immunosuppressive Microenvironment</b>: Production of inhibitory signals</li>
 * <li><b>Physical Barriers</b>: Limited T-cell infiltration due to tissue architecture</li>
 * <li><b>Clonal Selection</b>: Outgrowth of resistant subpopulations</li>
 * </ul>
 * </div>
 * 
 * @section model_parameters Key Model Parameters
 * 
 * <div class="model-parameters">
 * <p>
 * The CAR T-cell therapy simulation is controlled by several key parameters:
 * </p>
 * 
 * <h3>T-cell Parameters</h3>
 * <ul>
 * <li><b>Injection Timing</b>: When T-cells are introduced into the simulation</li>
 * <li><b>T-cell Count</b>: Initial number of CAR T-cells</li>
 * <li><b>Binding Affinity</b>: Strength of CAR-antigen interaction</li>
 * <li><b>Killing Efficiency</b>: Probability of inducing apoptosis after engagement</li>
 * <li><b>Exhaustion Rate</b>: How quickly T-cells lose functionality</li>
 * <li><b>Motility Parameters</b>: Speed and directional persistence of T-cell movement</li>
 * </ul>
 * 
 * <h3>Tumor Parameters</h3>
 * <ul>
 * <li><b>Antigen Expression Levels</b>: Distribution of target antigen across tumor cells</li>
 * <li><b>Heterogeneity Profile</b>: Specification of subpopulation characteristics</li>
 * <li><b>Growth Rate</b>: Proliferation speed of tumor cells</li>
 * <li><b>Resistance Probabilities</b>: Likelihood of developing resistance mechanisms</li>
 * </ul>
 * </div>
 * 
 * @section research_insights Research Insights
 * 
 * <div class="research-insights">
 * <p>
 * The simulation has provided several important insights into CAR T-cell therapy for solid tumors:
 * </p>
 * 
 * <ul>
 * <li><b>Heterogeneity Impact</b>: Tumor heterogeneity significantly affects treatment outcomes, with more heterogeneous tumors showing increased resistance</li>
 * <li><b>Infiltration Barriers</b>: Physical access of T-cells to tumor cells is a critical limiting factor in solid tumor treatment</li>
 * <li><b>Antigen Escape</b>: The model predicts patterns of antigen-loss variants emerging during therapy</li>
 * <li><b>Exhaustion Effects</b>: T-cell exhaustion substantially limits efficacy in densely packed tumors</li>
 * <li><b>Optimal Dosing</b>: Insights into ideal T-cell dosing strategies for different tumor compositions</li>
 * </ul>
 * </div>
 * 
 * @section future_directions Future Directions
 * 
 * <div class="future-directions">
 * <p>
 * Ongoing and future development of the CAR T-cell therapy simulation includes:
 * </p>
 * 
 * <ul>
 * <li><b>Combination Therapies</b>: Modeling CAR T-cells in combination with other treatment modalities</li>
 * <li><b>Advanced T-cell Designs</b>: Simulating next-generation CAR designs with enhanced properties</li>
 * <li><b>Personalized Medicine</b>: Adapting the model to patient-specific tumor characteristics</li>
 * <li><b>Immune Microenvironment</b>: Incorporating additional immune components like regulatory T-cells and myeloid cells</li>
 * <li><b>Validation Studies</b>: Correlating simulation predictions with experimental and clinical data</li>
 * </ul>
 * </div>
 */ 