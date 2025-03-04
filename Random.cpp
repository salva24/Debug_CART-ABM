/**
 * @file Random.cpp
 *
 * @author Luciana Melina Luque
 *
 * @brief Implementation of the RNG class for stochastic processes in cancer simulations
 *
 * @details
 * This file implements the functionality of the RNG (Random Number Generator) class,
 * which provides the stochastic foundation for modeling probabilistic cellular behaviors
 * in cancer simulations. The implementation supports both uniform and normal distributions,
 * essential for representing various sources of randomness in cancer biology:
 * 
 * - Cell cycle progression randomness
 * - Heterogeneous responses to therapies and microenvironmental conditions
 * - Stochastic cell migration and adhesion behaviors
 * - Probabilistic mutation and phenotypic switching
 * - Variable treatment responses across cells within a tumor
 * 
 * The implementation uses modern C++ random facilities (mt19937, distributions) to
 * ensure high-quality random sequences with appropriate statistical properties for
 * scientifically valid cancer research simulations.
 * 
 * Inbound Dependencies: Random.h
 * Outbound Dependencies: None
 *
 * @license: Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#include "Random.h"

/**
 * @brief Default constructor implementation
 * @details
 * Initializes the RNG with cancer research-relevant default values:
 * - Maximum value: 1.0 (appropriate for probability-based decisions)
 * - Seed: 13 (fixed for reproducibility)
 * - Mean: 1.0 (baseline for cell parameters)
 * - Standard deviation: 0.25 (moderate biological variability)
 * 
 * These values were selected to provide appropriate stochasticity for standard
 * cancer simulation scenarios, creating realistic heterogeneity without extreme values.
 */
RNG::RNG() : n_max(1.0), seed(13), media(1.0), std_desv(0.25) {}

/*
RNG::RNG(double NM, int SD) : n_max(NM), seed(SD),
    dist(std::uniform_real_distribution<double>(0.0, NM)) {
        generator.seed(seed);
}
*/

/**
 * @brief Parameterized constructor implementation
 * @param NM Maximum value for uniform distribution
 * @param SD Seed for the random generator
 * @param med Mean for normal distribution
 * @param desv Standard deviation for normal distribution
 * @details
 * Creates a customized RNG for specific cancer modeling requirements. The constructor:
 * 1. Initializes member variables with the provided parameters
 * 2. Creates appropriate distribution objects
 * 3. Seeds the generator for reproducible sequences
 * 
 * In cancer research, custom parameterization enables modeling specific tumor types
 * or experimental conditions with their unique stochastic characteristics. For example:
 * - Different cancer cell lines show unique variability in proliferation rates
 * - Metastatic cells typically exhibit greater migration randomness
 * - Treatment-resistant populations may show broader response distributions
 */
RNG::RNG(double NM, int SD, double med, double desv) : n_max(NM), seed(SD), media(med), std_desv(desv), 
dist(std::uniform_real_distribution<double>(0.0, NM)), d(std::normal_distribution<double>(media,std_desv)) {
        generator.seed(seed);
}


/**
 * @brief Implementation of uniform random number generation
 * @return Uniform random value between 0 and n_max
 * @details
 * Generates a uniformly distributed random number using the configured distribution and generator.
 * In cancer modeling, uniform distributions are particularly useful for:
 * - Stochastic cell fate decisions (using thresholds)
 * - Random spatial positioning of cells in initial tumor configurations
 * - Monte Carlo sampling of parameter spaces in model sensitivity analysis
 * - Equal-probability selection of cellular behaviors in phenotype switching models
 */
double RNG::RandomNumber() {
    return dist(generator);
};

/**
 * @brief Implementation of normal random number generation
 * @return Normally distributed random value with preset mean and standard deviation
 * @details
 * Generates a normally distributed random number using the configured normal distribution.
 * Normal distributions are essential in cancer research for modeling:
 * - Biological variation in cell cycle duration (proliferation heterogeneity)
 * - Probabilistic response to therapy (with mean representing average efficacy)
 * - Realistic variation in motility parameters across tumor cell populations
 * - Parameter variation when initializing heterogeneous tumor cell populations
 */
double RNG::NormalRandom(){
    return d(generator);
}


/**
 * @brief Implementation of normal random number generation with custom parameters
 * @param mean Center of the normal distribution
 * @param standard_deviation Spread of the normal distribution
 * @return Normally distributed random value with specified parameters
 * @details
 * Creates a temporary normal distribution with the specified parameters and generates a value.
 * In cancer modeling, this flexibility is crucial for context-dependent randomness:
 * - Different variation in hypoxic versus normoxic regions of a tumor
 * - Distinct stochasticity profiles for different cell types (cancer cells vs. stroma)
 * - Modeling variation in treatment response across heterogeneous tumor subpopulations
 * - Representing changes in parameter stochasticity during disease progression
 */
double RNG::NormalRandom_CM( double mean, double standard_deviation )
{
	std::normal_distribution<double> d2(mean,standard_deviation);
	return d2(generator);
}

/**
 * @brief Implementation of uniform random number generation in custom range
 * @param min Lower bound of the range
 * @param max Upper bound of the range
 * @return Uniform random value between min and max
 * @details
 * Creates a temporary uniform distribution with the specified range and generates a value.
 * This flexibility is important in cancer modeling for:
 * - Sampling parameters within biologically relevant bounds
 * - Implementing bounded stochastic processes (e.g., confined migration)
 * - Modeling decision processes with non-standard probability ranges
 * - Dynamic adjustment of randomness parameters based on tumor state
 */
double RNG::RandomNumber( double min, double max )
{
	std::uniform_real_distribution<double> dist2(min,max);
	return dist2(generator);
}
