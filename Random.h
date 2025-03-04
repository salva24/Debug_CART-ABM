/**
 * @file Random.h
 *
 * @author Luciana Melina Luque
 *
 * @brief Defines the RNG class for stochastic processes in cancer simulations
 *
 * @details
 * This file contains the definition of the RNG (Random Number Generator) class, which provides 
 * various probability distributions for implementing stochastic processes in cancer simulations.
 * Randomness is essential for realistic cancer modeling, as many cellular processes exhibit 
 * inherent stochasticity, including:
 * 
 * - Cell cycle phase transitions (G1→S, G2→M)
 * - Cell death decisions (apoptosis, necrosis)
 * - Cell migration and motility patterns
 * - Mutation and phenotypic variation
 * - Treatment response heterogeneity
 * 
 * The implementation supports both uniform and normal distributions with configurable parameters,
 * enabling researchers to model different levels of stochasticity in various cancer processes.
 * 
 * Inbound Dependencies: iostream, random
 * Outbound Dependencies: Multiple model components that require stochastic behavior
 *
 * @licence: Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */

#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <iostream>
#include <random>

/**
 * @brief Random Number Generator for stochastic processes in cancer modeling
 * @details
 * The RNG class provides controlled stochastic behavior for cancer simulations,
 * enabling the realistic modeling of probabilistic cellular processes. In cancer research,
 * stochasticity is critical for accurately representing biological variability,
 * from cell-to-cell differences in proliferation rates to the randomness in treatment responses.
 * 
 * Key cancer research applications include:
 * - Modeling cell cycle progression with realistic variability
 * - Simulating heterogeneous responses to microenvironmental conditions
 * - Implementing probabilistic cell fate decisions (division, death, differentiation)
 * - Generating diverse motility patterns for invasion and metastasis studies
 * - Creating realistic tumor heterogeneity through stochastic phenotype switching
 * 
 * The class uses modern C++ random number generation facilities (mt19937 engine)
 * to provide high-quality pseudorandom sequences with long periods and good
 * statistical properties, ensuring reliable simulation results for cancer research.
 */
class RNG {
    private:
        /**
         * @brief Maximum value for uniform distribution
         * @details Sets the upper bound for the uniform random distribution. In cancer
         * modeling, this parameter helps define the range of possible values for
         * stochastic processes like cell migration distance or division timing.
         */
        double n_max;
        
        /**
         * @brief Seed for the random number generator
         * @details Controls the initialization of the random sequence. Using consistent
         * seeds enables reproducible simulation runs, which is essential for systematic
         * studies of cancer progression and treatment response.
         */
        int seed;
        
        /**
         * @brief Mean value for normal distribution
         * @details Specifies the center of the normal distribution. In cancer modeling,
         * this typically represents the average value of biological parameters such as
         * cell cycle duration, apoptosis threshold, or migration speed.
         */
        double media;
        
        /**
         * @brief Standard deviation for normal distribution
         * @details Controls the spread of the normal distribution. In cancer research,
         * this parameter models biological variability between cells, allowing for
         * realistic heterogeneity in tumor populations.
         */
        double std_desv;
        
        /**
         * @brief Mersenne Twister random number engine
         * @details High-quality pseudorandom number generator with long period and
         * good statistical properties. Ensures reliable stochastic simulations for
         * cancer research where large numbers of random values are needed.
         */
        std::mt19937 generator;
        
        /**
         * @brief Uniform distribution over [0, n_max]
         * @details Distribution for generating uniformly distributed random numbers.
         * Used for processes where any value in a range is equally likely, such as
         * initial cell positioning or random parameter selection.
         */
        std::uniform_real_distribution<double> dist;
        
        /**
         * @brief Secondary uniform distribution (parameters set at call time)
         * @details Used by the RandomNumber(min,max) method for generating values
         * in specified ranges. Enables flexible stochastic modeling of different
         * cancer processes with varying parameters.
         */
        std::uniform_real_distribution<double> dist2;
        
        /**
         * @brief Normal distribution with preset mean and standard deviation
         * @details Distribution for generating normally distributed random numbers.
         * In cancer modeling, normal distributions often better represent biological
         * variation in parameters like cell cycle duration or response thresholds.
         */
        std::normal_distribution<double> d;
        
        /**
         * @brief Secondary normal distribution (parameters set at call time)
         * @details Used by the NormalRandom_CM method for generating values with
         * specified mean and standard deviation. Allows dynamic configuration of
         * stochastic processes throughout cancer simulations.
         */
        std::normal_distribution<double> d2;

    public:
         /**
          * @brief Default constructor
          * @details
          * Initializes the RNG with default parameters:
          * - Maximum value: 1.0
          * - Seed: 13
          * - Mean: 1.0
          * - Standard deviation: 0.25
          * 
          * These defaults provide reasonable stochastic behavior for many cancer
          * simulation scenarios without requiring explicit configuration.
          */
         RNG();

         /**
          * @brief Parameterized constructor
          * @param NM Maximum value for uniform distribution
          * @param SD Seed for the random generator
          * @param med Mean for normal distribution
          * @param desv Standard deviation for normal distribution
          * @details
          * Creates an RNG with custom parameters, allowing researchers to tailor
          * the stochastic behavior to specific cancer modeling requirements.
          * Different cancer types or experimental scenarios may require distinct
          * probability distributions to accurately represent biological variability.
          */
         RNG(double NM, int SD, double med, double desv);

        /**
         * @brief Generate a uniform random number between 0 and n_max
         * @return Random value in [0, n_max]
         * @details
         * Provides uniform randomness for cancer simulation processes where any
         * value in the range is equally likely. Common applications include random
         * cell selection, stochastic decision making, and uniform spatial sampling.
         */
        double RandomNumber();
        
        /**
         * @brief Generate a normally distributed random number
         * @return Random value from normal distribution with preset mean and standard deviation
         * @details
         * Provides normally distributed randomness, which better represents many
         * biological parameters in cancer. Used for modeling inherent variation in
         * cell properties, cycle times, and response thresholds, creating realistic
         * tumor heterogeneity.
         */
        double NormalRandom();
        
        /**
         * @brief Generate a normally distributed random number with custom parameters
         * @param mean Center of the normal distribution
         * @param standard_deviation Spread of the normal distribution
         * @return Random value from specified normal distribution
         * @details
         * Allows dynamic specification of normal distribution parameters for different
         * stochastic processes. In cancer modeling, this enables context-specific
         * randomness, such as different variability in hypoxic versus normoxic regions
         * or different cell populations within a heterogeneous tumor.
         */
        double NormalRandom_CM( double mean, double standard_deviation );
        
        /**
         * @brief Generate a uniform random number in a custom range
         * @param min Lower bound of the range
         * @param max Upper bound of the range
         * @return Random value in [min, max]
         * @details
         * Provides uniform randomness within a specified range. In cancer research,
         * this is used for processes with defined boundaries, such as constrained
         * migration distances, parameter sampling within biologically relevant ranges,
         * or probabilistic selection between discrete cellular states.
         */
        double RandomNumber( double min, double max );

};


#endif
