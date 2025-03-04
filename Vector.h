/*!
 * @file Vector.h
 *
 * @brief Vector mathematics definitions for 3D spatial calculations in cancer simulation.
 *
 * @author Luciana Melina Luque
 *
 * @details
 * This file defines the Vector class and related operations that provide the mathematical
 * foundation for all spatial calculations in the cancer simulation model. Vectors are
 * essential for representing cell positions, forces, movement directions, and gradients
 * within the tumor microenvironment. The implementation supports efficient 3D vector
 * operations needed for modeling cellular mechanics, movement, and interactions.
 * 
 * In cancer research applications, these vector operations enable:
 * - Precise tracking of tumor cell positions and movements
 * - Calculation of mechanical forces between cells that drive tumor morphology
 * - Representation of biochemical gradients that influence cell migration and phenotype
 * - Quantification of spatial heterogeneity within the tumor microenvironment
 * 
 * Inbound Dependencies: <cmath>, <iostream>, <vector>
 *
 * Outbound Dependencies: Used by Celula.h, Microambiente.h, Contenedor_de_Celulas.h
 * 
 * @license: Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#ifndef _VECTOR_H
#define _VECTOR_H

#include <iostream>
#include <vector>
#include <cmath>

/**
 * \brief 3D vector class for spatial calculations in cancer simulation.
 * 
 * The Vector class provides the mathematical foundation for all spatial operations 
 * in the cancer simulation model. It implements a standard 3D vector with x, y, z 
 * components and operations like addition, subtraction, multiplication, normalization,
 * and other vector algebra operations.
 * 
 * In cancer modeling, vectors represent key spatial properties:
 * - Cell positions in 3D tissue space
 * - Forces acting between cells (adhesion, repulsion)
 * - Cell movement directions
 * - Substrate gradients in the microenvironment
 * - Morphological measurements of tumor structure
 */
class Vector {
private:

public:
	double x, y, z; ///< The x, y, z components of the 3D vector

    /**
     * \brief Default constructor, initializes vector to (0,0,0).
     * 
     * Creates a zero vector, typically used as an initial state before accumulating
     * forces or positions in cancer simulations.
     */
    Vector();
    
    /**
     * \brief Constructor with specific x, y, z components.
     * \param x1 The x component
     * \param y1 The y component
     * \param z1 The z component
     * 
     * Creates a vector with specific values, used to initialize cell positions,
     * direction vectors, or forces in cancer modeling.
     */
    Vector(double x1, double y1, double z1);
    
    /**
     * \brief Vector addition operator.
     * \param v Vector to add
     * \return Sum vector
     * 
     * Adds two vectors, used for combining forces or displacements in tumor growth models.
     */
    Vector operator+(Vector v);
    
    /**
     * \brief Vector subtraction operator.
     * \param v Vector to subtract
     * \return Difference vector
     * 
     * Subtracts a vector, used for calculating displacement vectors between cells
     * or determining relative positions in tumor spatial analysis.
     */
    Vector operator-(Vector v);
    
    /**
     * \brief Vector dot product operator.
     * \param v Vector to dot with
     * \return Scalar dot product result
     * 
     * Computes the dot product, used for calculating projections and angles 
     * between cellular movement directions or interaction vectors.
     */
     double operator*(Vector v);
     
    /**
     * \brief Component-wise vector division.
     * \param v Vector to divide by
     * \return Result of component-wise division
     * 
     * Performs component-wise division, useful for specialized calculations in 
     * spatial analysis of tumor structures.
     */
     Vector operator/(Vector v);
    
    /**
     * \brief Vector scaling (multiplication by scalar).
     * \param d Scale factor
     * \return Scaled vector
     * 
     * Multiplies a vector by a scalar, essential for scaling forces, 
     * velocities, or gradients in cancer modeling.
     */
    Vector operator*(double d);
    
    /**
     * \brief Vector division by scalar.
     * \param d Divisor
     * \return Divided vector
     * 
     * Divides a vector by a scalar, useful for normalization or 
     * rescaling operations in cell mechanics calculations.
     */
    Vector operator/(double d);
    
    /**
     * \brief Component-wise addition with scalar.
     * \param d Value to add to each component
     * \return Vector with scalar added to each component
     * 
     * Adds a scalar to each component, used in specialized calculations
     * for tumor microenvironment modeling.
     */
    Vector operator+(double d);
    
    /**
     * \brief Component-wise subtraction with scalar.
     * \param d Value to subtract from each component
     * \return Vector with scalar subtracted from each component
     * 
     * Subtracts a scalar from each component, used in specialized
     * calculations for tumor microenvironment modeling.
     */
    Vector operator-(double d);
    
    /**
     * \brief Calculates the magnitude (length) of the vector.
     * \return The magnitude as a scalar value
     * 
     * Computes the length of the vector, critical for distance calculations
     * between cells, measuring tumor dimensions, or determining gradient strengths.
     */
    double modulo();
    
    /**
     * \brief Normalizes a vector (unit length).
     * \param v Vector to normalize
     * \return Normalized vector (unit length)
     * 
     * Creates a unit vector in the same direction, essential for representing
     * directional information without magnitude, such as cell migration directions
     * or substrate gradient directions in tumor models.
     */
    Vector normaliza(Vector& v);
    
    /**
     * \brief Normalizes a vector in place.
     * \param v Pointer to vector to normalize
     * 
     * Modifies the vector to have unit length, used for efficient normalization
     * in performance-critical code paths of cancer simulations.
     */
    void normalizame(Vector* v);
    
    /**
     * \brief Stream output operator for Vector.
     * \param out Output stream
     * \param v Vector to output
     * \return Modified output stream
     * 
     * Enables printing vector data for debugging and visualization of
     * cell positions, forces, or other spatial properties in cancer models.
     */
    friend std::ostream& operator<<(std::ostream& out, const Vector& v);
    
    /**
     * \brief Scalar-vector subtraction operator.
     * \param d Scalar value
     * \param v Vector to subtract
     * \return Result vector of d-v
     * 
     * Computes d-v element-wise, used in specialized micro-environment
     * calculations and boundary condition implementations.
     */
    friend Vector operator-(double d, Vector v); // Resta d - v
    
    /**
     * \brief Vector axpy operation (y = y + a*x).
     * \param v Pointer to y vector to be modified
     * \param a Scalar multiplier
     * \param vv x vector to be scaled and added
     * 
     * Performs the mathematical operation y = y + a*x in place,
     * a common operation in numerical methods for diffusion and mechanics
     * calculations in cancer modeling.
     */
	friend void axpy( Vector* v, double& a , Vector& vv );
};

/**
 * \brief Component-wise addition of std::vector<double> containers.
 * \param v1 First vector, modified in place (v1 = v1 + v2)
 * \param v2 Second vector to add
 * 
 * Adds two std::vector containers element-wise, used for operating on
 * substrate concentration vectors in the cancer microenvironment.
 */
void operator+=( std::vector<double>& v1, const std::vector<double>& v2 );

/**
 * \brief Component-wise subtraction of std::vector<double> containers.
 * \param v1 First vector, modified in place (v1 = v1 - v2)
 * \param v2 Second vector to subtract
 * 
 * Subtracts two std::vector containers element-wise, used for 
 * calculating changes in substrate concentrations over time.
 */
void operator-=( std::vector<double>& v1, const std::vector<double>& v2 );

/**
 * \brief Component-wise division of std::vector<double> containers.
 * \param v1 First vector, modified in place (v1 = v1 / v2)
 * \param v2 Second vector as divisor
 * 
 * Divides two std::vector containers element-wise, used for normalizing
 * substrate concentrations or calculating ratios of different factors.
 */
void operator/=( std::vector<double>& v1, const std::vector<double>& v2 );

/**
 * \brief Scalar multiplication of std::vector<double>.
 * \param v1 Vector to be scaled in place (v1 = v1 * a)
 * \param a Scalar multiplier
 * 
 * Multiplies each element by a scalar, used for scaling substrate
 * concentrations or adjusting diffusion coefficients in cancer models.
 */
void operator*=( std::vector<double>& v1, const double& a );

/**
 * \brief Component-wise multiplication of std::vector<double> containers.
 * \param v1 First vector, modified in place (v1 = v1 * v2)
 * \param v2 Second vector to multiply with
 * 
 * Multiplies two std::vector containers element-wise, used for
 * combining factors that have multiplicative effects in cancer biology.
 */
void operator*=( std::vector<double>& v1, const std::vector<double>& v2 );

/**
 * \brief Scalar division of std::vector<double>.
 * \param v1 Vector to be divided in place (v1 = v1 / a)
 * \param a Scalar divisor
 * 
 * Divides each element by a scalar, used for normalizing substrate
 * concentrations or calculating relative values in cancer simulations.
 */
void operator/=( std::vector<double>& v1, const double& a );

/**
 * \brief Vector axpy operation for std::vector<double> (y = y + a*x).
 * \param y Pointer to vector to be modified
 * \param a Scalar multiplier
 * \param x Vector to be scaled and added
 * 
 * Performs the mathematical operation y = y + a*x in place for substrate
 * vectors, used in diffusion solvers and biochemical reaction calculations.
 */
void axpy( std::vector<double>* y, double& a , std::vector<double>& x );

/**
 * \brief Vector component-wise axpy operation (y = y + a.*x).
 * \param y Pointer to vector to be modified
 * \param a Vector of multipliers
 * \param x Vector to be scaled and added
 * 
 * Performs the mathematical operation y = y + a.*x in place (element-wise),
 * used for complex substrate interaction calculations in cancer microenvironment.
 */
void axpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x );

/**
 * \brief Vector negative axpy operation (y = y - a*x).
 * \param y Pointer to vector to be modified
 * \param a Scalar multiplier
 * \param x Vector to be scaled and subtracted
 * 
 * Performs the mathematical operation y = y - a*x in place, used in
 * diffusion solvers and substrate consumption calculations.
 */
void naxpy( std::vector<double>* y, double& a , std::vector<double>& x );

/**
 * \brief Vector component-wise negative axpy operation (y = y - a.*x).
 * \param y Pointer to vector to be modified
 * \param a Vector of multipliers
 * \param x Vector to be scaled and subtracted
 * 
 * Performs the mathematical operation y = y - a.*x in place (element-wise),
 * used for modeling differential consumption of multiple substrates in cancer.
 */
void naxpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x );

/**
 * \brief Squared norm of a std::vector<double>.
 * \param v Vector to calculate squared norm
 * \return The squared magnitude of the vector
 * 
 * Calculates the sum of squares of elements, used for efficient distance
 * calculations and convergence testing in cancer simulation algorithms.
 */
double norm_squared( const std::vector<double>& v );

/**
 * \brief Squared norm of a Vector.
 * \param v Vector to calculate squared norm
 * \return The squared magnitude of the vector
 * 
 * Calculates the squared length of the vector, used for efficient distance
 * calculations when the exact magnitude isn't required.
 */
double norm_squared( const Vector& v );

/**
 * \brief Norm (magnitude) of a std::vector<double>.
 * \param v Vector to calculate norm
 * \return The magnitude of the vector
 * 
 * Calculates the Euclidean norm, used for measuring distances
 * between substrate states or determining gradient magnitudes.
 */
double norma( const std::vector<double>& v );

/**
 * \brief Norm (magnitude) of a Vector.
 * \param v Vector to calculate norm
 * \return The magnitude of the vector
 * 
 * Calculates the length of the vector, used for measuring distances between
 * cells or determining the strength of forces in cancer simulations.
 */
double norma( const Vector& v );

#endif
