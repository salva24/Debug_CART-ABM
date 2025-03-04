/*! @file Vector.cpp
 *
 * @brief Archivo de implementación de la clase Vector.
 *
 * @author Luciana Melina Luque
 *
 *
 * @details
 * Implementation of vector mathematics operations critical for 3D spatial calculations
 * in cancer simulations. These operations provide the foundation for representing cell
 * positions, mechanical interactions, movement, and microenvironmental gradients.
 * 
 * In cancer research, these vector operations enable:
 * - Precise tracking of tumor cell positions and movements in tissue
 * - Calculation of mechanical forces between cells that drive tumor morphology
 * - Representation of biochemical gradients that influence cell migration
 * - Quantification of spatial characteristics within the tumor microenvironment
 *
 * Inbound Dependencies: Vector.h
 *
 * Outbound Dependencies: Used by mechanical interaction, cell movement, and gradient calculations
 *
 * @license: Open Access - Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
 */
#include "Vector.h"
using namespace std;

/**
 * \brief Default constructor, initializes a zero vector (0,0,0).
 * 
 * Creates a vector with all components set to zero, typically used as an
 * initial state before accumulating forces or positions in cancer simulations.
 */
Vector::Vector(): x(0.0), y(0.0), z(0.0) {}

/**
 * \brief Constructor with specific x, y, z components.
 * \param x1 The x component
 * \param y1 The y component
 * \param z1 The z component
 * 
 * Creates a vector with specific values, commonly used to initialize cell positions,
 * direction vectors, or forces in cancer modeling.
 */
Vector::Vector(double x1, double y1, double z1) : x(x1), y(y1), z(z1) {}

/**
 * \brief Vector addition operator.
 * \param v Vector to add
 * \return Sum vector
 * 
 * Adds two vectors component-wise, essential for combining multiple forces
 * acting on cells or calculating cumulative displacements in tumor growth models.
 */
Vector Vector::operator+(Vector v)
{
    double x1, y1, z1;
    x1 = x + v.x;
    y1 = y + v.y;
    z1 = z + v.z;
    return Vector(x1, y1, z1);
}

/**
 * \brief Vector subtraction operator.
 * \param v Vector to subtract
 * \return Difference vector
 * 
 * Subtracts a vector component-wise, critical for calculating displacement
 * vectors between cells in contact or determining relative positions in
 * tumor spatial analysis.
 */
Vector Vector::operator-(Vector v)
{
    double x1, y1, z1;
    x1 = x - v.x;
    y1 = y - v.y;
    z1 = z - v.z;
    return Vector(x1, y1, z1);
}

/**
 * \brief Vector dot product operator.
 * \param v Vector to dot with
 * \return Scalar dot product result
 * 
 * Computes the dot product, used for calculating projections and angles
 * between cellular movement directions or interaction vectors. Essential for
 * determining alignment of cells and forces in tumor models.
 */
double Vector::operator*(Vector v)
{
    double x1, y1, z1;
    x1 = x * v.x;
    y1 = y * v.y;
    z1 = z * v.z;
    return (x1 + y1 + z1);
}

/**
 * \brief Component-wise vector division.
 * \param v Vector to divide by
 * \return Result of component-wise division
 * 
 * Performs component-wise division, useful for specialized calculations in
 * spatial analysis of tumor structures and normalizing gradients.
 */
Vector Vector::operator/(Vector v)
{
    double x1, y1, z1;
    x1 = x / v.x;
    y1 = y / v.y;
    z1 = z / v.z;
    return Vector(x1, y1, z1);
}

/**
 * \brief Vector scaling (multiplication by scalar).
 * \param d Scale factor
 * \return Scaled vector
 * 
 * Multiplies a vector by a scalar, essential for scaling forces, velocities,
 * or gradients in cancer modeling. Used extensively in mechanical interaction
 * and diffusion calculations.
 */
Vector Vector::operator*(double d)
{
    double x1, y1, z1;
    x1 = x * d;
    y1 = y * d;
    z1 = z * d;
    return Vector(x1, y1, z1);
}

/**
 * \brief Vector division by scalar.
 * \param d Divisor
 * \return Divided vector
 * 
 * Divides a vector by a scalar, useful for normalization or rescaling operations
 * in cell mechanics calculations. Critical for creating unit vectors that
 * represent direction without magnitude.
 */
Vector Vector::operator/(double d)
{
    double x1, y1, z1;
    x1 = x / d;
    y1 = y / d;
    z1 = z / d;
    return Vector(x1, y1, z1);
}

/**
 * \brief Calculates the magnitude (length) of the vector.
 * \return The magnitude as a scalar value
 * 
 * Computes the length of the vector, critical for distance calculations between cells,
 * measuring tumor dimensions, or determining gradient strengths in the microenvironment.
 */
double Vector::modulo(){
    return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}

/**
 * \brief Component-wise addition with scalar.
 * \param d Value to add to each component
 * \return Vector with scalar added to each component
 * 
 * Adds a scalar to each component, used in specialized calculations for
 * tumor microenvironment modeling, particularly for boundary adjustments
 * and offset calculations.
 */
Vector Vector::operator+(double d)
{
    double x1, y1, z1;
    x1 = x + d;
    y1 = y + d;
    z1 = z + d;
    return Vector(x1, y1, z1);
}

/**
 * \brief Component-wise subtraction with scalar.
 * \param d Value to subtract from each component
 * \return Vector with scalar subtracted from each component
 * 
 * Subtracts a scalar from each component, used in specialized calculations
 * for tumor microenvironment modeling and cell position adjustments.
 */
Vector Vector::operator-(double d)
{
    double x1, y1, z1;
    x1 = x - d;
    y1 = y - d;
    z1 = z - d;
    return Vector(x1, y1, z1);
}

/**
 * \brief Normalizes a vector (creates unit vector).
 * \param v Vector to normalize
 * \return Normalized vector (unit length)
 * 
 * Creates a unit vector in the same direction as the input, essential for
 * representing directional information without magnitude, such as cell migration
 * directions or substrate gradient directions in tumor models.
 */
Vector Vector::normaliza(Vector& v)
{
	Vector output = v;
    double norm = 0.0;

    norm += v.x * v.x;
    norm += v.y * v.y;
    norm += v.z * v.z;
    norm = sqrt( norm );

    output.x = v.x / norm;
    output.y = v.y / norm;
    output.z = v.z / norm;

    //  Si la norma es muy chiquita, no tiene sentido normalizar.
    //  En tal caso, seteo todo el vector a 0.0.
    static bool te_lo_adverti = false;
    if( norm <= 1e-16 ){

        if( te_lo_adverti == false ){

            std::cout << "Advertencia: El vector es muy chiquito por lo que se \
            lo normaliz� a 0" << std::endl << std::endl;

            te_lo_adverti = true;
        }

        output.x = 0.0;
        output.y = 0.0;
        output.z = 0.0;
    }

    return output;
}

/**
 * \brief Normalizes a vector in place.
 * \param v Pointer to vector to normalize
 * 
 * Modifies the vector to have unit length, used for efficient normalization
 * in performance-critical code paths of cancer simulations, particularly when
 * calculating directional forces and cell movements.
 */
void Vector::normalizame(Vector* v)
{
    double norm = 1e-32;

    norm += v->x * v->x;
    norm += v->y * v->y;
    norm += v->z * v->z;
    norm = sqrt( norm );

    v->x = v->x / norm;
    v->y = v->y / norm;
    v->z = v->z / norm;

    //  Si la norma es muy chiquita, no tiene sentido normalizar.
    //  En tal caso, seteo todo el vector a 0.0.
    static bool te_lo_adverti = false;
    if( norm <= 1e-16 ){

        if( te_lo_adverti == false ){

            std::cout << "Advertencia: El vector es muy chiquito por lo que se \
            lo normalizo a 0" << std::endl << std::endl;

            te_lo_adverti = true;
        }

        v->x = 0.0;
        v->y = 0.0;
        v->z = 0.0;
    };

}

/**
 * \brief Stream output operator for Vector.
 * \param out Output stream
 * \param v Vector to output
 * \return Modified output stream
 * 
 * Enables printing vector data in a readable format, useful for debugging and
 * visualization of cell positions, forces, or other spatial properties in
 * cancer simulation outputs.
 */
std::ostream& operator<<(std::ostream& out, const Vector& v)
{
    out << v.x << ", ";
    out << v.y << ", ";
    out << v.z << " ";
    out << endl;
    return out;
}

/**
 * \brief Scalar-vector subtraction operator.
 * \param d Scalar value
 * \param v Vector to subtract
 * \return Result vector of d-v
 * 
 * Computes d-v element-wise, used in specialized micro-environment calculations
 * and boundary condition implementations for tumor models.
 */
Vector operator-(double d, Vector v)
{
    double x1, y1, z1;
    x1 = d - v.x;
    y1 = d - v.y;
    z1 = d - v.z;
    return Vector(x1, y1, z1);
}

/**
 * \brief Vector axpy operation (y = y + a*x).
 * \param v Pointer to y vector to be modified
 * \param a Scalar multiplier
 * \param vv x vector to be scaled and added
 * 
 * Performs the operation y = y + a*x in place, a common operation in numerical
 * methods for diffusion, mechanics, and cell movement calculations in cancer models.
 */
void axpy( Vector* v, double& a , Vector& vv )
{
	v->x += a * vv.x;
	v->y += a * vv.y;
	v->z += a * vv.z;
// for( unsigned int i=0; i < (*y).size() ; i++ )
// {
//  (*y)[i] += a * x[i] ;
// }
 return ;
}

//Vector Vector::operator=(Vector v){
//
//    double x1, y1, z1;
//    x1 = v.x;
//    y1 = v.y;
//    z1 = v.z;
//    return Vector(x1, y1, z1);
//	return;
//
//}

/**
 * \brief Component-wise addition of std::vector<double> containers.
 * \param v1 First vector, modified in place (v1 = v1 + v2)
 * \param v2 Second vector to add
 * 
 * Adds two std::vector containers element-wise, used for operating on
 * substrate concentration vectors in tumor microenvironment calculations.
 */
void operator+=( std::vector<double>& v1, const std::vector<double>& v2 )
{
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v1[i] += v2[i]; }
 return;
}

/**
 * \brief Component-wise subtraction of std::vector<double> containers.
 * \param v1 First vector, modified in place (v1 = v1 - v2)
 * \param v2 Second vector to subtract
 * 
 * Subtracts two std::vector containers element-wise, used for 
 * calculating changes in substrate concentrations over time in
 * tumor growth simulations.
 */
void operator-=( std::vector<double>& v1, const std::vector<double>& v2 )
{
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v1[i] -= v2[i]; }
 return;
}

/**
 * \brief Component-wise division of std::vector<double> containers.
 * \param v1 First vector, modified in place (v1 = v1 / v2)
 * \param v2 Second vector as divisor
 * 
 * Divides two std::vector containers element-wise, used for normalizing
 * substrate concentrations or calculating ratios of different biochemical
 * factors in cancer modeling.
 */
void operator/=( std::vector<double>& v1, const std::vector<double>& v2 )
{
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v1[i] /= v2[i]; }
 return;
}

/**
 * \brief Scalar multiplication of std::vector<double>.
 * \param v1 Vector to be scaled in place (v1 = v1 * a)
 * \param a Scalar multiplier
 * 
 * Multiplies each element by a scalar, used for scaling substrate
 * concentrations or adjusting diffusion coefficients in cancer simulations.
 */
void operator*=( std::vector<double>& v1, const double& a )
{
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v1[i] *= a; }
 return;
}

/**
 * \brief Component-wise multiplication of std::vector<double> containers.
 * \param v1 First vector, modified in place (v1 = v1 * v2)
 * \param v2 Second vector to multiply with
 * 
 * Multiplies two std::vector containers element-wise, used for combining
 * factors that have multiplicative effects in cancer biology, such as
 * reaction rates or synergistic factors.
 */
void operator*=( std::vector<double>& v1, const std::vector<double>& v2 )
{
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v1[i] *= v2[i]; }
 return;
}

/**
 * \brief Scalar division of std::vector<double>.
 * \param v1 Vector to be divided in place (v1 = v1 / a)
 * \param a Scalar divisor
 * 
 * Divides each element by a scalar, used for normalizing substrate
 * concentrations or calculating relative values in cancer simulations.
 */
void operator/=( std::vector<double>& v1, const double& a )
{
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v1[i] /= a; }
 return;
}

/**
 * \brief Vector axpy operation for std::vector<double> (y = y + a*x).
 * \param y Pointer to vector to be modified
 * \param a Scalar multiplier
 * \param x Vector to be scaled and added
 * 
 * Performs the operation y = y + a*x in place for substrate vectors,
 * essential in diffusion solvers and biochemical reaction calculations
 * in tumor microenvironment modeling.
 */
void axpy( std::vector<double>* y, double& a , std::vector<double>& x )
{
 for( unsigned int i=0; i < (*y).size() ; i++ )
 {
  (*y)[i] += a * x[i] ;
 }
 return ;
}

/**
 * \brief Vector component-wise axpy operation (y = y + a.*x).
 * \param y Pointer to vector to be modified
 * \param a Vector of multipliers
 * \param x Vector to be scaled and added
 * 
 * Performs the operation y = y + a.*x in place (element-wise multiplication),
 * used for complex substrate interaction calculations in tumor microenvironment
 * models where each component has a different scaling factor.
 */
void axpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x )
{
 for( unsigned int i=0; i < (*y).size() ; i++ )
 {
  (*y)[i] += a[i] * x[i] ;
 }
 return;
}

/**
 * \brief Vector negative axpy operation (y = y - a*x).
 * \param y Pointer to vector to be modified
 * \param a Scalar multiplier
 * \param x Vector to be scaled and subtracted
 * 
 * Performs the operation y = y - a*x in place, used in diffusion solvers and
 * substrate consumption calculations when modeling cancer cell metabolism.
 */
void naxpy( std::vector<double>* y, double& a , std::vector<double>& x )
{
 for( unsigned int i=0; i < (*y).size() ; i++ )
 {
  (*y)[i] -= a * x[i] ;
 }
 return ;
}

/**
 * \brief Vector component-wise negative axpy operation (y = y - a.*x).
 * \param y Pointer to vector to be modified
 * \param a Vector of multipliers
 * \param x Vector to be scaled and subtracted
 * 
 * Performs the operation y = y - a.*x in place (element-wise multiplication),
 * used for modeling differential consumption of multiple substrates in tumor
 * metabolism with varying consumption rates.
 */
void naxpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x )
{
 for( unsigned int i=0; i < (*y).size() ; i++ )
 {
  (*y)[i] -= a[i] * x[i] ;
 }
 return;
}

/**
 * \brief Squared norm of a std::vector<double>.
 * \param v Vector to calculate squared norm
 * \return The squared magnitude of the vector
 * 
 * Calculates the sum of squares of elements, used for efficient distance
 * calculations and convergence testing in cancer simulation algorithms
 * without the computational cost of a square root operation.
 */
double norm_squared( const std::vector<double>& v )
{
 double out = 0.0;
 for( unsigned int i=0 ; i < v.size() ; i++ )
 { out += ( v[i] * v[i] ); }
 return out;
}

/**
 * \brief Squared norm of a Vector.
 * \param v Vector to calculate squared norm
 * \return The squared magnitude of the vector
 * 
 * Calculates the squared length of the vector, used for efficient distance
 * calculations between cells when the exact magnitude isn't required,
 * reducing computational cost in tumor simulations.
 */
double norm_squared( const Vector& v )
{
 double out = 0.0;

 out += v.x * v.x;
 out += v.y * v.y;
 out += v.z * v.z;

 return out;
}

/**
 * \brief Norm (magnitude) of a Vector.
 * \param v Vector to calculate norm
 * \return The magnitude of the vector
 * 
 * Calculates the length of the vector, used for measuring exact distances
 * between cells or determining the strength of forces in cancer simulations.
 */
double norma(const Vector& v){

    return sqrt( norm_squared( v ) );

}

/**
 * \brief Norm (magnitude) of a std::vector<double>.
 * \param v Vector to calculate norm
 * \return The magnitude of the vector
 * 
 * Calculates the Euclidean norm, used for measuring distances between
 * substrate states or determining gradient magnitudes in tumor
 * microenvironment modeling.
 */
double norma( const std::vector<double>& v ){

    return sqrt( norm_squared( v ) );

}
