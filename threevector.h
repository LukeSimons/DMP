/** @file threevector.h
 *  @brief Class defining a threevector object for use in physics models
 *
 *  A class which deals with linear algebra of vectors with three dimensions.
 *
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug No known bugs.
 */

#ifndef __THREEVECTOR_H_INCLUDED__
#define __THREEVECTOR_H_INCLUDED__

#include <math.h>
#include <fstream>
#include <iostream>

#include "Constants.h"

using namespace dimplconsts;

/** @class threevector
 *  @brief Class defining a threevector object for use in physics models
 *
 *  A class which deals with linear algebra of vectors with three dimensions.
 */
class threevector{
    protected:
        /** @name Member Data
         *  @brief Private data defining elements of vector
         *
         *  private, three coordinates defining the three vector
         */
        ///@{
        double xcoord;
        double ycoord;
        double zcoord;
        ///@}

    public:
        /** @name Constructors
         *  @brief functions to construct threevector class
         */
        ///@{
        /** @brief Default constructor.
         *
         *  Set all member data to zero
         */
        threevector();
        /** @brief Parameterised Cartesian constructor.
         *  @param x defines the first coordinate
         *  @param y defines the second coordinate
         *  @param z defines the third coordinate
         */
        threevector(double x, double y, double z);
        /** @brief Parameterised Polar constructor.
         *  @param r defines the first coordinate
         *  @param theta defines the second coordinate
         *  @param phiorz defines the third coordinate
         *  @param type Possible values: 's'pherical or 'c'ylindrical
         */
        threevector(double r, double theta, double phiorz, char type);
        ///@}

        /** @name Print Functions
         *  @brief functions to print data in threevector
         */
        ///@{
        /** @brief Print coordinates to in vectorised format to screen
         */
        void print();

        /** @brief Print coordinates to in vectorised format to file
         *  @param fout ofstream file stream to print to
         */
        void print(std::ofstream &fout);
        ///@}

        /** @name Mutator Functions
         *  @brief change value of coordinates
         *  @param xvalue Defines the value to overwrite the first coordinate
         *  @param yvalue Defines the value to overwrite the second coordinate
         *  @param zvalue Defines the value to overwrite the third coordinate
         */
        ///@{
        inline void setx(double xvalue){ xcoord = xvalue; };
        inline void sety(double yvalue){ ycoord = yvalue; };
        inline void setz(double zvalue){ zcoord = zvalue; };
        ///@}

        /** @name Accessor functions
         * get the value of a particular coordinate
         */
        ///@{
        //* @var xcoord get the first coordinate */
        inline double getx()const{ return xcoord; };
        //* @var ycoord get the second coordinate */
        inline double gety()const{ return ycoord; };
        //* @var zcoord get the third coordinate */
        inline double getz()const{ return zcoord; };
        /**
         *  get value of third coordinate converted to polar coordinates
         *  @see mag3()
         *  @return zcoord/mag3() the third coordinate in polar coordinates
         */
        inline double gettheta()const{ return acos(zcoord/mag3()); }

        /**
         *  get value of second polar coordinate from cartesian coordinates
         *  @return Value of second polar coordinate, phi
         */
        inline double getphi()const{
            double phi = atan2(ycoord,xcoord);
            if (phi<0.0) {phi += 2.0*PI;} // atan2 range is -pi to +pi
            return phi;
        };

        /**
         *  get value of cylindrical coordinate rho from cartesian coordinates
         *  @return Value of cylindrical coordinate rho
         */
        inline double getrho()const{
            return sqrt(xcoord*xcoord + ycoord*ycoord);
        }
        ///@}

        /** @brief Square of vector
         *  Return square of the vector
         *  @return The squared sum of three coordinates
         */
        inline double square()const{
            return xcoord*xcoord + ycoord*ycoord + zcoord*zcoord;
        }

        /** @brief Return magnitude of vector
         *  Return magnitude of the vector implementing square() function
         *  @see square()
         *  @return The squared sum of three coordinates
         */
        inline double mag3()const{
            return sqrt(square());
        }

        /** @brief Return unitvector
         *
         *  Return a vector of magnitude 1 in same direction as vector
         *  @see mag3()
         *  @return the unit vector of this vector
         */
        inline threevector getunit()const{

            double Magnitude = mag3();
            if( Magnitude != 0.0 ){
                threevector unitvec(xcoord/Magnitude,ycoord/Magnitude,zcoord/Magnitude);
                return unitvec;
            }else{
                threevector unitvec(0.0,0.0,0.0);
                return unitvec;
            }
        };

        /** @name Overload Operators
         *  @brief Overload the arithmetic operators to perform linear algebra
         *  operations on the threevector class.
         *  @see getx()
         *  @see gety()
         *  @see getz()
         */
        ///@{
        /**@brief Overload addition operator
         *
         * Return threevector sum of this vector and v_old
         * @param v_old threevector to be added to this vector
         * @return the sum of this vector and v_old using linear algebra
         */
        inline threevector operator+(const threevector v_old)const{
            threevector v_sum(v_old.getx() + xcoord, v_old.gety() + ycoord,
                v_old.getz() + zcoord);
            return v_sum;
        };

        /**@brief Overload subtraction operator
         *
         * Return threevector difference of this vector and v_old
         * @param v_old threevector to be subtracted from this vector
         * @return the difference of this vector and v_old using linear algebra
         */
        inline threevector operator-(const threevector v_old)const{
            threevector v_sum(xcoord - v_old.getx(), ycoord - v_old.gety(),
                zcoord - v_old.getz());
            return v_sum;
        };

        /**@brief Overload addition assignment operator
         *
         * Mutate this threevector as addition of this vector and v_old
         * @param v_old threevector to be added to this vector
         */
        inline void operator+=(threevector v_old){
            xcoord += v_old.getx();
            ycoord += v_old.gety();
            zcoord += v_old.getz();
        };

        /**@brief Overload subtraction assignment operator
         *
         * Mutate this threevector as subtraction of this vector and v_old
         * @param v_old threevector to be subtracted from this vector
         */
        inline void operator-=(threevector v_old){
            xcoord -= v_old.getx();
            ycoord -= v_old.gety();
            zcoord -= v_old.getz();
        };

        /**@brief Overload assignment operator
         *
         * Assign this threevector as v_old
         * @param v_old threevector to be assigned to this vector
         */
        inline void operator=(threevector v_old){
            xcoord = v_old.getx();
            ycoord = v_old.gety();
            zcoord = v_old.getz();
        };

        /**@brief Overload multiplication operator, scalar
         *
         * return threevector as the product of this vector and a scalar
         * @param scalar the multiple for this vector
         * @return threevector scalar multiple of this vector and scalar
         */
        inline threevector operator*(double scalar)const{
            threevector v_new(xcoord*scalar,ycoord*scalar,zcoord*scalar);
            return v_new;
        };

        /**@brief Overload multiplication operator, threevector
         *
         * return the inner product of this vector and v_old threevector
         * @param v_old threevector to be used for inner product
         * @return the inner product of this vector and v_old
         */
        inline double operator*(threevector v_old)const{
            return xcoord*v_old.getx() + ycoord*v_old.gety()
            + zcoord*v_old.getz();
        };

        /**@brief Overload binary operator
         *
         * return the outer product of this vector and v_old threevector
         * @param v_old threevector to be used for outer product
         * @return threevector defines outer product of this vector and v_old
         */
        inline threevector operator^(threevector v_old)const{
            threevector v_new((ycoord * v_old.getz()) - (zcoord * v_old.gety()),
                (zcoord * v_old.getx()) - (xcoord * v_old.getz()),
                (xcoord * v_old.gety()) - (ycoord * v_old.getx()));
            return v_new;
        };

        /**@brief Overload bitwise leftshit operator
         *
         * return write the data of the three vector to ostream
         * @param os out stream to send threevector data to
         * @param v threevector to be sent to ostream
         */
        friend std::ostream& operator<<(std::ostream& os,
            const threevector& v);
        ///@}
};

/** @brief Overload multiplication operator, scalar
 *
 *  return the scalar product of threevector v_old and scalar
 *  @param scalar the scalar multiple of the three vector
 *  @param v_old constant reference to threevector v_old
 *  @return threevector defines scalar product of scalar and v_old
 */
inline threevector operator*(double scalar, const threevector& v_old){
    return v_old*scalar;
};

/** @brief Overload multiplication operator, bool
 *
 *  return the scalar product of threevector v_old and bool
 *  @param scalar the scalar multiple of the three vector
 *  @param v_old constant reference to threevector v_old
 *  @return threevector defines scalar product of scalar and v_old
 */
inline threevector operator*(const bool scalar, const threevector& v_old){
    threevector v_new(v_old.getx()*scalar,v_old.gety()*scalar,
        v_old.getz()*scalar);
    return v_new;
};

#endif /* __THREEVECTOR_H_INCLUDED__ */
