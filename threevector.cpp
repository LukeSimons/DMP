#include "threevector.h"
#include "Constants.h"

// Default constructor
threevector::threevector() 
{
	xcoord = 0.0; ycoord = 0.0; zcoord = 0.0;
}

// Cartesian constructor
threevector::threevector(double x, double y, double z)
{
	xcoord = x; ycoord = y; zcoord = z;
}

// Polar constructor
threevector::threevector(double r, double theta, double phiorz, char type)
{
	// Spherical
	if(type=='s')
	{
		xcoord = r*sin(theta)*cos(phiorz);
		ycoord = r*sin(theta)*sin(phiorz);
		zcoord = r*cos(theta);
	}

	// Cylindrical
	if(type=='c')
	{
		xcoord = r*cos(theta);
		ycoord = r*sin(theta);
		zcoord = phiorz;
	}
	else std::cout << "Error: Unknown type of polar threevector" << std::endl;
}

// Method to print out contents to screen
void threevector::print()
{
	std::cout << "(" << xcoord << ", " << ycoord << ", " << zcoord << ")" << std::endl;
}

// Method to print out contents to file
void threevector::print(std::ofstream &fout)
{
	fout << "(" << xcoord << ", " << ycoord << ", " << zcoord << ")" << '\t';
}

// Modifier method to change a component of a threevector
void threevector::setx(double xvalue)
{
	xcoord = xvalue;
}
void threevector::sety(double yvalue)
{
	ycoord = yvalue;
}
void threevector::setz(double zvalue)
{
	zcoord = zvalue;
}

// Access method for coordinates
double threevector::getx()const
{
	return xcoord;
}
double threevector::gety()const
{
	return ycoord;
}
double threevector::getz()const
{
	return zcoord;
}

// Calculating the square of a threevector
double threevector::square()const
{
	return xcoord*xcoord + ycoord*ycoord + zcoord*zcoord;
}

// Calculating the magnitude
double threevector::mag3()const
{
	return sqrt(square());
}

// Calculate theta
double threevector::gettheta()const
{
	return acos(zcoord/mag3());
}

// Calculate phi in range 0 to 2pi
double threevector::getphi()const
{
	double phi = atan2(ycoord,xcoord);
	if(phi < 0.0) phi += 2.0*PI; // atan2 range is -pi to +pi
	return phi;
}

// Calculate a unit threevector
threevector threevector::getunit()const
{
	threevector unitvec(0.0,0.0,0.0);
	if(mag3()!=0.0){
		unitvec.setx(xcoord/mag3());
		unitvec.sety(ycoord/mag3());
		unitvec.setz(zcoord/mag3());
	}
	return unitvec;
}

// Overload the + operator
threevector threevector::operator+(const threevector v_old)const
{
	threevector v_sum(v_old.getx() + xcoord, v_old.gety() + ycoord, v_old.getz() + zcoord);
	return v_sum;
}

// Overload the - operator
threevector threevector::operator-(const threevector v_old)const
{
	threevector v_sum(xcoord - v_old.getx(), ycoord - v_old.gety(), zcoord - v_old.getz());
	return v_sum;
}

// Overload the += operator
void threevector::operator+=(threevector v_old)
{
	xcoord += v_old.getx();
	ycoord += v_old.gety();
	zcoord += v_old.getz();
}

// Overload the -= operator
void threevector::operator-=(threevector v_old)
{
	xcoord -= v_old.getx();
	ycoord -= v_old.gety();
	zcoord -= v_old.getz();
}

// Overload the = operator
void threevector::operator=(threevector v_old)
{
	xcoord = v_old.getx();
	ycoord = v_old.gety();
	zcoord = v_old.getz();
}

// Overload the * operator to multiply by a scalar
threevector threevector::operator*(double scalar)const
{
	threevector v_new(xcoord*scalar,ycoord*scalar,zcoord*scalar);
	return v_new;
}

// The dot product
double threevector::operator*(threevector v_old)const
{
	return xcoord*v_old.getx() + ycoord*v_old.gety() +  zcoord*v_old.getz();
}

// The threevector product
threevector threevector::operator^(threevector v_old)const
{
	threevector v_new((ycoord * v_old.getz()) - (zcoord * v_old.gety()),
				 (zcoord * v_old.getx()) - (xcoord * v_old.getz()),
				 (xcoord * v_old.gety()) - (ycoord * v_old.getx()));
	return v_new;
}

// The threevector product
std::ostream& operator<<(std::ostream& os, const threevector &v)
{
	os << v.getx() << " " << v.gety() << " " << v.getz();
	return os;
}

threevector operator*(double scalar, const threevector& y){
	return y*scalar;
}

threevector operator*(const bool scalar, const threevector& y){
	threevector v_new(y.getx()*scalar,y.gety()*scalar,y.getz()*scalar);
	return v_new;
}

