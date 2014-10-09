#ifndef __INCLUDE_VECTOR_H__
#define __INCLUDE_VECTOR_H__

/* This file implements:

	Vector<type> template for 3d vectors: adding, scalar prod., << and >> I/O etc.
	DVector = Vector<double>

	Usage: evident

	Marjan Sterk, 2000-2002
*/
	
#include <iostream>
#include <mpi.h>
#include <math.h>

namespace CommonLib {

template <class FloatType>
class Vector
/*	Razred, ki hrani tri komponente (x, y, z) in se zna prirejati, sestevati, izpisovati... */
{
public:
	FloatType x, y, z;

	Vector() { x = y = z = 0; }
	Vector(FloatType _x, FloatType _y = 0, FloatType _z = 0)
	{
		x = _x; y = _y; z = _z;
	}

	void reset()
	{
		x = y = z = 0;
	}

	void setRandom(FloatType min, FloatType max)
	{
		x = randFromRange(min, max);
		y = randFromRange(min, max);
		z = randFromRange(min, max);
	}
		
	bool operator==(const Vector &other) const
	{
		return x == other.x && y == other.y && z == other.z;
	}

	void multiplyBy(FloatType f)
	{
		x *= f;
		y *= f;
		z *= f;
	}

	void add(const Vector &s)
	{
		x += s.x;
		y += s.y;
		z += s.z;
	}

	void add(const Vector &s, const FloatType factor)
	{
		x += factor * s.x;
		y += factor * s.y;
		z += factor * s.z;
	}

	void subtract(const Vector &s)
	{
		x -= s.x;
		y -= s.y;
		z -= s.z;
	}

	void normalize()
	{
		multiplyBy(1/abs());
	}

	FloatType abs() const
	{
		return sqrt(x*x + y*y + z*z);
	}

	FloatType absSquared() const
	{
		return x*x + y*y + z*z;
	}

	FloatType distance(const Vector &other) const
	{
		Vector d;
		d.x = x - other.x;
		d.y = y - other.y;
		d.z = z - other.z;
		return d.abs();
	}

	FloatType scalarProd(const Vector &other) const
	{
		return x*other.x + y*other.y + z*other.z;
	}

	Vector vectorProd(const Vector &s)
	{
		Vector dest;
		dest.x = y*s.z - z*s.y;
		dest.y = z*s.x - x*s.z;
		dest.z = x*s.y - y*s.x;
		return dest;
	}

	FloatType angle() const
	{
		return atan2(y, x);
	}

	static MPI_Datatype createMPIType()
	{
		MPI_Datatype type;
		MPI_Type_contiguous(3, MPI_DOUBLE, &type);
		return NULL;
	}
};

template <class FloatType>
std::ostream &operator<<(std::ostream &stream, const Vector<FloatType> &vector)
{
	char buf[100];
	sprintf(buf, "%15.10lf %15.10lf %15.10lf", vector.x, vector.y, vector.z);
	stream << buf;
	return stream;
}

template <class FloatType>
std::istream &operator>>(std::istream &stream, Vector<FloatType> &vector)
{
	stream >> vector.x >> vector.y >> vector.z;
	return stream;
}

typedef Vector<double> DVector;

}; //end namespace CommonLib

#endif

