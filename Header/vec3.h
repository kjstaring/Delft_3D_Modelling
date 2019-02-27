#pragma once
#pragma once

#ifndef GEO1004_VEC3_H
#define GEO1004_VEC3_H	// TODO -- why is this here? Be prepared to explain to the TA.

#include<iostream>
#include<cfloat>
#include<cmath>
#include<cassert>

/// three-dimensional vector class

class vec3 {
public:
	// default constructor: initialize all elements to be zero
	vec3();

	// initialized constructor
	vec3(float x, float y, float z);

	// copy constructor
	vec3(const vec3& other);

	// destructor
	~vec3();

	// sqared length of the vector
	float length2() const;

	// length of the vector
	float length() const;

	// normalize to a unit vector
	void normalize();

	// compute the dotproduct between this and other
	float dot(const vec3& other) const;

	// compute the crossproduct between this and other
	vec3 cross(const vec3& other) const;

	// operators -- vector-vector
	const vec3& operator+=(const vec3& other);		// cumulative addition
	const vec3& operator-=(const vec3& other);		// cumulative subtraction
	vec3 operator-() const;							// negation (unary operator)
	vec3 operator+(const vec3& other) const;		// addition
	vec3 operator-(const vec3& other) const;		// subtraction

	// operators -- vector-scalar
	const vec3& operator*=(float scalar);			// multiplication by scalar
	const vec3& operator/=(float scalar);			// division by scalar
	vec3 operator*(float scalar);					// vector times scalar
	vec3 operator/(float scalar);					// vector divided by scalar
	const vec3 operator=(const vec3& other);		// assignment operator

	// element access
	float& operator()(int index);       // RW access to element: return the index_th element
	float operator()(int index) const;	// RO access to element: return the index_th element

protected:
	float	m_data[3];	// data array
};

//DONE
// TODO -- why is the following needed? Be prepared to explain to the TA
// It is called overloading operators, it is needed because C++ cannot handle multiplication of vectors
// only if we describe, what exactly needs to be done when calling "*", C++ can handle this case

//DONE
inline vec3 operator*(float scalar, const vec3& v) {
	// TODO -- multiply each component of v with scalar, in a new vector. return new vector
	return vec3(v(0) * scalar, v(1) * scalar, v(2)*scalar);
}

//PREDONE
inline std::ostream& operator<<(std::ostream& out, const vec3& v) {
	for (int i = 0; i < 3; ++i)
		out << v(i) << " ";
	return out;
}

//DONE NOT SURE IF RIGHT
inline std::istream& operator>>(std::istream& in, vec3& v) {
	// TODO: read a vector component-wise from the "in" stream
	for (int i = 0; i < 3; ++i)
		in >> v(i);
	return in;
}

//DONE
inline vec3::vec3() {
	// TODO -- initialize m_data with 0s
	m_data[0] = 0;
	m_data[1] = 0;
	m_data[2] = 0;
}

//DONE
inline vec3::vec3(float x, float y, float z) {
	// TODO -- initialize m_data with x,y,z What does the "= float(0)" in the class definition do? Be prepared to explain to the TA.
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
}

//DONE
inline vec3::vec3(const vec3& other) {
	// TODO -- copy contents of other to this vector
	m_data[0] = other(0);
	m_data[1] = other(1);
	m_data[2] = other(2);
}

//DONE
inline vec3::~vec3() {
	// TODO -- is there anything to do here? Be prepared to explain to the TA.
	//Nope, nothing to do here, the deconstructor take no arguments and have no returns, so nothing to process
}

//DONE
inline float vec3::length2() const {
	// TODO -- compute squared length of the vector
	return float(length()*length()); // replace this line
}

//DONE
inline float vec3::length() const {
	// TODO -- compute the length of the vector
	return float(sqrt(m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2]));
}

//DONE
inline void vec3::normalize() {
	// TODO -- if length()==0, do nothing, otherwise normalize the vector
	if (length()!=0)
	{
		m_data[0] = m_data[0] / length();
		m_data[1] = m_data[1] / length();
		m_data[2] = m_data[2] / length();
	}
}

//DONE
inline float vec3::dot(const vec3& other) const {
	// TODO -- compute dot product between this vector and other
	return float(m_data[0] * other(0) + m_data[1] * other(1) + m_data[2] * other(2)); // replace this line.
}

//DONE
inline vec3 vec3::cross(const vec3& other) const {
	// TODO -- compute the crossproduct between the first three components of this and other
	return vec3(m_data[1] * other(2) - m_data[2] * other(1), m_data[2] * other(0) - m_data[0] * other(2), m_data[0] * other(1) - m_data[1] * other(0)); // replace this line
}

//DONE - but answer for question must be improved
inline const vec3& vec3::operator+=(const vec3& other) {
	// TODO -- add other to this vector component-wise, store in this vector
	
	m_data[0] = m_data[0] + other(0);
	m_data[1] = m_data[1] + other(1);
	m_data[2] = m_data[2] + other(2);

	return *this; // TODO -- why would you return a reference to *this? Be prepared to explain to the TA.
	//It returns a temporary copy of the object

}

//DONE
inline const vec3& vec3::operator-=(const vec3& other) {
	// TODO -- subtract other from this vector component-wise, store in this vector

	m_data[0] = m_data[0] - other(0);
	m_data[1] = m_data[1] - other(1);
	m_data[2] = m_data[2] - other(2);

	return *this;
}

//DONE - but answer for question must be improved
inline vec3 vec3::operator-() const {
	// TODO -- why can't we return a reference, here? Be prepared to explain to the TA.
	//Because we're not using the equal sign here

	// TODO -- return a new vector in which each component is the negated component from this vector
	return vec3(-m_data[0], -m_data[1], -m_data[2]); // replace this line
}

//DONE
inline vec3 vec3::operator+(const vec3& other) const {
	// TODO -- return a new vector in which each component is the sum of the components of this vector and other

	return vec3(m_data[0] + other(0), m_data[1] + other(1), m_data[2] + other(2)); // replace this line
}

//DONE
inline vec3 vec3::operator-(const vec3& other) const {
	// TODO -- return a new vector in which each component is the difference between the components of this vector and other		
	return vec3(abs(m_data[0] - other(0)), abs(m_data[1] - other(1)), abs(m_data[2] - other(2))); // replace this line
}

//DONE
inline const vec3& vec3::operator*=(float scalar) {
	// TODO -- replace each component of this with the matching component of this multiplied with the scalar
	//         Make sure to convert the scalar from S to float
	m_data[0] = m_data[0] * scalar;
	m_data[1] = m_data[1] * scalar;
	m_data[2] = m_data[2] * scalar;
	return *this;
}

//DONE
inline const vec3& vec3::operator/=(float scalar) {
	// TODO -- why do we convert scalar and 0 to type float and not compare it as type S? Be prepared to explain to the TA.
	//Because the values in m_data are also float
	assert("vec3::operator/= -- invalid argument" && float(scalar) != float(0));
	// TODO -- replace each component of this with the matching component of this multiplied with the scalar
	//         Make sure to convert the scalar from S to float

	m_data[0] = m_data[0] / scalar;
	m_data[1] = m_data[1] / scalar;
	m_data[2] = m_data[2] / scalar;

	return *this;
}

//DONE
inline vec3 vec3::operator*(float scalar) {
	// TODO -- return a new vector in which each component equals the matching component of this vector times the scalar
	//		   Make sure to convert the scalar from S to float
	
	return vec3(m_data[0] * scalar, m_data[1] * scalar, m_data[2] * scalar); // replace this line
}

//DONE
inline vec3 vec3::operator/(float scalar) {
	assert("vec3::operator/ -- invalid argument" && float(scalar) != float(0));
	// TODO -- return a new vector in which each component equals the matching component of this vector divided by the scalar

	return vec3(m_data[0] / scalar, m_data[1] / scalar, m_data[2] / scalar); // replace this line
}

//DONE
inline const vec3 vec3::operator=(const vec3& other) {
	// TODO -- overwrite each component in this vector with the matching component of other.
	m_data[0] = other(0);
	m_data[1] = other(1);
	m_data[2] = other(2);
	return *this;
}

//PREDONE
inline float& vec3::operator()(int n) {
	assert("vec3::operator() -- invalid argument" && n < 4);
	return m_data[n];
}

//PREDONE
inline float vec3::operator()(int n) const {
	assert("vec3::operator() const -- invalid argument" && n < 4);
	return m_data[n];
}

#endif
