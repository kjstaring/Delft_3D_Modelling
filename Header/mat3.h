#pragma once
#ifndef GEO1004_MAT3_H
#define GEO1004_MAT3_H

#include "vec3.h"

class mat3
{
public:
	// default constructor
	mat3();

	// initialized constructor
	mat3(
		float m00, float m01, float m02,
		float m10, float m11, float m12,
		float m20, float m21, float m22);

	// copy constructor
	mat3(const mat3& other);

	// destructor
	~mat3();

	// overwrite this mat3 with all zero entries
	void zeros();

	// overwrite this mat3 with identity
	void identity();

	// compute the transposed mat3
	mat3 transpose() const;

	// operators -- mat3-mat3
	mat3 operator+(const mat3& other) const; // addition
	mat3 operator-(const mat3& other) const; // subtraction
	mat3 operator*(const mat3& other) const; // multiplication
	mat3 operator-() const;			 // negation

	const mat3& operator+=(const mat3& other); // cumulative addition
	const mat3& operator-=(const mat3& other); // cumulative subtraction
	const mat3& operator*=(const mat3& other); // cumulative multiplication

	// operators -- mat3-vector
	vec3 operator*(const vec3& v);	// mat3-vector product

	// operators -- mat3-scalar
	const mat3& operator*=(float scalar);	// mat3-scalar product
	const mat3& operator/=(float scalar);	// mat3-scalar division
	mat3 operator*(float scalar);		// mat3-scalar product
	mat3 operator/(float scalar);		// mat3-scalar division

// assignment operator
	const mat3& operator=(const mat3& other);	// assignment

// access components
	float& operator()(int i, int j);	// RW access to components
	float operator()(int i, int j) const;	// RO access to components										// cast to float*

protected:
	float	m_data[9];
};

//DONE
inline mat3 operator*(float scalar, const mat3& M) {
	// TODO -- multiply each component of M with scalar, in a new mat3. return new mat3

	return mat3(M(0,0)*scalar, M(0,1)*scalar, M(0,2)*scalar,
		M(1, 0)*scalar, M(1, 1)*scalar, M(1, 2)*scalar,
		M(2, 0)*scalar, M(2, 1)*scalar, M(2, 2)*scalar); // replace this line
}

//PREDONE
inline std::ostream& operator<<(std::ostream& out, const mat3& M) {
	// output a mat3 row-by-row to the "out" stream
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			out << M(i, j) << " ";
		}
		out << std::endl;
	}
	return out;
}

//DONE NOT SURE IF RIGHT
inline std::istream& operator>>(std::istream& in, mat3& M) {
	// TODO: read a mat3 row-by-row from the "in" stream
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			in >> M(i, j);
		}
	}
	return in;
}

//DONE
inline mat3::mat3() {
	// TODO -- initialize m_data with 0s
	m_data[0] = 0;
	m_data[1] = 0;
	m_data[2] = 0;
	m_data[3] = 0;
	m_data[4] = 0;
	m_data[5] = 0;
	m_data[6] = 0;
	m_data[7] = 0;
	m_data[8] = 0;
}

//DONE
inline mat3::mat3(float m00, float m01, float m02,
	float m10, float m11, float m12,
	float m20, float m21, float m22)
{
	// TODO -- initialize m_data with the provided components.
	m_data[0] = m00;
	m_data[1] = m01;
	m_data[2] = m02;
	m_data[3] = m10;
	m_data[4] = m11;
	m_data[5] = m12;
	m_data[6] = m20;
	m_data[7] = m21;
	m_data[8] = m22;
}

//DONE
inline mat3::mat3(const mat3& other) {
	// TODO -- copy other to (*this) component by component
	m_data[0] = other(0, 0);
	m_data[1] = other(0, 1);
	m_data[2] = other(0, 2);
	m_data[3] = other(1, 0);
	m_data[4] = other(1, 1);
	m_data[5] = other(1, 2);
	m_data[6] = other(2, 0);
	m_data[7] = other(2, 1);
	m_data[8] = other(2, 2);
}

//PREDONE
inline mat3::~mat3() {
}

//DONE
inline void mat3::zeros() {
	//overwrite this mat3 with all zero entries
	m_data[0] = 0;
	m_data[1] = 0;
	m_data[2] = 0;
	m_data[3] = 0;
	m_data[4] = 0;
	m_data[5] = 0;
	m_data[6] = 0;
	m_data[7] = 0;
	m_data[8] = 0;
}

//DONE
inline void mat3::identity() {
	// TODO -- overwrite this mat3 with an identity mat3.
	m_data[0] = 1;
	m_data[1] = 0;
	m_data[2] = 0;
	m_data[3] = 0;
	m_data[4] = 1;
	m_data[5] = 0;
	m_data[6] = 0;
	m_data[7] = 0;
	m_data[8] = 1;
}

//DONE
inline mat3 mat3::transpose() const {
	// TODO -- compute the transpose of this mat3 in a new mat3 and return.
	
	return mat3(m_data[0], m_data[3], m_data[6], 
		m_data[1], m_data[4], m_data[7], 
		m_data[2], m_data[5], m_data[8]); // replace this line
}

//DONE
inline mat3 mat3::operator+(const mat3& other) const {
	// TODO -- compute a new mat3 (*this)+other, return the new mat3
	return mat3(m_data[0] + other(0,0), m_data[1] + other(0,1), m_data[2] + other(0,2), 
		m_data[3] + other(1,0), m_data[4] + other(1,1), m_data[5] + other(1,2), 
		m_data[6] + other(2,0), m_data[7] + other(2,1), m_data[8] + other(2,2)); // replace this line
}

//DONE
inline mat3 mat3::operator-(const mat3& other) const {
	// TODO -- compute a new mat3 (*this)-other, return the new mat3
	return mat3(m_data[0] - other(0, 0), m_data[1] - other(0, 1), m_data[2] - other(0, 2),
		m_data[3] - other(1, 0), m_data[4] - other(1, 1), m_data[5] - other(1, 2),
		m_data[6] - other(2, 0), m_data[7] - other(2, 1), m_data[8] - other(2, 2)); // replace this line
}

//DONE
inline mat3 mat3::operator*(const mat3& other) const {
	// TODO -- compute a new mat3 (*this) * other, return the new mat3
	
	return mat3(m_data[0] * other(0, 0) + m_data[1] * other(1, 0) + m_data[2] * other(2, 0),
		m_data[0] * other(0, 1) + m_data[2] * other(1, 1) + m_data[2] * other(2, 1),
		m_data[0] * other(0, 2) + m_data[2] * other(1, 2) + m_data[2] * other(2, 2),
		m_data[3] * other(0, 0) + m_data[4] * other(1, 0) + m_data[5] * other(2, 0),
		m_data[3] * other(0, 1) + m_data[4] * other(1, 1) + m_data[5] * other(2, 1),
		m_data[3] * other(0, 2) + m_data[4] * other(1, 2) + m_data[5] * other(2, 2),
		m_data[6] * other(0, 1) + m_data[7] * other(1, 1) + m_data[8] * other(2, 1),
		m_data[6] * other(0, 2) + m_data[7] * other(1, 2) + m_data[8] * other(2, 2),
		m_data[6] * other(0, 2) + m_data[7] * other(1, 2) + m_data[8] * other(2, 2)
	); // replace this line
}

//DONE
inline mat3 mat3::operator-() const {
	// TODO -- compute a new mat3 -(*this), return the new mat3
	return mat3(-m_data[0], -m_data[1], -m_data[2],
		-m_data[3], -m_data[4], -m_data[5],
		-m_data[6], -m_data[7], -m_data[8]); // replace this line
}

//DONE
inline const mat3& mat3::operator+=(const mat3& other) {
	// TODO -- add other to this mat3
	m_data[0] = m_data[0] - other(0, 0);
	m_data[1] = m_data[1] - other(0, 1);
	m_data[2] = m_data[2] - other(0, 2);
	m_data[3] = m_data[3] - other(1, 0);
	m_data[4] = m_data[4] - other(1, 1);
	m_data[5] = m_data[5] - other(1, 2);
	m_data[6] = m_data[6] - other(2, 0);
	m_data[7] = m_data[7] - other(2, 1);
	m_data[8] = m_data[8] - other(2, 2);
	return *this;
}

//DONE
inline const mat3& mat3::operator-=(const mat3& other) {
	// TODO -- subtract other from this mat3
	m_data[0] = m_data[0] - other(0, 0);
	m_data[1] = m_data[1] - other(0, 1);
	m_data[2] = m_data[2] - other(0, 2);
	m_data[3] = m_data[3] - other(1, 0);
	m_data[4] = m_data[4] - other(1, 1);
	m_data[5] = m_data[5] - other(1, 2);
	m_data[6] = m_data[6] - other(2, 0);
	m_data[7] = m_data[7] - other(2, 1);
	m_data[8] = m_data[8] - other(2, 2);
	return *this;
}

//DONE
inline const mat3& mat3::operator*=(const mat3& other) {
	// TODO -- replace this mat3 by (*this) * other. Make sure you do not overwrite elements that you still need.
	//		   You may use mat3::operator*()
	
	float var00 = m_data[0] * other(0, 0) + m_data[1] * other(1, 0) + m_data[2] * other(2, 0);
	float var01 = m_data[0] * other(0, 1) + m_data[2] * other(1, 1) + m_data[2] * other(2, 1);
	float var02 = m_data[0] * other(0, 2) + m_data[2] * other(1, 2) + m_data[2] * other(2, 2);
	float var10 = m_data[3] * other(0, 0) + m_data[4] * other(1, 0) + m_data[5] * other(2, 0);
	float var11 = m_data[3] * other(0, 1) + m_data[4] * other(1, 1) + m_data[5] * other(2, 1);
	float var12 = m_data[3] * other(0, 2) + m_data[4] * other(1, 2) + m_data[5] * other(2, 2);
	float var20 = m_data[6] * other(0, 1) + m_data[7] * other(1, 1) + m_data[8] * other(2, 1);
	float var21 = m_data[6] * other(0, 2) + m_data[7] * other(1, 2) + m_data[8] * other(2, 2);
	float var22 = m_data[6] * other(0, 2) + m_data[7] * other(1, 2) + m_data[8] * other(2, 2);

	m_data[0] = var00;
	m_data[1] = var01;
	m_data[2] = var02;
	m_data[3] = var10;
	m_data[4] = var11;
	m_data[5] = var12;
	m_data[6] = var20;
	m_data[7] = var21;
	m_data[8] = var22;

	return *this;
}

//DONE
inline vec3 mat3::operator*(const vec3& v) {
	// TODO -- compute the mat3-vector product (*this) * v and return the result

	return vec3(v(0) * m_data[0] + v(1)* m_data[1] + v(2)*m_data[2],
		v(0) * m_data[3] + v(1)* m_data[4] + v(2)*m_data[5],
		v(0) * m_data[6] + v(1)* m_data[7] + v(2)*m_data[8]); // replace this line
}

//DONE
inline const mat3& mat3::operator*=(float scalar) {
	// TODO -- multiply each mat3 component by scalar.
	m_data[0] = m_data[0] * scalar;
	m_data[1] = m_data[1] * scalar;
	m_data[2] = m_data[2] * scalar;
	m_data[3] = m_data[3] * scalar;
	m_data[4] = m_data[4] * scalar;
	m_data[5] = m_data[5] * scalar;
	m_data[6] = m_data[6] * scalar;
	m_data[7] = m_data[7] * scalar;
	m_data[8] = m_data[8] * scalar;
	return *this;
}

//DONE
inline const mat3& mat3::operator/=(float scalar) {
	assert("mat3::operator/= -- invalid argument" && scalar != 0);
	// TODO -- divide each mat3 component by scalar.
	m_data[0] = m_data[0] / scalar;
	m_data[1] = m_data[1] / scalar;
	m_data[2] = m_data[2] / scalar;
	m_data[3] = m_data[3] / scalar;
	m_data[4] = m_data[4] / scalar;
	m_data[5] = m_data[5] / scalar;
	m_data[6] = m_data[6] / scalar;
	m_data[7] = m_data[7] / scalar;
	m_data[8] = m_data[8] / scalar;

	return *this;
}

//DONE
inline mat3 mat3::operator*(float scalar) {
	// TODO -- compute a new mat3 (*this) * scalar.
	return mat3(m_data[0] * scalar, m_data[1] * scalar, m_data[2] * scalar,
		m_data[3] * scalar, m_data[4] * scalar, m_data[5] * scalar,
		m_data[6] * scalar, m_data[7] * scalar, m_data[8] * scalar); // replace this line
}

//DONE
inline mat3 mat3::operator/(float scalar) {
	assert("mat3::operator/ -- invalid argument" && scalar != 0);
	// TODO -- divide each mat3 component by scalar and store in a new mat3. return the new mat3.
	return mat3(m_data[0] / scalar, m_data[1] / scalar, m_data[2] / scalar,
		m_data[3] / scalar, m_data[4] / scalar, m_data[5] / scalar,
		m_data[6] / scalar, m_data[7] / scalar, m_data[8] / scalar); // replace this line
}

//DONE
inline const mat3& mat3::operator=(const mat3& other) {
	// TODO -- overwrite each component in this mat3 by the matching component in other
	m_data[0] = other(0, 0);
	m_data[1] = other(0, 1);
	m_data[2] = other(0, 2);
	m_data[3] = other(1, 0);
	m_data[4] = other(1, 1);
	m_data[5] = other(1, 2);
	m_data[6] = other(2, 0);
	m_data[7] = other(2, 1);
	m_data[8] = other(2, 2);
	return *this;
}

//PREDONE
inline float& mat3::operator()(int i, int j) {
	assert("mat3::operator() -- invalid arguments" && i < 4 && j < 4);
	return m_data[i + 4 * j];
}

//PREDONE
inline float mat3::operator()(int i, int j) const {
	assert("mat3::operator() const -- invalid arguments" && i < 4 && j < 4);
	return m_data[i + j * 4];
}



#endif
