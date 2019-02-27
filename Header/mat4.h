#pragma once
#ifndef GEO1004_MAT_H
#define GEO1004_MAT_H

#include "vec4.h"

class mat4
{
public:
	// default constructor
	mat4();

	// initialized constructor
	mat4(
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33);

	// copy constructor
	mat4(const mat4& other);

	// destructor
	~mat4();

	// overwrite this mat4 with all zero entries
	void zeros();

	// overwrite this mat4 with identity
	void identity();

	// compute the transposed mat4
	mat4 transpose() const;

	// operators -- mat4-mat4
	mat4 operator+(const mat4& other) const; // addition
	mat4 operator-(const mat4& other) const; // subtraction
	mat4 operator*(const mat4& other) const; // multiplication
	mat4 operator-() const;			 // negation

	const mat4& operator+=(const mat4& other); // cumulative addition
	const mat4& operator-=(const mat4& other); // cumulative subtraction
	const mat4& operator*=(const mat4& other); // cumulative multiplication

	// operators -- mat4-vector
	vec4 operator*(const vec4& v);	// mat4-vector product

	// operators -- mat4-scalar
	const mat4& operator*=(float scalar);	// mat4-scalar product
	const mat4& operator/=(float scalar);	// mat4-scalar division
	mat4 operator*(float scalar);		// mat4-scalar product
	mat4 operator/(float scalar);		// mat4-scalar division

// assignment operator
	const mat4& operator=(const mat4& other);	// assignment

// access components
	float& operator()(int i, int j);	// RW access to components
	float operator()(int i, int j) const;	// RO access to components										// cast to float*

protected:
	float	m_data[16];
};

//DONE
inline mat4 operator*(float scalar, const mat4& M) {
	// TODO -- multiply each component of M with scalar, in a new mat4. return new mat4
	return mat4(M(0, 0) * scalar, M(0, 1) * scalar, M(0, 2) * scalar, M(0, 3)*scalar,
		M(1, 0) * scalar, M(1, 1) * scalar, M(1, 2) * scalar, M(1, 3)*scalar,
		M(2, 0) * scalar, M(2, 1) * scalar, M(2, 2) * scalar, M(2, 3)*scalar,
		M(3, 0) * scalar, M(3, 1) * scalar, M(3, 2) * scalar, M(3, 3)*scalar
	); // replace this line
}

//PREDONE
inline std::ostream& operator<<(std::ostream& out, const mat4& M) {
	// output a mat4 row-by-row to the "out" stream
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			out << M(i, j) << " ";
		}
		out << std::endl;
	}
	return out;
}

//DONE NOT SURE IF RIGHT
inline std::istream& operator>>(std::istream& in, mat4& M) {
	// TODO: read a mat4 row-by-row from the "in" stream
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			in >> M(i, j);
		}
	}
	return in;
}

//DONE
inline mat4::mat4() {
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
	m_data[9] = 0;
	m_data[10] = 0;
	m_data[11] = 0;
	m_data[12] = 0;
	m_data[13] = 0;
	m_data[14] = 0;
	m_data[15] = 0;
}

//DONE
inline mat4::mat4(float m00, float m01, float m02, float m03,
	float m10, float m11, float m12, float m13,
	float m20, float m21, float m22, float m23,
	float m30, float m31, float m32, float m33)
{
	// TODO -- initialize m_data with the provided components.
	m_data[0] = m00;
	m_data[1] = m01;
	m_data[2] = m02;
	m_data[3] = m03;
	m_data[4] = m10;
	m_data[5] = m11;
	m_data[6] = m12;
	m_data[7] = m13;
	m_data[8] = m20;
	m_data[9] = m21;
	m_data[10] = m22;
	m_data[11] = m23;
	m_data[12] = m30;
	m_data[13] = m31;
	m_data[14] = m32;
	m_data[15] = m33;

}

//DONE
inline mat4::mat4(const mat4& other) {
	// TODO -- copy other to (*this) component by component
	m_data[0] = other(0, 0);
	m_data[1] = other(0, 1);
	m_data[2] = other(0, 2);
	m_data[3] = other(0, 3);
	m_data[4] = other(1, 0);
	m_data[5] = other(1, 1);
	m_data[6] = other(1, 2);
	m_data[7] = other(1, 3);
	m_data[8] = other(2, 0);
	m_data[9] = other(2, 1);
	m_data[10] = other(2, 2);
	m_data[11] = other(2, 3);
	m_data[12] = other(3, 0);
	m_data[13] = other(3, 1);
	m_data[14] = other(3, 2);
	m_data[15] = other(3, 3);

}

//PREDONE
inline mat4::~mat4() {
}

//DONE
inline void mat4::zeros() {
	// TODO -- overwrite this mat4 with an empty matrix
	m_data[0] = 0;
	m_data[1] = 0;
	m_data[2] = 0;
	m_data[3] = 0;
	m_data[4] = 0;
	m_data[5] = 0;
	m_data[6] = 0;
	m_data[7] = 0;
	m_data[8] = 0;
	m_data[9] = 0;
	m_data[10] = 0;
	m_data[11] = 0;
	m_data[12] = 0;
	m_data[13] = 0;
	m_data[14] = 0;
	m_data[15] = 0;
}

//DONE
inline void mat4::identity() {
	// TODO -- overwrite this mat4 with an identity mat4.
	m_data[0] = 1;
	m_data[1] = 0;
	m_data[2] = 0;
	m_data[3] = 0;
	m_data[4] = 0;
	m_data[5] = 1;
	m_data[6] = 0;
	m_data[7] = 0;
	m_data[8] = 0;
	m_data[9] = 0;
	m_data[10] = 1;
	m_data[11] = 0;
	m_data[12] = 0;
	m_data[13] = 0;
	m_data[14] = 0;
	m_data[15] = 1;
}

//DONE
inline mat4 mat4::transpose() const {
	// TODO -- compute the transpose of this mat4 in a new mat4 and return.

	return mat4(m_data[0], m_data[4], m_data[8], m_data[12],
		m_data[1], m_data[5], m_data[9], m_data[13],
		m_data[2], m_data[6], m_data[10], m_data[14],
		m_data[3], m_data[7], m_data[11], m_data[15]); // replace this line
}

//DONE
inline mat4 mat4::operator+(const mat4& other) const {
	// TODO -- compute a new mat4 (*this)+other, return the new mat4

	return mat4(m_data[0] + other(0,0),
		m_data[1] + other(0,1), 
		m_data[2] + other(0,2), 
		m_data[3] + other(0,3), 
		m_data[4] + other(1,0), 
		m_data[5] + other(1,1), 
		m_data[6] + other(1,2), 
		m_data[7] + other(1,3), 
		m_data[8] + other(2,0), 
		m_data[9] + other(2,1), 
		m_data[10] + other(2,2), 
		m_data[11] + other(2,3), 
		m_data[12] + other(3,0), 
		m_data[13] + other(3,1),
		m_data[14] + other(3,2),
		m_data[15] + other(3,3)); // replace this line
}

//DONE
inline mat4 mat4::operator-(const mat4& other) const {
	// TODO -- compute a new mat4 (*this)-other, return the new mat4

	return mat4(m_data[0] - other(0, 0),
		m_data[1] - other(0, 1),
		m_data[2] - other(0, 2),
		m_data[3] - other(0, 3),
		m_data[4] - other(1, 0),
		m_data[5] - other(1, 1),
		m_data[6] - other(1, 2),
		m_data[7] - other(1, 3),
		m_data[8] - other(2, 0),
		m_data[9] - other(2, 1),
		m_data[10] - other(2, 2),
		m_data[11] - other(2, 3),
		m_data[12] - other(3, 0),
		m_data[13] - other(3, 1),
		m_data[14] - other(3, 2),
		m_data[15] - other(3, 3)); // replace this line
}

//DONE
inline mat4 mat4::operator*(const mat4& other) const {
	// TODO -- compute a new mat4 (*this) * other, return the new mat4

	return mat4(m_data[0] * other(0, 0) + m_data[1] * other(1, 0) + m_data[2] * other(2, 0) + m_data[3] * other(3, 0),
		m_data[0] * other(0, 1) + m_data[1] * other(1, 1) + m_data[2] * other(2, 1) + m_data[3] * other(3, 1),
		m_data[0] * other(0, 2) + m_data[1] * other(1, 2) + m_data[2] * other(2, 2) + m_data[3] * other(3, 2),
		m_data[0] * other(0, 3) + m_data[1] * other(1, 3) + m_data[2] * other(2, 3) + m_data[3] * other(3, 3),
		m_data[4] * other(0, 0) + m_data[5] * other(1, 0) + m_data[6] * other(2, 0) + m_data[7] * other(3, 0),
		m_data[4] * other(0, 1) + m_data[5] * other(1, 1) + m_data[6] * other(2, 1) + m_data[7] * other(3, 1),
		m_data[4] * other(0, 2) + m_data[5] * other(1, 2) + m_data[6] * other(2, 2) + m_data[7] * other(3, 2),
		m_data[4] * other(0, 3) + m_data[5] * other(1, 3) + m_data[6] * other(2, 3) + m_data[7] * other(3, 3),
		m_data[8] * other(0, 0) + m_data[9] * other(1, 0) + m_data[10] * other(2, 0) + m_data[11] * other(3, 0),
		m_data[8] * other(0, 1) + m_data[9] * other(1, 1) + m_data[10] * other(2, 1) + m_data[11] * other(3, 1),
		m_data[8] * other(0, 2) + m_data[9] * other(1, 2) + m_data[10] * other(2, 2) + m_data[11] * other(3, 2),
		m_data[8] * other(0, 3) + m_data[9] * other(1, 3) + m_data[10] * other(2, 3) + m_data[11] * other(3, 3),
		m_data[12] * other(0, 0) + m_data[13] * other(1, 0) + m_data[14] * other(2, 0) + m_data[15] * other(3, 0),
		m_data[12] * other(0, 1) + m_data[13] * other(1, 1) + m_data[14] * other(2, 1) + m_data[15] * other(3, 1),
		m_data[12] * other(0, 2) + m_data[13] * other(1, 2) + m_data[14] * other(2, 2) + m_data[15] * other(3, 2),
		m_data[12] * other(0, 3) + m_data[13] * other(1, 3) + m_data[14] * other(2, 3) + m_data[15] * other(3, 3)
	); // replace this line
}

//DONE
inline mat4 mat4::operator-() const {
	// TODO -- compute a new mat4 -(*this), return the new mat4

	return mat4(-m_data[0], -m_data[1], -m_data[2], -m_data[3],
		-m_data[4], -m_data[5], -m_data[6], -m_data[7],
		-m_data[8], -m_data[9], -m_data[10], -m_data[11],
		-m_data[12], -m_data[13], -m_data[14], -m_data[15]); // replace this line
}

//DONE
inline const mat4& mat4::operator+=(const mat4& other) {
	// TODO -- add other to this mat4
	m_data[0] = m_data[0] + other(0, 0);
	m_data[1] = m_data[1] + other(0, 1);
	m_data[2] = m_data[2] + other(0, 2);
	m_data[3] = m_data[3] + other(0, 3);
	m_data[4] = m_data[4] + other(1, 0);
	m_data[5] = m_data[5] + other(1, 1);
	m_data[6] = m_data[6] + other(1, 2);
	m_data[7] = m_data[7] + other(1, 3);
	m_data[8] = m_data[8] + other(2, 0);
	m_data[9] = m_data[9] + other(2, 1);
	m_data[10] = m_data[10] + other(2, 2);
	m_data[11] = m_data[11] + other(2, 3);
	m_data[12] = m_data[12] + other(3, 0);
	m_data[13] = m_data[13] + other(3, 1);
	m_data[14] = m_data[14] + other(3, 2);
	m_data[15] = m_data[15] + other(3, 3);
	return *this;
}

//DONE
inline const mat4& mat4::operator-=(const mat4& other) {
	// TODO -- subtract other from this mat4
	m_data[0] = m_data[0] - other(0, 0);
	m_data[1] = m_data[1] - other(0, 1);
	m_data[2] = m_data[2] - other(0, 2);
	m_data[3] = m_data[3] - other(0, 3);
	m_data[4] = m_data[4] - other(1, 0);
	m_data[5] = m_data[5] - other(1, 1);
	m_data[6] = m_data[6] - other(1, 2);
	m_data[7] = m_data[7] - other(1, 3);
	m_data[8] = m_data[8] - other(2, 0);
	m_data[9] = m_data[9] - other(2, 1);
	m_data[10] = m_data[10] - other(2, 2);
	m_data[11] = m_data[11] - other(2, 3);
	m_data[12] = m_data[12] - other(3, 0);
	m_data[13] = m_data[13] - other(3, 1);
	m_data[14] = m_data[14] - other(3, 2);
	m_data[15] = m_data[15] - other(3, 3);

	return *this;
}

//DONE
inline const mat4& mat4::operator*=(const mat4& other) {
	// TODO -- replace this mat4 by (*this) * other. Make sure you do not overwrite elements that you still need.
	//		   You may use mat4::operator*()
	float var00 = m_data[0] * other(0, 0) + m_data[1] * other(1, 0) + m_data[2] * other(2, 0) + m_data[3] * other(3, 0);
	float var01 = m_data[0] * other(0, 1) + m_data[1] * other(1, 1) + m_data[2] * other(2, 1) + m_data[3] * other(3, 1);
	float var02 = m_data[0] * other(0, 2) + m_data[1] * other(1, 2) + m_data[2] * other(2, 2) + m_data[3] * other(3, 2);
	float var03 = m_data[0] * other(0, 3) + m_data[1] * other(1, 3) + m_data[2] * other(2, 3) + m_data[3] * other(3, 3);
	float var10 = m_data[4] * other(0, 0) + m_data[5] * other(1, 0) + m_data[6] * other(2, 0) + m_data[7] * other(3, 0);
	float var11 = m_data[4] * other(0, 1) + m_data[5] * other(1, 1) + m_data[6] * other(2, 1) + m_data[7] * other(3, 1);
	float var12 = m_data[4] * other(0, 2) + m_data[5] * other(1, 2) + m_data[6] * other(2, 2) + m_data[7] * other(3, 2);
	float var13 = m_data[4] * other(0, 3) + m_data[5] * other(1, 3) + m_data[6] * other(2, 3) + m_data[7] * other(3, 3);
	float var20 = m_data[8] * other(0, 0) + m_data[9] * other(1, 0) + m_data[10] * other(2, 0) + m_data[11] * other(3, 0);
	float var21 = m_data[8] * other(0, 1) + m_data[9] * other(1, 1) + m_data[10] * other(2, 1) + m_data[11] * other(3, 1);
	float var22 =	m_data[8] * other(0, 2) + m_data[9] * other(1, 2) + m_data[10] * other(2, 2) + m_data[11] * other(3, 2);
	float var23 = m_data[8] * other(0, 3) + m_data[9] * other(1, 3) + m_data[10] * other(2, 3) + m_data[11] * other(3, 3);
	float var30 = m_data[12] * other(0, 0) + m_data[13] * other(1, 0) + m_data[14] * other(2, 0) + m_data[15] * other(3, 0);
	float var31 = m_data[12] * other(0, 1) + m_data[13] * other(1, 1) + m_data[14] * other(2, 1) + m_data[15] * other(3, 1);
	float var32 = m_data[12] * other(0, 2) + m_data[13] * other(1, 2) + m_data[14] * other(2, 2) + m_data[15] * other(3, 2);
	float var33 = m_data[12] * other(0, 3) + m_data[13] * other(1, 3) + m_data[14] * other(2, 3) + m_data[15] * other(3, 3);

	m_data[0] = var00;
	m_data[1] = var01;
	m_data[2] = var02;
	m_data[3] = var03;
	m_data[4] = var10;
	m_data[5] = var11;
	m_data[6] = var12;
	m_data[7] = var13;
	m_data[8] = var10;
	m_data[9] = var21;
	m_data[10] = var22;
	m_data[11] = var23;
	m_data[12] = var30;
	m_data[13] = var31;
	m_data[14] = var32;
	m_data[15] = var33;

	return *this;
}

//DONE
inline vec4 mat4::operator*(const vec4& v) {
	// TODO -- compute the mat4-vector product (*this) * v and return the result

	return vec4(v(0) * m_data[0] + v(1)* m_data[1] + v(2)*m_data[2] + v(3)*m_data[3],
		v(0) * m_data[4] + v(1)* m_data[5] + v(2)*m_data[6] + v(3)*m_data[7],
		v(0) * m_data[8] + v(1)* m_data[9] + v(2)*m_data[10] + v(3)*m_data[11],
		v(0) * m_data[12] + v(1)* m_data[13] + v(2)*m_data[14] + v(3)*m_data[15]); // replace this line
}

//DONE
inline const mat4& mat4::operator*=(float scalar) {
	// TODO -- multiply each mat4 component by scalar.
	m_data[0] = m_data[0] * scalar;
	m_data[1] = m_data[1] * scalar;
	m_data[2] = m_data[2] * scalar;
	m_data[3] = m_data[3] * scalar;
	m_data[4] = m_data[4] * scalar;
	m_data[5] = m_data[5] * scalar;
	m_data[6] = m_data[6] * scalar;
	m_data[7] = m_data[7] * scalar;
	m_data[8] = m_data[8] * scalar;
	m_data[9] = m_data[9] * scalar;
	m_data[10] = m_data[10] * scalar;
	m_data[11] = m_data[11] * scalar;
	m_data[12] = m_data[12] * scalar;
	m_data[13] = m_data[13] * scalar;
	m_data[14] = m_data[14] * scalar;
	m_data[15] = m_data[15] * scalar;

	return *this;
}

//DONE
inline const mat4& mat4::operator/=(float scalar) {
	assert("mat4::operator/= -- invalid argument" && scalar != 0);
	// TODO -- divide each mat4 component by scalar.
	m_data[0] = m_data[0] / scalar;
	m_data[1] = m_data[1] / scalar;
	m_data[2] = m_data[2] / scalar;
	m_data[3] = m_data[3] / scalar;
	m_data[4] = m_data[4] / scalar;
	m_data[5] = m_data[5] / scalar;
	m_data[6] = m_data[6] / scalar;
	m_data[7] = m_data[7] / scalar;
	m_data[8] = m_data[8] / scalar;
	m_data[9] = m_data[9] / scalar;
	m_data[10] = m_data[10] / scalar;
	m_data[11] = m_data[11] / scalar;
	m_data[12] = m_data[12] / scalar;
	m_data[13] = m_data[13] / scalar;
	m_data[14] = m_data[14] / scalar;
	m_data[15] = m_data[15] / scalar;

	return *this;
}

//DONE
inline mat4 mat4::operator*(float scalar) {
	// TODO -- compute a new mat4 (*this) * scalar.
	return mat4(m_data[0] * scalar,
		m_data[1] * scalar,
		m_data[2] * scalar,
		m_data[3] * scalar,
		m_data[4] * scalar,
		m_data[5] * scalar,
		m_data[6] * scalar,
		m_data[7] * scalar,
		m_data[8] * scalar,
		m_data[9] * scalar,
		m_data[10] * scalar,
		m_data[11] * scalar,
		m_data[12] * scalar,
		m_data[13] * scalar,
		m_data[14] * scalar,
		m_data[15] * scalar
	); // replace this line
}

//DONE
inline mat4 mat4::operator/(float scalar) {
	assert("mat4::operator/ -- invalid argument" && scalar != 0);
	// TODO -- divide each mat4 component by scalar and store in a new mat4. return the new mat4.
	return mat4(m_data[0] / scalar,
		m_data[1] / scalar,
		m_data[2] / scalar,
		m_data[3] / scalar,
		m_data[4] / scalar,
		m_data[5] / scalar,
		m_data[6] / scalar,
		m_data[7] / scalar,
		m_data[8] / scalar,
		m_data[9] / scalar,
		m_data[10] / scalar,
		m_data[11] / scalar,
		m_data[12] / scalar,
		m_data[13] / scalar,
		m_data[14] / scalar,
		m_data[15] / scalar
	); // replace this line
}

//DONE
inline const mat4& mat4::operator=(const mat4& other) {
	// TODO -- overwrite each component in this mat4 by the matching component in other
	m_data[0] = other(0, 0);
	m_data[1] = other(0, 1); 
	m_data[2] = other(0, 2); 
	m_data[3] = other(0, 3); 
	m_data[4] = other(1, 0); 
	m_data[5] = other(1, 1); 
	m_data[6] = other(1, 2); 
	m_data[7] = other(1, 3); 
	m_data[8] = other(2, 0); 
	m_data[9] = other(2, 1); 
	m_data[10] = other(2, 2); 
	m_data[11] = other(2, 3); 
	m_data[12] = other(3, 0); 
	m_data[13] = other(3, 1); 
	m_data[14] = other(3, 2); 
	m_data[15] = other(3, 3); 
	return *this;
}

//PREDONE
inline float& mat4::operator()(int i, int j) {
	assert("mat4::operator() -- invalid arguments" && i < 4 && j < 4);
	return m_data[i + 4 * j];
}

//PREDONE
inline float mat4::operator()(int i, int j) const {
	assert("mat4::operator() const -- invalid arguments" && i < 4 && j < 4);
	return m_data[i + j * 4];
}

#endif
