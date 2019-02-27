
#include "vec3.h"
#include "vec4.h"
#include "mat3.h"
#include "mat4.h"

#include <iostream>

int main(int argc, char** argv) {
	std::cout << "Hello!!!" << std::endl << std::endl;

	 /*
	 TODO -- test every vector and matrix function. Submit your code with tests in here.
	 * For example:
	 *
	 * mat4 M;
	 * std::cout << "the matrix is: " << std::endl << M << std::endl;
	 *
	 * vec4 v(1.2f, 2.3f, -1.7f);
	 * std::cout << "the vector is: " << std::endl << v << std::endl;
	 *
	 * std::cout << "M * p: " << std::endl << M * v << std::endl;
	 */

	/////////
	//VEC 3//
	/////////

	//Constructors
	vec3 v3_1 = vec3();
	vec3 v3_2 = vec3(1, 2, 3);
	vec3 v3_3 = vec3(v3_1);
	vec3 v3_4 = vec3(v3_2);

	//Deconstructors
	v3_3.~vec3();

	//Functions
	//squared length
	float sqlv3_1 = v3_1.length2();
	float sqlv3_2 = v3_2.length2();

	//length
	float lv3_1 = v3_1.length();
	float lv3_2 = v3_2.length();

	//normalize
	v3_1.normalize();
	v3_2.normalize();

	//dotproduct
	float dp3_1 = v3_1.dot(v3_2);
	float dp3_2 = v3_2.dot(v3_4);

	//crossproduct
	vec3 v3_cp1 = v3_1.cross(v3_2);
	vec3 v3_cp2 = v3_1.cross(v3_4);

	//operators vector vector
	vec3 v3_op1 = vec3(2, 2, 2);
	vec3 v3_op2 = vec3(1, 1, 1);

	v3_op1 += v3_op2;
	v3_op1 -= v3_op2;
	v3_op1 = v3_op1 - v3_op2;
	v3_op1 = v3_op1 + v3_op1;
	v3_op1 = -v3_op1;

	//operators vector scalar
	v3_op1 *= 2;
	v3_op1 /= 2;
	v3_op1 = v3_op1 * 2;
	v3_op1 = v3_op1 / 2;
	v3_op1 = v3_op2;

	/////////
	//VEC 4//
	/////////

		//Constructors
	vec4 v4_1 = vec4();
	vec4 v4_2 = vec4(1, 2, 3);
	vec4 v4_3 = vec4(v3_1);
	vec4 v4_4 = vec4(v3_2);
	vec4 v4_5 = vec4(1, 2, 3, 4);
	vec4 v4_6 = vec4(v3_1);

	//Deconstructors
	v4_3.~vec4();

	//Functions
	//squared length
	float sqlv4_1 = v4_1.length2();
	float sqlv4_2 = v4_2.length2();

	//length
	float lv4_1 = v4_1.length();
	float lv4_2 = v4_2.length();

	//normalize
	v4_1.normalize();
	v4_2.normalize();

	//dotproduct
	float dp4_1 = v4_1.dot(v4_2);
	float dp4_2 = v4_2.dot(v4_4);

	//crossproduct
	vec4 v4_cp1 = v4_1.cross(v4_2);
	vec4 v4_cp2 = v4_1.cross(v4_4);

	//operators vector vector
	vec4 v4_op1 = vec4(2, 2, 2, 2);
	vec4 v4_op2 = vec4(1, 1, 1, 1);

	v4_op1 += v4_op2;
	v4_op1 -= v4_op2;
	v4_op1 = v4_op1 - v4_op2;
	v4_op1 = v4_op1 + v4_op1;
	v4_op1 = -v4_op1;

	//operators vector scalar
	v4_op1 *= 2;
	v4_op1 /= 2;
	v4_op1 = v4_op1 * 2;
	v4_op1 = v4_op1 / 2;
	v4_op1 = v4_op2;

	return 0;
}