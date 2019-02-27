#include "eigen_solver.h"
#include "eigen_symmetric.h"

EigenSolver::EigenSolver()
{

}


EigenSolver::~EigenSolver()
{

}



// compute the eigen values and eigen vectors of matrix m.
// after computation:
//  - the eigen values are stored in the member "m_eigen_values"
//  - the eigen vectors are stored in the member "m_eigen_vectors"

void EigenSolver::solve(const mat3& m) {
	// TODO -> call the function defined in "eigen_symmetric.h"
	// Please read carefully the manual of the function.
	
	int myvector[m.size()];
	int pos = 0;

	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix.size(); j++) {
			if (j > i) myvector[pos++] = matrix[i][j];
			else myvector[pos++] = matrix[j][i];
		}
	}


	eigen::eigen_symmetric(myvector, 3, m_eigen_vectors, m_eigen_values);
	//m_eigen_values = get_eigen_value(0);
	//m_eigen_vectors = get_eigen_vector(0);
}
