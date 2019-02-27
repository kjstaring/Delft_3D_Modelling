#pragma once
#ifndef EIGEN_SEMI_DEFINITE_SYMMEATRIC_H
#define EIGEN_SEMI_DEFINITE_SYMMEATRIC_H


namespace eigen {

	/**
	* computes the eigen values and eigen vectors of a semi definite symmetric matrix.
	*
	* @param  matrix is stored in column symmetric storage (only the upper triangle),
	*	  i.e. matrix = { m11, m12, m22, m13, m23, m33, m14, m24, m34, m44 ... },
	*     size = n(n+1)/2
	*
	* @param eigen_values are in decreasing order. size = n, must be allocated by caller.
	*
	* @param eigen_vectors = { v1, v2, v3, ..., vn } are sorted in accordance.
	*     where vk = vk0, vk1, ..., vkn ( 1 <= k <= n).
	*     size = n^2, must be allocated by caller.
	*/
	template <class FT>
	static inline void eigen_symmetric(const FT *matrix,
		const int n,
		FT *eigen_vectors,
		FT *eigen_values,
		const int MAX_ITER = 100,
		const FT EPSILON = (FT)1e-016)
	{
		// number of entries in mat
		int nn = (n*(n + 1)) / 2;

		// Step 1: copy matrix into a
		FT *a = new FT[nn];
		int ij;
		for (ij = 0; ij < nn; ij++)
			a[ij] = matrix[ij];
		// Ugly Fortran-porting trick: indices for a are between 1 and n
		a--;

		// Step 2: init diagonalization matrix as the unit matrix
		FT *v = new FT[n*n];
		ij = 0;
		int i;
		for (i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				if (i == j)
					v[ij++] = 1.0;
				else
					v[ij++] = 0.0;
		// Ugly Fortran-porting trick: indices for a are between 1 and n
		v--;

		// Step 3: compute weight of the non diagonal terms 
		ij = 1;
		FT a_norm = 0.0;
		for (i = 1; i <= n; i++)
			for (int j = 1; j <= i; j++)
			{
				if (i != j)
				{
					FT a_ij = a[ij];
					a_norm += a_ij * a_ij;
				}
				ij++;
			}

		if (a_norm != 0.0)
		{
			FT a_normEPS = a_norm * EPSILON;
			FT thr = a_norm;

			// Step 4: rotations
			int nb_iter = 0;
			while (thr > a_normEPS && nb_iter < MAX_ITER)
			{
				nb_iter++;
				FT thr_nn = thr / nn;

				for (int l = 1; l < n; l++)
				{
					for (int m = l + 1; m <= n; m++)
					{
						// compute sinx and cosx 
						int lq = (l*l - l) / 2;
						int mq = (m*m - m) / 2;

						int lm = l + mq;
						FT a_lm = a[lm];
						FT a_lm_2 = a_lm * a_lm;

						if (a_lm_2 < thr_nn)
							continue;

						int ll = l + lq;
						int mm = m + mq;
						FT a_ll = a[ll];
						FT a_mm = a[mm];

						FT delta = a_ll - a_mm;

						FT x;
						if (delta == 0.0)
							x = (FT)-M_PI / 4;
						else
							x = (FT)(-std::atan((a_lm + a_lm) / delta) / 2.0);

						FT sinx = std::sin(x);
						FT cosx = std::cos(x);
						FT sinx_2 = sinx * sinx;
						FT cosx_2 = cosx * cosx;
						FT sincos = sinx * cosx;

						// rotate L and M columns 
						int ilv = n * (l - 1);
						int imv = n * (m - 1);

						int i;
						for (i = 1; i <= n; i++)
						{
							if ((i != l) && (i != m))
							{
								int iq = (i*i - i) / 2;

								int im;
								if (i < m)
									im = i + mq;
								else
									im = m + iq;
								FT a_im = a[im];

								int il;
								if (i < l)
									il = i + lq;
								else
									il = l + iq;
								FT a_il = a[il];

								a[il] = a_il * cosx - a_im * sinx;
								a[im] = a_il * sinx + a_im * cosx;
							}

							ilv++;
							imv++;

							FT v_ilv = v[ilv];
							FT v_imv = v[imv];

							v[ilv] = cosx * v_ilv - sinx * v_imv;
							v[imv] = sinx * v_ilv + cosx * v_imv;
						}

						x = a_lm * sincos;
						x += x;

						a[ll] = a_ll * cosx_2 + a_mm * sinx_2 - x;
						a[mm] = a_ll * sinx_2 + a_mm * cosx_2 + x;
						a[lm] = 0.0;

						thr = std::fabs(thr - a_lm_2);
					}
				}
			}
		}

		// Step 5: index conversion and copy eigen values
		a++; // back from Fortran to C++

		for (i = 0; i < n; i++)
		{
			int k = i + (i*(i + 1)) / 2;
			eigen_values[i] = a[k];
		}
		delete[] a;

		// Step 6: sort the eigen values and eigen vectors 
		int *index = new int[n];
		for (i = 0; i < n; i++)
			index[i] = i;

		for (i = 0; i < (n - 1); i++)
		{
			FT x = eigen_values[i];
			int k = i;

			for (int j = i + 1; j < n; j++)
				if (x < eigen_values[j])
				{
					k = j;
					x = eigen_values[j];
				}

			eigen_values[k] = eigen_values[i];
			eigen_values[i] = x;

			int jj = index[k];
			index[k] = index[i];
			index[i] = jj;
		}


		// Step 7: save eigen vectors 
		v++; // back from Fortran to to C++

		ij = 0;
		for (int k = 0; k < n; k++)
		{
			int ik = index[k] * n;
			for (int i = 0; i < n; i++)
				eigen_vectors[ij++] = v[ik++];
		}

		delete[] v;
		delete[] index;
	}
}


#endif
