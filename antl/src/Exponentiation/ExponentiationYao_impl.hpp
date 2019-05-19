/**
 * @file ExponentiationYao_impl.hpp
 * @author Anton Mosunov
 * @brief implementation of Yao's algorithm (Yao, "On evaluation of powers", 1974)
 */

using namespace ANTL;

#include <NTL/vector.h>

template < class T >
void
ExponentiationYao<T>::power (T &C, const T &A, const ZZ &n)
{
	unsigned long b = NumBits(n);

	//cout << "b = " << b << endl;

	powersOf2.reserve(b - 1);

    long i = powersOf2.size();

	unsigned long k = NumBits(b);
	k = (k >> 1) + ((k & 1) == 1);
	unsigned long maxpow = k * ((b - 1) / k);

	//cout << "k = " << k << endl;

	if (i == 0)
	{
		//cout << i << ": " << A << endl;
		i++;
		powersOf2.push_back(A);
	}

	T B;
	assign(B,powersOf2[i - 1]);

	for (; i <= maxpow; i++)
	{
		sqr(B,B);
		powersOf2.push_back(B);
		//cout << i << ": " << B << endl;
	}

	// exponentiation base
	unsigned long D = (1L << k);

	vector<T> v;
	v.reserve(D - 1);

	for (i = 0; i < D - 1; i++)
	{
		v.push_back(conv<T>(NULL));
	}

	ZZ m = n;

	unsigned long digit;
	unsigned long mask = D - 1;

	i = 0;

	// compute d(1), d(2), ... , d(D-1)
	while (m != 0)
	{
		digit = conv<unsigned long>(m & mask);

		if (digit > 0)
		{
			if (v[digit - 1] == NULL)
			{
				assign(v[digit - 1], powersOf2[i]);
			}
			else
			{
				mul(v[digit - 1], v[digit - 1], powersOf2[i]);
			}
		}

		m >>= k;
		i += k;
	}

	//for (i = 0; i < D - 1; i++) cout << "d(" << (i + 1) << ") = " << v[i] << endl;

	// compute d(D-1), d(D-1) + d(D-2), ... , d(D-1) + ... + d(2) + d(1)
	for (i = D - 3; i >= 0; i--)
	{
		if (v[i+1] != NULL)
		{
			if (v[i] == NULL)
			{
				assign(v[i], v[i+1]);
			}
			else
			{
				mul(v[i], v[i], v[i+1]);
			}
		}
	}

	//for (i = 0; i < D - 1; i++) cout << "D(" << (i + 1) << ") = " << v[i] << endl;

	i = 0;

	while (v[i] == NULL) i++;

	assign(C, v[i++]);

	for (; i < D - 1; i++)
	{
		if (v[i] != NULL)
		{
			mul(C, C, v[i]);
		}
	}

	//cout << "RESULT: " << C << endl;

	return;
}
