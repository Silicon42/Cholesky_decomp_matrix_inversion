#include "cholesky_inv_sym_5.c"//"manual_unroll.c"//
#include <stdio.h>
/*
Correct output should be approximately this if using doubles
[0.00007770612, ],
[-0.54479757557, 0.00075790637, ],
[0.59134647493, -1.81843357380, 0.00108826553, ],
[-0.15633158700, 0.08094687885, -0.07096465742, 0.14228091784, ],
[-0.02055195534, 0.06828680765, -0.17002660150, -0.21531729686, 0.18793051518, ],

or this if using floats
[0.00007770611, ],
[-0.54479753971, 0.00075790618, ],
[0.59134626389, -1.81843328476, 0.00108826533, ],
[-0.15633150935, 0.08094671369, -0.07096459717, 0.14228069782, ],
[-0.02055184543, 0.06828679144, -0.17002657056, -0.21531811357, 0.18793073297, ],

*/
int main()
{
	float A[15] = {
		12869, 
		 7011, 5139, 
		 5139, 5199, 7334, 
		 1809, 1049,  903,  269, 
		 1049,  903, 1192,  177,  206
	};
	
	cholesky_inv_sym_5(A);
	
	for(int i = 0; i < 5; ++i)
	{
		putchar('[');
		for(int j = 0; j <= i; ++j)
		{
			printf("%.11f, ", A[TRI_INDEX(i,j)]);
		}
		puts("],");
	}

}