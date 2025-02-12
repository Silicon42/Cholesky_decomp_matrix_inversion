#include "manual_unroll.c"//"cholesky_inv_sym_5.c"//
#include <stdio.h>

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