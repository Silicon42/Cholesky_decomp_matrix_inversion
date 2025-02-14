/* MIT License: Copyright (c) 2025 Curtis Stofer a.k.a. Silicon42
| 
| Permission is hereby granted, free of charge, to any person obtaining a copy
| of this software and associated documentation files (the "Software"), to deal
| in the Software without restriction, including without limitation the rights
| to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
| copies of the Software, and to permit persons to whom the Software is
| furnished to do so, subject to the following conditions:
| 
| The above copyright notice and this permission notice shall be included in all
| copies or substantial portions of the Software.
| 
| THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
| IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
| FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
| AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
| LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
| OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
| SOFTWARE.
\_____________________________________________________________________________*/

/*NOTE: all operations and intermediate results are done in-place on the packed
lower triangular 5x5 matrix. This means that the contents should be nearly 
guaranteed to stay in local memory since the cache lines will be frequently 
accessed. For floats the whole contents are 60 bytes ie just under the typical 
cache line size.
*/

/*NOTE: due to the way that the general conic is parameterized without a 6th 
constant component, there may be some instability for situations where the 
constant part would be very near 0, namely if the curve would pass extremely 
close to or through (0,0) so if you're using this for fitting conics like I am, 
be careful to avoid that
*/

#define TRI_INDEX(i,j)	(((int)(i)*((int)(i) + 1))/2 + (int)(j))

//easy type changing if you want another width of float
typedef float generic_float;

/*
Expects a representation of a symmetric positive definite matrix lower half 
packed in a linear array such that the index is as follows.
 0
 1  2
 3  4  5
 6  7  8  9
10 11 12 13 14

Does not check to verify the matrix satisfies these conditions.

Uses the L*D*L^T version of the Cholesky decomposition where D is a diagonal 
matrix and L is a lower unit triangular matrix.

Returns the coefficients of L**-1 in the lower triangle and the coefficiensts of
D**-1 on the diagonal
*/
void cholesky_inv_sym_5(generic_float A[15])
{
	// GCC only unrolls deepest loops but for this unrolling is both faster
	// and results in a shorter program
	#pragma GCC unroll 9
	// for each row
	for(int i = 1; i < 5; ++i)
	{	// calculate the D**-1 coeff, not stored as D coeff as a speed/accuracy 
		// tradeoff, doing this only uses 5 reciprocals + 15 multiplies and doing
		// the division normally uses 15 divides, if we assume a divide is 2x the
		// cycle cost of a multiply (or worse), the total cost is 5*2 + 15 = 25 
		// multiplies this way vs 15*2 = 30 multiplies when using as part of a 
		// solver for systems of linear equations, with 5 of the operations going
		// toward applying the diagonal to a column vector
		A[TRI_INDEX(i,-1)] = 1 / A[TRI_INDEX(i,-1)];
		// for each element in the current row before the diagonal
		for(int j = 0; j < i; ++j)
		{	// calculate the -L coefficient, negated b/c it's only ever needed as a negative
			// stored in temp variable so as to not overwrite the cell before using its S element
			generic_float L_temp = -A[TRI_INDEX(i,j)] * A[TRI_INDEX(j,j)];
			for(int k = i; k < 5; ++k)
				A[TRI_INDEX(k,i)] += L_temp * A[TRI_INDEX(k,j)];

			A[TRI_INDEX(i,j)] = L_temp;
		}
	}
	// get reciprocal of the last D coeff to convert to D**-1 coeff for consistency
	// could leave it and manually have a separate divison to slightly mitigate
	// the compounding accuracy loss, however any vectorizing wouldn't be able 
	// to handle that
	A[TRI_INDEX(4,4)] = 1 / A[TRI_INDEX(4,4)];

	// now that the L coefficients are fully calculated, they can be inverted by
	// row operations of subtracting multiples of lower rows. Since this would
	// zero out the element at that position, we can use the same element to 
	// store the resulting value that corresponds to the right side of the
	// augmented matrix

	#pragma GCC unroll 9
	// for each row, bottom first
	for(int i = 3; i > 0; --i)
		for(int j = i - 1; j >= 0; --j)	// for each element below the diagonal, right-most first
			for(int k = i + 1; k < 5; ++k)	// for each element below [i][j]
				A[TRI_INDEX(k,j)] += A[TRI_INDEX(i,j)] * A[TRI_INDEX(k,i)];

}