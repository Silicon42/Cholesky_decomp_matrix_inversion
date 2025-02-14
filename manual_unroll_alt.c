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

//NOTE: This is intended as a reference program for what the unrolled code should
// look like. Some compilers may not like this (looking at you MSVC) and may even
// do worse with this despite all loop indexing logic being stripped out. It's 
// mostly here as a control for benchmarking / correctness or just in case I run
// into a compiler that doesn't have good unrolling support.
#define TRI_INDEX(i,j)	(((int)(i)*((int)(i) + 1))/2 + (int)(j))

void cholesky_inv_sym_5(float A[15])
{
	//encouraging the compiler to prefetch
	float A00,A10,A11,A20,A21,A22,A30,A31,A32,A33,A40,A41,A42,A43,A44, L;
	A00 = A[0];
	A10 = A[1];
	A11 = A[2];
	A20 = A[3];
	A21 = A[4];
	A22 = A[5];
	A30 = A[6];
	A31 = A[7];
	A32 = A[8];
	A33 = A[9];
	A40 = A[10];
	A41 = A[11];
	A42 = A[12];
	A43 = A[13];
	A44 = A[14];

	//cholesky decomposition
	A00 = 1/A00;	// 1
	L = -A10*A00;	// 2*
	A11 += L*A10;	// 3
	A21 += L*A20;	// 3
	A31 += L*A30;	// 3
	A41 += L*A40;	// 3
	A10 = L;		// 2

	A11 = 1/A11;	// 4
	L = -A20*A00;	// 5*
	A22 += L*A20;	// 6
	A20 = L;		// 5
	L = -A21*A11;	// 7*
	A22 += L*A21;	// 8
	A21 = L;		// 7
	A32 += A20*A30 + A21*A31;	// 6, 8
	A42 += A20*A40 + A21*A41;	// 6, 8
	
	A22 = 1/A22;	// 9
	L = -A30*A00;	//10*
	A33 += L*A30;	//11
	A30 = L;		//10
	L = -A31*A11;	//12*
	A33 += L*A31;	//13
	A31 = L;		//12
	L = -A32*A22;	//14*
	A33 += L*A32;	//15
	A32 = L;		//14
	A43 = (A43 + A30*A40) + (A31*A41 + A32*A42);	//11, 13, 15

	A33 = 1/A33;	//16
	L = -A40*A00;	//17*
	A44 += L*A40;	//18
	A40 = L;		//17
	L = -A41*A11;	//19*
	A44 += L*A41;	//20
	A41 = L;		//19
	L = -A42*A22;	//21*
	A44 += L*A42;	//22
	A42 = L;		//21
	L = -A43*A33;	//23*
	A44 += L*A43;	//24
	A43 = L;		//23

	A44 = 1/A44;	//25

	//invert L matrix
	A42 += A32 * A43;
	A41 += A31 * A43;
	A40 += A30 * A43;
	A31 += A21 * A32;
	A41 += A21 * A42;
	A30 += A20 * A32;
	A40 += A20 * A42;
	A20 += A10 * A21;
	A30 += A10 * A31;
	A40 += A10 * A41;

	//store back to array
	A[0]  = A00;
	A[1]  = A10;
	A[2]  = A11;
	A[3]  = A20;
	A[4]  = A21;
	A[5]  = A22;
	A[6]  = A30;
	A[7]  = A31;
	A[8]  = A32;
	A[9]  = A33;
	A[10] = A40;
	A[11] = A41;
	A[12] = A42;
	A[13] = A43;
	A[14] = A44;
}