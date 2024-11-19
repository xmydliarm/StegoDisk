#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>
#include <jpeglib.h>
#include <sys/stat.h>

using namespace std;

const int Q1 = 85;
const int Q2 = 70;

#pragma region ICDT/FDCT TRANSFORM

#define DCTSIZE 8
#define DCTSIZE2 64
#define CONST_BITS 13
#define PASS1_BITS 2
#define INT32 long
#define JSAMPLE unsigned char
#define ONE ((INT32)1)

#define FIX_0_298631336 ((INT32)2446)
#define FIX_0_390180644 ((INT32)3196)
#define FIX_0_541196100 ((INT32)4433)
#define FIX_0_765366865 ((INT32)6270)
#define FIX_0_899976223 ((INT32)7373)
#define FIX_1_175875602 ((INT32)9633)
#define FIX_1_501321110 ((INT32)12299)
#define FIX_1_847759065 ((INT32)15137)
#define FIX_1_961570560 ((INT32)16069)
#define FIX_2_053119869 ((INT32)16819)
#define FIX_2_562915447 ((INT32)20995)
#define FIX_3_072711026 ((INT32)25172)

#define MULTIPLY(var, const) ((var) * (const))
#define DEQUANTIZE(coef, quantval) (((int)(coef)) * (quantval))
#define RIGHT_SHIFT(x, shft) ((x) >> (shft))
#define DESCALE(x, n) RIGHT_SHIFT((x) + (ONE << ((n)-1)), n)

int range_limit(int x) {
    x = x + 128;
    if (x < 0) return 0;
    else if (x > 255) return 255;
    else return x;
}

void jpeg_idct_islow (double out[], const double coef_block[], const double qt[])
{
  INT32 tmp0, tmp1, tmp2, tmp3;
  INT32 tmp10, tmp11, tmp12, tmp13;
  INT32 z1, z2, z3, z4, z5;
  short * inptr, *intp;
  int * quantptr, *qtp;
  int * wsptr;
  short *outptr, *outbuf;
  int ctr, i;
  int workspace[DCTSIZE2];


  /* Pass 1: process columns from input, store into work array. */
  /* Note results are scaled up by sqrt(8) compared to a true IDCT; */
  /* furthermore, we scale the results by 2**PASS1_BITS. */

  /*inptr = coef_block;
  quantptr = (ISLOW_MULT_TYPE *) compptr->dct_table;*/
	inptr = (short*) malloc (DCTSIZE2 * sizeof(short));
	quantptr = (int*) malloc (DCTSIZE2 * sizeof(int));
	intp = inptr;
	qtp = quantptr;
	for (i=0; i<DCTSIZE2; i++) inptr[i] = (short) coef_block[i];
	for (i=0; i<DCTSIZE2; i++) quantptr[i] = (int) qt[i];

  wsptr = workspace;
  for (ctr = DCTSIZE; ctr > 0; ctr--) {
    /* Due to quantization, we will usually find that many of the input
     * coefficients are zero, especially the AC terms.  We can exploit this
     * by short-circuiting the IDCT calculation for any column in which all
     * the AC terms are zero.  In that case each output is equal to the
     * DC coefficient (with scale factor as needed).
     * With typical images and quantization tables, half or more of the
     * column DCT calculations can be simplified this way.
     */

    if (inptr[DCTSIZE*1] == 0 && inptr[DCTSIZE*2] == 0 &&
	inptr[DCTSIZE*3] == 0 && inptr[DCTSIZE*4] == 0 &&
	inptr[DCTSIZE*5] == 0 && inptr[DCTSIZE*6] == 0 &&
	inptr[DCTSIZE*7] == 0) {
      /* AC terms all zero */
      int dcval = DEQUANTIZE(inptr[DCTSIZE*0], quantptr[DCTSIZE*0]) << PASS1_BITS;

      wsptr[DCTSIZE*0] = dcval;
      wsptr[DCTSIZE*1] = dcval;
      wsptr[DCTSIZE*2] = dcval;
      wsptr[DCTSIZE*3] = dcval;
      wsptr[DCTSIZE*4] = dcval;
      wsptr[DCTSIZE*5] = dcval;
      wsptr[DCTSIZE*6] = dcval;
      wsptr[DCTSIZE*7] = dcval;

      inptr++;			/* advance pointers to next column */
      quantptr++;
      wsptr++;
      continue;
    }

    /* Even part: reverse the even part of the forward DCT. */
    /* The rotator is sqrt(2)*c(-6). */

    z2 = DEQUANTIZE(inptr[DCTSIZE*2], quantptr[DCTSIZE*2]);
    z3 = DEQUANTIZE(inptr[DCTSIZE*6], quantptr[DCTSIZE*6]);

    z1 = MULTIPLY(z2 + z3, FIX_0_541196100);
    tmp2 = z1 + MULTIPLY(z3, - FIX_1_847759065);
    tmp3 = z1 + MULTIPLY(z2, FIX_0_765366865);

    z2 = DEQUANTIZE(inptr[DCTSIZE*0], quantptr[DCTSIZE*0]);
    z3 = DEQUANTIZE(inptr[DCTSIZE*4], quantptr[DCTSIZE*4]);

    tmp0 = (z2 + z3) << CONST_BITS;
    tmp1 = (z2 - z3) << CONST_BITS;

    tmp10 = tmp0 + tmp3;
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;

    /* Odd part per figure 8; the matrix is unitary and hence its
     * transpose is its inverse.  i0  i3 are y7,y5,y3,y1 respectively.
     */

    tmp0 = DEQUANTIZE(inptr[DCTSIZE*7], quantptr[DCTSIZE*7]);
    tmp1 = DEQUANTIZE(inptr[DCTSIZE*5], quantptr[DCTSIZE*5]);
    tmp2 = DEQUANTIZE(inptr[DCTSIZE*3], quantptr[DCTSIZE*3]);
    tmp3 = DEQUANTIZE(inptr[DCTSIZE*1], quantptr[DCTSIZE*1]);

    z1 = tmp0 + tmp3;
    z2 = tmp1 + tmp2;
    z3 = tmp0 + tmp2;
    z4 = tmp1 + tmp3;
    z5 = MULTIPLY(z3 + z4, FIX_1_175875602); /* sqrt(2) * c3 */

    tmp0 = MULTIPLY(tmp0, FIX_0_298631336); /* sqrt(2) * (-c1+c3+c5-c7) */
    tmp1 = MULTIPLY(tmp1, FIX_2_053119869); /* sqrt(2) * ( c1+c3-c5+c7) */
    tmp2 = MULTIPLY(tmp2, FIX_3_072711026); /* sqrt(2) * ( c1+c3+c5-c7) */
    tmp3 = MULTIPLY(tmp3, FIX_1_501321110); /* sqrt(2) * ( c1+c3-c5-c7) */
    z1 = MULTIPLY(z1, - FIX_0_899976223); /* sqrt(2) * (c7-c3) */
    z2 = MULTIPLY(z2, - FIX_2_562915447); /* sqrt(2) * (-c1-c3) */
    z3 = MULTIPLY(z3, - FIX_1_961570560); /* sqrt(2) * (-c3-c5) */
    z4 = MULTIPLY(z4, - FIX_0_390180644); /* sqrt(2) * (c5-c3) */

    z3 += z5;
    z4 += z5;

    tmp0 += z1 + z3;
    tmp1 += z2 + z4;
    tmp2 += z2 + z3;
    tmp3 += z1 + z4;

    /* Final output stage: inputs are tmp10..tmp13, tmp0..tmp3 */

    wsptr[DCTSIZE*0] = (int) DESCALE(tmp10 + tmp3, CONST_BITS-PASS1_BITS);
    wsptr[DCTSIZE*7] = (int) DESCALE(tmp10 - tmp3, CONST_BITS-PASS1_BITS);
    wsptr[DCTSIZE*1] = (int) DESCALE(tmp11 + tmp2, CONST_BITS-PASS1_BITS);
    wsptr[DCTSIZE*6] = (int) DESCALE(tmp11 - tmp2, CONST_BITS-PASS1_BITS);
    wsptr[DCTSIZE*2] = (int) DESCALE(tmp12 + tmp1, CONST_BITS-PASS1_BITS);
    wsptr[DCTSIZE*5] = (int) DESCALE(tmp12 - tmp1, CONST_BITS-PASS1_BITS);
    wsptr[DCTSIZE*3] = (int) DESCALE(tmp13 + tmp0, CONST_BITS-PASS1_BITS);
    wsptr[DCTSIZE*4] = (int) DESCALE(tmp13 - tmp0, CONST_BITS-PASS1_BITS);

    inptr++;			/* advance pointers to next column */
    quantptr++;
    wsptr++;
  }

  /* Pass 2: process rows from work array, store into output array. */
  /* Note that we must descale the results by a factor of 8 == 2**3, */
  /* and also undo the PASS1_BITS scaling. */
	outptr = (short*) malloc (DCTSIZE * sizeof (short));
	outbuf = (short*) malloc (DCTSIZE2 * sizeof (short));
  wsptr = workspace;
  for (ctr = 0; ctr < DCTSIZE; ctr++) {
    /*outptr = output_buf[ctr] + output_col;*/
    /* Rows of zeroes can be exploited in the same way as we did with columns.
     * However, the column calculation has created many nonzero AC terms, so
     * the simplification applies less often (typically 5% to 10% of the time).
     * On machines with very fast multiplication, it's possible that the
     * test takes more time than it's worth.  In that case this section
     * may be commented out.
     */

#ifndef NO_ZERO_ROW_TEST
    if (wsptr[1] == 0 && wsptr[2] == 0 && wsptr[3] == 0 && wsptr[4] == 0 &&
	wsptr[5] == 0 && wsptr[6] == 0 && wsptr[7] == 0) {
      /* AC terms all zero */
      JSAMPLE dcval = range_limit((int) DESCALE((INT32) wsptr[0], PASS1_BITS+3)
				  );

      outptr[0] = dcval;
      outptr[1] = dcval;
      outptr[2] = dcval;
      outptr[3] = dcval;
      outptr[4] = dcval;
      outptr[5] = dcval;
      outptr[6] = dcval;
      outptr[7] = dcval;
			memcpy(&outbuf[ctr*DCTSIZE], outptr, DCTSIZE*sizeof(short));
      wsptr += DCTSIZE;		/* advance pointer to next row */
      continue;
    }
#endif

    /* Even part: reverse the even part of the forward DCT. */
    /* The rotator is sqrt(2)*c(-6). */

    z2 = (INT32) wsptr[2];
    z3 = (INT32) wsptr[6];

    z1 = MULTIPLY(z2 + z3, FIX_0_541196100);
    tmp2 = z1 + MULTIPLY(z3, - FIX_1_847759065);
    tmp3 = z1 + MULTIPLY(z2, FIX_0_765366865);

    tmp0 = ((INT32) wsptr[0] + (INT32) wsptr[4]) << CONST_BITS;
    tmp1 = ((INT32) wsptr[0] - (INT32) wsptr[4]) << CONST_BITS;

    tmp10 = tmp0 + tmp3;
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;

    /* Odd part per figure 8; the matrix is unitary and hence its
     * transpose is its inverse.  i0..i3 are y7,y5,y3,y1 respectively.
     */

    tmp0 = (INT32) wsptr[7];
    tmp1 = (INT32) wsptr[5];
    tmp2 = (INT32) wsptr[3];
    tmp3 = (INT32) wsptr[1];

    z1 = tmp0 + tmp3;
    z2 = tmp1 + tmp2;
    z3 = tmp0 + tmp2;
    z4 = tmp1 + tmp3;
    z5 = MULTIPLY(z3 + z4, FIX_1_175875602); /* sqrt(2) * c3 */

    tmp0 = MULTIPLY(tmp0, FIX_0_298631336); /* sqrt(2) * (-c1+c3+c5-c7) */
    tmp1 = MULTIPLY(tmp1, FIX_2_053119869); /* sqrt(2) * ( c1+c3-c5+c7) */
    tmp2 = MULTIPLY(tmp2, FIX_3_072711026); /* sqrt(2) * ( c1+c3+c5-c7) */
    tmp3 = MULTIPLY(tmp3, FIX_1_501321110); /* sqrt(2) * ( c1+c3-c5-c7) */
    z1 = MULTIPLY(z1, - FIX_0_899976223); /* sqrt(2) * (c7-c3) */
    z2 = MULTIPLY(z2, - FIX_2_562915447); /* sqrt(2) * (-c1-c3) */
    z3 = MULTIPLY(z3, - FIX_1_961570560); /* sqrt(2) * (-c3-c5) */
    z4 = MULTIPLY(z4, - FIX_0_390180644); /* sqrt(2) * (c5-c3) */

    z3 += z5;
    z4 += z5;

    tmp0 += z1 + z3;
    tmp1 += z2 + z4;
    tmp2 += z2 + z3;
    tmp3 += z1 + z4;

    /* Final output stage: inputs are tmp10 tmp13, tmp0 tmp3 */

    outptr[0] = range_limit((int) DESCALE(tmp10 + tmp3,
					  CONST_BITS+PASS1_BITS+3)
			    );
    outptr[7] = range_limit((int) DESCALE(tmp10 - tmp3,
					  CONST_BITS+PASS1_BITS+3)
			    );
    outptr[1] = range_limit((int) DESCALE(tmp11 + tmp2,
					  CONST_BITS+PASS1_BITS+3)
			    );
    outptr[6] = range_limit((int) DESCALE(tmp11 - tmp2,
					  CONST_BITS+PASS1_BITS+3)
			    );
    outptr[2] = range_limit((int) DESCALE(tmp12 + tmp1,
					  CONST_BITS+PASS1_BITS+3)
			    );
    outptr[5] = range_limit((int) DESCALE(tmp12 - tmp1,
					  CONST_BITS+PASS1_BITS+3)
			    );
    outptr[3] = range_limit((int) DESCALE(tmp13 + tmp0,
					  CONST_BITS+PASS1_BITS+3)
			    );
    outptr[4] = range_limit((int) DESCALE(tmp13 - tmp0,
					  CONST_BITS+PASS1_BITS+3)
			    );
    memcpy(&outbuf[ctr*DCTSIZE], outptr, DCTSIZE*sizeof(short));
    wsptr += DCTSIZE;		/* advance pointer to next row */
  }
	for (i=0; i<DCTSIZE2; i++) out[i] = (double) outbuf[i];
	free (intp);
	free (qtp);
	free (outptr);
	free (outbuf);
}

#define DCTELEM	int

void jpeg_fdct_islow (double out[], const double in[])
{
  INT32 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  INT32 tmp10, tmp11, tmp12, tmp13;
  INT32 z1, z2, z3, z4, z5;
  DCTELEM *dataptr, *data;
  int ctr, i;

  /* Pass 1: process rows. */
  /* Note results are scaled up by sqrt(8) compared to a true DCT; */
  /* furthermore, we scale the results by 2**PASS1_BITS. */
	data = (DCTELEM*) malloc (DCTSIZE * DCTSIZE * sizeof (DCTELEM));
	for (i=0; i<DCTSIZE * DCTSIZE; i++) data[i] = (int)in[i];

  dataptr = data;
  for (ctr = DCTSIZE-1; ctr >= 0; ctr--) {
    tmp0 = dataptr[0] + dataptr[7];
    tmp7 = dataptr[0] - dataptr[7];
    tmp1 = dataptr[1] + dataptr[6];
    tmp6 = dataptr[1] - dataptr[6];
    tmp2 = dataptr[2] + dataptr[5];
    tmp5 = dataptr[2] - dataptr[5];
    tmp3 = dataptr[3] + dataptr[4];
    tmp4 = dataptr[3] - dataptr[4];

    /* Even part per LL&M figure 1 --- note that published figure is faulty;
     * rotator "sqrt(2)*c1" should be "sqrt(2)*c6".
     */

    tmp10 = tmp0 + tmp3;
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;

    dataptr[0] = (DCTELEM) ((tmp10 + tmp11) << PASS1_BITS);
    dataptr[4] = (DCTELEM) ((tmp10 - tmp11) << PASS1_BITS);

    z1 = MULTIPLY(tmp12 + tmp13, FIX_0_541196100);
    dataptr[2] = (DCTELEM) DESCALE(z1 + MULTIPLY(tmp13, FIX_0_765366865),
				   CONST_BITS-PASS1_BITS);
    dataptr[6] = (DCTELEM) DESCALE(z1 + MULTIPLY(tmp12, - FIX_1_847759065),
				   CONST_BITS-PASS1_BITS);

    /* Odd part per figure 8 --- note paper omits factor of sqrt(2).
     * cK represents cos(K*pi/16).
     * i0..i3 in the paper are tmp4..tmp7 here.
     */

    z1 = tmp4 + tmp7;
    z2 = tmp5 + tmp6;
    z3 = tmp4 + tmp6;
    z4 = tmp5 + tmp7;
    z5 = MULTIPLY(z3 + z4, FIX_1_175875602); /* sqrt(2) * c3 */

    tmp4 = MULTIPLY(tmp4, FIX_0_298631336); /* sqrt(2) * (-c1+c3+c5-c7) */
    tmp5 = MULTIPLY(tmp5, FIX_2_053119869); /* sqrt(2) * ( c1+c3-c5+c7) */
    tmp6 = MULTIPLY(tmp6, FIX_3_072711026); /* sqrt(2) * ( c1+c3+c5-c7) */
    tmp7 = MULTIPLY(tmp7, FIX_1_501321110); /* sqrt(2) * ( c1+c3-c5-c7) */
    z1 = MULTIPLY(z1, - FIX_0_899976223); /* sqrt(2) * (c7-c3) */
    z2 = MULTIPLY(z2, - FIX_2_562915447); /* sqrt(2) * (-c1-c3) */
    z3 = MULTIPLY(z3, - FIX_1_961570560); /* sqrt(2) * (-c3-c5) */
    z4 = MULTIPLY(z4, - FIX_0_390180644); /* sqrt(2) * (c5-c3) */

    z3 += z5;
    z4 += z5;

    dataptr[7] = (DCTELEM) DESCALE(tmp4 + z1 + z3, CONST_BITS-PASS1_BITS);
    dataptr[5] = (DCTELEM) DESCALE(tmp5 + z2 + z4, CONST_BITS-PASS1_BITS);
    dataptr[3] = (DCTELEM) DESCALE(tmp6 + z2 + z3, CONST_BITS-PASS1_BITS);
    dataptr[1] = (DCTELEM) DESCALE(tmp7 + z1 + z4, CONST_BITS-PASS1_BITS);

    dataptr += DCTSIZE;		/* advance pointer to next row */
  }

  /* Pass 2: process columns.
   * We remove the PASS1_BITS scaling, but leave the results scaled up
   * by an overall factor of 8.
   */

  dataptr = data;
  for (ctr = DCTSIZE-1; ctr >= 0; ctr--) {
    tmp0 = dataptr[DCTSIZE*0] + dataptr[DCTSIZE*7];
    tmp7 = dataptr[DCTSIZE*0] - dataptr[DCTSIZE*7];
    tmp1 = dataptr[DCTSIZE*1] + dataptr[DCTSIZE*6];
    tmp6 = dataptr[DCTSIZE*1] - dataptr[DCTSIZE*6];
    tmp2 = dataptr[DCTSIZE*2] + dataptr[DCTSIZE*5];
    tmp5 = dataptr[DCTSIZE*2] - dataptr[DCTSIZE*5];
    tmp3 = dataptr[DCTSIZE*3] + dataptr[DCTSIZE*4];
    tmp4 = dataptr[DCTSIZE*3] - dataptr[DCTSIZE*4];

    /* Even part per LL&M figure 1 --- note that published figure is faulty;
     * rotator "sqrt(2)*c1" should be "sqrt(2)*c6".
     */

    tmp10 = tmp0 + tmp3;
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;

    dataptr[DCTSIZE*0] = (DCTELEM) DESCALE(tmp10 + tmp11, PASS1_BITS);
    dataptr[DCTSIZE*4] = (DCTELEM) DESCALE(tmp10 - tmp11, PASS1_BITS);

    z1 = MULTIPLY(tmp12 + tmp13, FIX_0_541196100);
    dataptr[DCTSIZE*2] = (DCTELEM) DESCALE(z1 + MULTIPLY(tmp13, FIX_0_765366865),
					   CONST_BITS+PASS1_BITS);
    dataptr[DCTSIZE*6] = (DCTELEM) DESCALE(z1 + MULTIPLY(tmp12, - FIX_1_847759065),
					   CONST_BITS+PASS1_BITS);

    /* Odd part per figure 8 --- note paper omits factor of sqrt(2).
     * cK represents cos(K*pi/16).
     * i0..i3 in the paper are tmp4..tmp7 here.
     */

    z1 = tmp4 + tmp7;
    z2 = tmp5 + tmp6;
    z3 = tmp4 + tmp6;
    z4 = tmp5 + tmp7;
    z5 = MULTIPLY(z3 + z4, FIX_1_175875602); /* sqrt(2) * c3 */

    tmp4 = MULTIPLY(tmp4, FIX_0_298631336); /* sqrt(2) * (-c1+c3+c5-c7) */
    tmp5 = MULTIPLY(tmp5, FIX_2_053119869); /* sqrt(2) * ( c1+c3-c5+c7) */
    tmp6 = MULTIPLY(tmp6, FIX_3_072711026); /* sqrt(2) * ( c1+c3+c5-c7) */
    tmp7 = MULTIPLY(tmp7, FIX_1_501321110); /* sqrt(2) * ( c1+c3-c5-c7) */
    z1 = MULTIPLY(z1, - FIX_0_899976223); /* sqrt(2) * (c7-c3) */
    z2 = MULTIPLY(z2, - FIX_2_562915447); /* sqrt(2) * (-c1-c3) */
    z3 = MULTIPLY(z3, - FIX_1_961570560); /* sqrt(2) * (-c3-c5) */
    z4 = MULTIPLY(z4, - FIX_0_390180644); /* sqrt(2) * (c5-c3) */

    z3 += z5;
    z4 += z5;

    dataptr[DCTSIZE*7] = (DCTELEM) DESCALE(tmp4 + z1 + z3,
					   CONST_BITS+PASS1_BITS);
    dataptr[DCTSIZE*5] = (DCTELEM) DESCALE(tmp5 + z2 + z4,
					   CONST_BITS+PASS1_BITS);
    dataptr[DCTSIZE*3] = (DCTELEM) DESCALE(tmp6 + z2 + z3,
					   CONST_BITS+PASS1_BITS);
    dataptr[DCTSIZE*1] = (DCTELEM) DESCALE(tmp7 + z1 + z4,
					   CONST_BITS+PASS1_BITS);

    dataptr++;			/* advance pointer to next column */
  }
	for (i=0; i<DCTSIZE*DCTSIZE; i++) out[i] = (double)data[i];
	free (data);
}

#pragma endregion ICDT/FDCT TRANSFORM

// Write 1D vector
template<typename T>
void writeVectorToFile(const std::vector<T>& vec1D, const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    for (size_t row = 0 ; row < vec1D.size() ; row++) {
        outfile << vec1D[row] << "\n";
    }

    outfile.close();
}

// Write 2D vector
template<typename T>
void writeVectorToFile(const std::vector<std::vector<T>>& vec2D, const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    for (size_t row = 0 ; row < vec2D.size() ; row++) {
        for (size_t column = 0 ; column < vec2D[0].size() ; column++) {
            outfile << vec2D[row][column] << "\n";
        }
    }

    outfile.close();
}

// Write 3D vector
template<typename T>
void writeVectorToFile(const std::vector<std::vector<std::vector<T>>>& vec3D, const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    for (size_t depth = 0 ; depth < vec3D[0][0].size() ; depth++) {
        for (size_t row = 0 ; row < vec3D.size() ; row++) {
            for (size_t column = 0 ; column < vec3D[0].size() ; column++) {
                outfile << vec3D[row][column][depth] << "\n";
            }
        }
    }

    outfile.close();
}

void print_quantization_table(jpeg_decompress_struct* cinfo) {
    for (int tblno = 0; tblno < NUM_QUANT_TBLS; tblno++) {
        JQUANT_TBL* qtbl = cinfo->quant_tbl_ptrs[tblno];
        if (qtbl != nullptr) {
            printf("Quantization table %d:\n", tblno + 1);
            for (int i = 0; i < DCTSIZE; i++) {
                for (int j = 0; j < DCTSIZE; j++) {
                    printf("%5d ", qtbl->quantval[i * DCTSIZE + j]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }
}

void Qmatrix(int quality, vector<vector<int>>& Q) {
    // Default 50% quality quantization matrix
    const vector<vector<int>> Q50 = {
        {16, 11, 10, 16, 24, 40, 51, 61},
        {12, 12, 14, 19, 26, 58, 60, 55},
        {14, 13, 16, 24, 40, 57, 69, 56},
        {14, 17, 22, 29, 51, 87, 80, 62},
        {18, 22, 37, 56, 68, 109, 103, 77},
        {24, 35, 55, 64, 81, 104, 113, 92},
        {49, 64, 78, 87, 103, 121, 120, 101},
        {72, 92, 95, 98, 112, 100, 103, 99}
    };

    // Ensure quality is between 1 and 100
    if (quality < 1) quality = 1;
    if (quality > 100) quality = 100;

    // Compute quantization matrix based on quality
    if (quality >= 50) {
        // High quality case (quality >= 50)
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                Q[i][j] = fmax(1, round(2 * Q50[i][j] * (1 - quality / 100.0)));
            }
        }
    } else {
        // Low quality case (quality < 50)
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                Q[i][j] = fmin(255, round(Q50[i][j] * 50.0 / quality));
            }
        }
    }
}

void extract_dct_coefficients(j_decompress_ptr cinfo, vector<vector<double>>& dct_coefficients) {
    // Read DCT coefficient arrays
    jvirt_barray_ptr *coef_arrays = jpeg_read_coefficients(cinfo);

    // Determine the total image dimensions (height, width) to initialize the 2D matrix
    const size_t total_height = cinfo->image_height;
    const size_t total_width = cinfo->image_width;

    // Resize the dct_coefficients to hold the full image (total_height x total_width)
    dct_coefficients.resize(total_height, vector(total_width, 0.0));

    // Loop over each component (Y, Cb, Cr)
    int current_height_offset = 0;  // Keep track of vertical offset for each component
    for (int ci = 0; ci < cinfo->num_components; ci++) {
        jpeg_component_info *compptr = &cinfo->comp_info[ci];

        // Determine the height and width of the DCT block grid for the current component
        auto block_rows = compptr->height_in_blocks;
        auto block_cols = compptr->width_in_blocks;
        auto c_height = block_rows * DCTSIZE;  // Full height in pixels

        // Loop through each block row (y direction)
        for (unsigned int blk_y = 0; blk_y < block_rows; blk_y++) {
            // Access the current row of blocks from the virtual array
            const JBLOCKARRAY buffer = (cinfo->mem->access_virt_barray)(reinterpret_cast<j_common_ptr>(cinfo), coef_arrays[ci], blk_y, 1, FALSE);

            // Loop through each block column (x direction)
            for (unsigned int blk_x = 0; blk_x < block_cols; blk_x++) {
                // Get the current DCT block (8x8)
                const JBLOCK &block = buffer[0][blk_x];

                // Copy DCT coefficients from the block to the correct position in the 2D matrix
                for (int i = 0; i < DCTSIZE; i++) {         // For each row in the block
                    for (int j = 0; j < DCTSIZE; j++) {     // For each column in the block
                        // Place the coefficient into the correct spot in the 2D matrix
                        const size_t row = blk_y * DCTSIZE + i + current_height_offset;
                        const size_t col = blk_x * DCTSIZE + j;
                        if (row < total_height && col < total_width) {
                            dct_coefficients[row][col] = static_cast<double>(block[i * DCTSIZE + j]);
                        }
                    }
                }
            }
        }

        // Increment the height offset for the next component
        current_height_offset += c_height;
    }
}

void print_dct_coefficients(vector<vector<double>>& dct_coefficients) {
    cout << "Full image DCT coefficients:\n";
    int row_count = 0;
    for (const auto& row : dct_coefficients) {
        cout << "Row " << row_count << ": ";
        row_count++;
        for (const auto& val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
}

void read_jpeg_info(FILE *infile, jpeg_decompress_struct& cinfo, vector<vector<double>>& dct_coefficients) {
    struct jpeg_error_mgr jerr{};

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    jpeg_read_header(&cinfo, TRUE);

    extract_dct_coefficients(&cinfo, dct_coefficients);
}

void printMatrix(const vector<vector<int>>& Q) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            cout << Q[i][j] << " "; // Dereference the shared_ptr and access elements
        }
        cout << endl;
    }
}

void DC_pairs(vector<vector<int>>& E) {
    const int MaxQ = 121; // Maximal quantization step

    // Calculate contributing pairs of quantization steps
    for (int i = 2; i < MaxQ; i++) {
        for (int j = i + 1; j < MaxQ; j++) {
            int g = gcd(i, j);
            if ((j / g) % 2 == 0) {
                E[i - 1][j - 1] = j / g;
            }
        }
    }
}

void PlaneToVec(const vector<vector<double>>& plane, const int Y, const int X, vector<vector<vector<int>>>& D1) {
    const int M = X / 8;  // Count of columns
    const int N = Y / 8;  // Count of rows

    // Resize the D1 vector to hold the reshaped blocks
    D1.resize(N, vector(M, vector(64, 0)));

    // Fill the D1 array with reshaped blocks from the plane in column-major order
    for (int i = 0; i < N; ++i) { // Rows of blocks
        for (int j = 0; j < M; ++j) { // Columns of blocks
            int idx = 0;
            // Iterate in column-major order for an 8x8 block
            for (int n = 0; n < 8; ++n) { // Columns inside block (MATLAB-style: first move column)
                for (int m = 0; m < 8; ++m) { // Rows inside block
                    D1[i][j][idx++] = static_cast<int>(plane[i * 8 + m][j * 8 + n]);
                }
            }
        }
    }
}

void VecToPlane(vector<vector<double>>& plane, const int MB, const int NB, const vector<vector<vector<double>>>& cube) {
    // Define the dimensions for rows and columns based on Cube size
    const int M = NB * 8; // Count of columns in Plane
    const int N = MB * 8; // Count of rows in Plane

    // Resize the Plane vector to match the required dimensions
    plane.resize(N, vector<double>(M, 0.0));

    // Fill the Plane array with data from Cube in 8x8 blocks
    for (int i = 0; i < MB; ++i) { // Rows of 8x8 blocks
        for (int j = 0; j < NB; ++j) { // Columns of 8x8 blocks
            int idx = 0;
            // Copy values from Cube(i,j,:) into an 8x8 block in Plane
            for (int n = 0; n < 8; ++n) { // Rows within 8x8 block
                for (int m = 0; m < 8; ++m) { // Columns within 8x8 block
                    plane[i * 8 + m][j * 8 + n] = cube[i][j][idx++];
                }
            }
        }
    }
}

vector<vector<double>> showD(const vector<vector<vector<int>>>& D, const vector<vector<int>>& Q) {
    int B = 8;               // Block size
    int MD = D.size();       // Number of 8x8 blocks in the vertical direction
    int ND = D[0].size();    // Number of 8x8 blocks in the horizontal direction
    int M = B * MD;          // Full image height in pixels
    int N = B * ND;          // Full image width in pixels

    // Decompressed image matrix
    vector X(M, vector<double>(N, 0.0));

    // Temporary block to store the spatial representation of a DCT block
    vector Dblock(8, vector<int>(8));

    // Loop over each block
    for (int i = 0; i < MD; i++) {
        for (int j = 0; j < ND; j++) {
            // Copy the DCT coefficients of the (i, j)-th block into Dblock
            for (int k = 0; k < B; k++) {
                for (int l = 0; l < B; l++) {
                    Dblock[k][l] = D[i][j][k * B + l]; // Convert flat 64 to 8x8 matrix
                }
            }

            // Calculate the starting pixel positions for the block in the decompressed image
            int ib = i * B;
            int jb = j * B;

            vector<double> out = vector<double>(64, 0);
            vector<double> quant = vector<double>(64, 0);
            vector<double> coef_block_linear = vector<double>(64, 0);

            int counter = 0;
            for (size_t pp = 0 ; pp < Dblock.size(); pp++){
                for (size_t pq = 0 ; pq < Dblock[pp].size(); pq++)
                {
                    coef_block_linear[counter] = Dblock[pq][pp];
                    counter++;
                }
            }

            counter = 0;
            for (size_t k = 0 ; k < Q.size(); k++){
                for (size_t l = 0 ; l < Q[k].size(); l++)
                {
                    quant[counter] = Q[k][l];
                    counter++;
                }
            }

            // All are row-first matrices
            double out_array[64];
            double quant_array[64];
            double coefs_array[64];

            for (int k = 0; k < 64; ++k) {
                out_array[k] = out[k]; // Example values
                quant_array[k] = quant[k]; // Example values
                coefs_array[k] = coef_block_linear[k]; // Example values
            }

            //FUNCTION HEADER: jpeg_idct_islow (double out[], double coef_block[], double qt[])
            jpeg_idct_islow(out_array, coefs_array, quant_array);

            // Place the spatial block into the decompressed image X
            counter = 0;
            for (int k = 0; k < B; k++)
            {
                for (int l = 0; l < B; l++) {
                    X[ib + k][jb + l] = out_array[counter];
                    counter++;
                }
            }
        }
    }

    return X;  // Return the decompressed image
}

void save_jpeg(const char* filename, const vector<vector<double>>& image, int quality) {
    // Define image dimensions based on X1
    int width = image[0].size();
    int height = image.size();

    // Initialize JPEG compression object
    jpeg_compress_struct cinfo{};
    jpeg_error_mgr jerr{};

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    // Open file for output
    FILE* outfile;
    if ((outfile = fopen(filename, "wb")) == nullptr) {
        cerr << "Error: Can't open output file " << filename << endl;
        return;
    }

    jpeg_stdio_dest(&cinfo, outfile);

    // Set image parameters
    cinfo.image_width = width;
    cinfo.image_height = height;
    cinfo.input_components = 1; // 1 for grayscale or 3 for RGB
    cinfo.in_color_space = JCS_GRAYSCALE; // or JCS_RGB for color images

    // Set default compression parameters and set quality
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);

    // Start compression
    jpeg_start_compress(&cinfo, TRUE);

    // Allocate memory for row buffer
    JSAMPROW row_pointer[1];
    int row_stride = width; // Number of bytes in a row (width of image)

    // Write each row of the image to the output JPEG
    while (cinfo.next_scanline < cinfo.image_height) {
        // Convert the vector data to unsigned char (expected by libjpeg)
        vector<unsigned char> row(width);
        for (int i = 0; i < width; ++i) {
            row[i] = static_cast<unsigned char>(image[cinfo.next_scanline][i]);
        }
        row_pointer[0] = &row[0];  // Row to write
        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    // Finish compression
    jpeg_finish_compress(&cinfo);

    // Clean up and close the file
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);

    cout << "JPEG image saved successfully: " << filename << endl;
}

// Function to count non-zeros in a 3D vector for specified range
int count_nonzeros(const vector<vector<vector<int>>>& matrix, const int start, const int end) {
    int non_zero_count = 0;
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            for (int k = start; k <= end; ++k) {
                if (matrix[i][j][k] != 0) {
                    non_zero_count++;
                }
            }
        }
    }
    return non_zero_count;
}

// Function to calculate capacity and number of non-zero coefficients
void DC_capacity_indiv(const vector<vector<int>>& Qm1, const vector<vector<int>>& Qm2,
                       const vector<vector<vector<int>>>& D, const vector<vector<int>>& E,
                       size_t& Cap, int& N0) {
    Cap = 0;

    // Iterate over the 64 DCT modes (8x8 blocks)
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            int q1 = Qm1[j][i];  // First quantization step
            int q2 = Qm2[j][i];  // Second quantization step

            // If q1, q2 is a contributing pair
            if (E[q1-1][q2-1] != 0 && E[q1-1][q2-1] % 2 == 0) {
                // Calculate capacity by determining whether the koeficient is multiple of quantization step
                size_t start = E[q1-1][q2-1] * 0.5;
                size_t step = E[q1-1][q2-1];   // Step size

                int mode_index = i * 8 + j;  // Calculate mode index within 8x8 matrix
                // Iterate over all blocks in D for the current mode
                for (size_t x = 0; x < D.size(); ++x) {
                   for (size_t y = 0; y < D[0].size(); ++y) {
                       int coef = D[x][y][mode_index];  // Get the coefficient for current block and mode

                       // Check if the coefficient is a multiple of step size
                       if (coef != 0) {
                           N0++;  // Count non-zero coefficients
                           if (abs(coef - static_cast<int>(start)) % step == 0) {
                               Cap++;  // Increment capacity if coefficient is contributing
                           }
                       }
                   }
                }
            }
        }
    }

    // Count the total number of non-zero coefficients in D
    N0 = 0;
    for (const auto & i : D) {
        for (const auto & j : i) {
            for (const int k : j) {
                if (k != 0) {
                    N0++;
                }
            }
        }
    }

    cout << "Capacity: " << Cap << endl;
    cout << "Nonzero coefs of image: " << N0 << endl;
}

// raw2jpg02raw function for quantized DCT coefficients
vector<vector<double>> raw2jpg02raw(const vector<vector<double>>& Z, const vector<vector<double>>& Qf) {
    std::vector Z_shifted(8, std::vector(8, 0.0));

    // Shift matrix values by subtracting 128 from each element
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            Z_shifted[i][j] = Z[i][j] - 128;
        }
    }

    double out[64];
    double in[64];

    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
            in[8 * i + j] = Z_shifted[i][j];

    jpeg_fdct_islow(out, in);

    // Column-first
    std::vector D(8, vector(8, 0.0));
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
            D[i][j] = out[ 8 * i + j];

    // Quantize the DCT coefficients
    std::vector QD(8, std::vector(8, 0.0));
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            QD[i][j] = D[i][j] / (Qf[i][j] * 8);
        }
    }

    return QD;
}

void DCTcutraw(const vector<vector<double>>& X, vector<vector<vector<double>>>& D_raw) {
    const int B = 8;  // Block size
    const int B2 = B * B;

    const int M = X.size();  // Number of rows
    const int N = X[0].size();  // Number of columns

    const int MB = M / B;  // Number of blocks along rows
    const int NB = N / B;  // Number of blocks along columns

    // Resize the D matrix to hold the DCT coefficients
    D_raw.resize(MB, std::vector(NB, std::vector<double>(B2, 0)));

    // Loop over all blocks
    for (int i = 0; i < MB; ++i) {
        for (int j = 0; j < NB; ++j) {
            int ib = i * B;
            int jb = j * B;

            // Extract the 8x8 block from the original matrix X
            std::vector Block(B, std::vector<double>(B, 0));
            for (int m = 0; m < B; ++m) {
                for (int n = 0; n < B; ++n) {
                    Block[m][n] = static_cast<double>(X[ib + m][jb + n]);
                }
            }

            // Perform DCT (raw2jpg02raw function to be implemented later)
            std::vector<std::vector<double>> Dblock = raw2jpg02raw(Block, std::vector(B, std::vector<double>(B, 1)));

            // Flatten the DCT coefficients and assign them to the D_raw matrix (correcting the indexing)
            int k = 0;
            for (int m = 0; m < B; ++m) {
                for (int n = 0; n < B; ++n) {
                    D_raw[i][j][k] = Dblock[n][m];
                    ++k;
                }
            }
        }
    }
}

void generateMessage(vector<int>& Message, const size_t Cap, const size_t AbsoluteMsgLength) {
    // Set the message length to the minimum of Cap and AbsoluteMsgLength
    size_t MessageLength = std::min(Cap, AbsoluteMsgLength);

    // Bias for generating random 1's and -1's
    const double Bias = 0.5;

    // Initialize the random number generator with a fixed seed
    //std::mt19937 engine(12345); // Fixed seed for reproducible results
    //std::uniform_real_distribution<double> dist(0.0, 1.0); // Uniform distribution between [0, 1]

    // Resize the Message vector to hold the message
    Message.resize(MessageLength);

    // Generate random message of 1's and -1's
    for (size_t i = 0; i < MessageLength; ++i) {
        // Generate a random value in the range [0,1]
        //double randomValue = dist(engine);
        //Message[i] = (randomValue > Bias) ? 1 : -1;
        Message[i] = -1;
    }
}

void sortRowsEmulation(vector<double>& diff_vector, vector<int>& indices) {
    const size_t n = diff_vector.size();

    for (size_t i = 0; i < n; ++i) {
        indices[i] = i;
    }

    auto comparator = [](const std::pair<double, size_t>& a, const std::pair<double, size_t>& b) {
        return a.first < b.first;
    };

    vector<std::pair<double, size_t>> paired_vector(n);
    for (size_t i = 0; i < n; ++i) {
        paired_vector[i] = {diff_vector[i], indices[i]};
    }

    std::stable_sort(paired_vector.begin(), paired_vector.end(), comparator);

    for (size_t i = 0; i < n; ++i) {
        diff_vector[i] = paired_vector[i].first; // Update sorted values
        indices[i] = paired_vector[i].second;    // Update indices to original order
    }
}

int roundToInt(const double val) {
    return static_cast<int>(std::round(val));
}

int sign(const int val) {
    return (val > 0) - (val < 0);
}

vector<vector<vector<double>>> embedMessage(const vector<vector<int>>& Qm1, const vector<vector<int>>& Qm2,
                  const vector<vector<int>>& E, const vector<vector<vector<int>>>& D1,
                  const vector<double>& Diff_vector, const vector<int>& Message,
                  const vector<vector<vector<double>>>& D2raw,
                  const string& nonzero_spec, const int MB, const int NB) {
    int c = -1;                   // Message bit counter
    int d = -1;                   // Contributing multiples counter
    int changes = 0;             // Number of changes against cover
    int NonZeros = 0, Zeros = 0; // Zero and non-zero counters

    std::vector D2n(MB, std::vector(NB, std::vector<double>(64)));

    for (int i = 0; i < MB; ++i) {
        for (int j = 0; j < NB; ++j) {
            for (int k = 0; k < 64; ++k) {
                D2n[i][j][k] = static_cast<double>(D1[i][j][k]);
            }
        }
    }

    for (int k = 0; k < 64; ++k) {
        int ii = k / 8; // row in 8x8 block
        int jj = k % 8; // column in 8x8 block

        int q1 = Qm1[jj][ii];  // First quantization step
        int q2 = Qm2[jj][ii];  // Second quantization step

        if (E[q1-1][q2-1] != 0 && E[q1-1][q2-1] % 2 == 0) {
            for (int i = 0; i < MB; ++i) {
                for (int j = 0; j < NB; ++j) {
                    if (static_cast<int>((D1[i][j][k] - E[q1 - 1][q2 - 1] * 0.5)) % E[q1 - 1][q2 - 1] == 0) {
                        // D1[i][j][k] is the contributing multiple
                        d++;
                        if (Diff_vector[d] == 1) {
                            // Embed one message bit
                            c++;
                            D2n[i][j][k] = (D1[i][j][k] * q1 + sign(D1[i][j][k]) * Message[c] * q2 / 2) / q2;

                            int cover_value = roundToInt(D2raw[i][j][k] / q2);
                            if (std::abs(q2 * cover_value - q1 * roundToInt(q2 * cover_value / q1)) == q2 / 2) {
                                cover_value += sign(D2raw[i][j][k] / q2 - cover_value);
                            }

                            if (D2n[i][j][k] != cover_value) {
                                changes++;
                            }
                        } else {
                            // D1[i][j][k] is not suitable for embedding
                            D2n[i][j][k] = roundToInt(D2raw[i][j][k] / q2);
                            if (std::abs(q2 * D2n[i][j][k] - q1 * roundToInt(q2 * D2n[i][j][k] / q1)) == q2 / 2) {
                                D2n[i][j][k] += sign(D2raw[i][j][k] / q2 - D2n[i][j][k]);
                            }
                        }
                    } else {
                        // D1[i][j][k] is not a contributing multiple
                        D2n[i][j][k] = roundToInt(D2raw[i][j][k] / q2);

                        if (std::abs(q2 * D2n[i][j][k] - q1 * roundToInt(q2 * D2n[i][j][k] / q1)) == q2 / 2) {
                            D2n[i][j][k] += sign(D2raw[i][j][k] / q2 - D2n[i][j][k]);

                        }
                    }

                    // Update Zero/NonZero counters
                    if (nonzero_spec == "DC-DC" || nonzero_spec == "DC-AC") {
                        if (D2n[i][j][k]) NonZeros++;
                        else Zeros++;
                    } else {
                        if (D2n[i][j][k] && k > 0) NonZeros++;
                        else Zeros++;
                    }
                }
            }
        } else {
            // Non-contributing pair, round to closest multiple of q2
            for (int i = 0; i < MB; ++i) {
                for (int j = 0; j < NB; ++j) {
                    D2n[i][j][k] = roundToInt(D2raw[i][j][k] / q2);
                    if (nonzero_spec == "DC-DC" || nonzero_spec == "DC-AC") {
                        if (D2n[i][j][k]) NonZeros++;
                        else Zeros++;
                    } else {
                        if (D2n[i][j][k] && k > 0) NonZeros++;
                        else Zeros++;
                    }
                }
            }
        }
    }

    cout << "\n ---- Message embedded! ----" << endl;
    std::cout << "Changes: " << changes << "\n";
    std::cout << "Zeros: " << Zeros << ", NonZeros: " << NonZeros << "\n";

    return D2n;
}

int save_perturbed_jpeg(vector<vector<double>>& Plane) {
    FILE *infile = fopen("/home/mato/CLionProjects/StegoDisk/cover.jpg", "rb");
    if (!infile) {
        fprintf(stderr, "Can't open %s\n", "/home/mato/CLionProjects/StegoDisk/cover.jpg");
        return -1;
    }

    jpeg_decompress_struct cinfo;
    jpeg_error_mgr jerr;
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    jpeg_read_header(&cinfo, TRUE);

    jvirt_barray_ptr* coef_arrays = jpeg_read_coefficients(&cinfo);

    jpeg_compress_struct cinfo_out;
    jpeg_error_mgr jerr_out;

    cinfo_out.err = jpeg_std_error(&jerr_out);
    jpeg_create_compress(&cinfo_out);
    jpeg_copy_critical_parameters(&cinfo, &cinfo_out);

    FILE *outfile = fopen("/home/mato/CLionProjects/StegoDisk/output.jpg", "wb");
    if (outfile == nullptr) {
        fprintf(stderr, "Can't open %s\n", "/home/mato/CLionProjects/StegoDisk/output.jpg");
        jpeg_destroy_compress(&cinfo_out);
        return -1;
    }
    jpeg_stdio_dest(&cinfo_out, outfile);

    for (int ci = 0; ci < cinfo.num_components; ci++) {
        jpeg_component_info *compptr = &cinfo.comp_info[ci];
        auto block_rows = compptr->height_in_blocks;
        auto block_cols = compptr->width_in_blocks;

        assert(block_rows * DCTSIZE <= Plane.size());
        assert(block_cols * DCTSIZE <= Plane[0].size());

        for (unsigned int blk_y = 0; blk_y < block_rows; blk_y++) {
            JBLOCKARRAY buffer = (cinfo.mem->access_virt_barray)
                (reinterpret_cast<j_common_ptr>(&cinfo), coef_arrays[ci], blk_y, 1, TRUE);

            for (unsigned int blk_x = 0; blk_x < block_cols; blk_x++) {
                for (int i = 0; i < 8; i++) {
                    for (int j = 0; j < 8; j++) {
                        int plane_row = blk_y * 8 + i;
                        int plane_col = blk_x * 8 + j;
                        buffer[0][blk_x][i * DCTSIZE + j] = static_cast<JCOEF>(Plane[plane_row][plane_col]);
                    }
                }
            }
        }
    }

    jpeg_write_coefficients(&cinfo_out, coef_arrays);

    jpeg_finish_compress(&cinfo_out);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo_out);

    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);

    return 0;
}

int main() {
    cout << "\n ---- Algorithm started successfully! ----" << endl;

    std::srand(12345);

    auto filename = "/home/mato/CLionProjects/StegoDisk/img_1200.jpg";

    FILE *infile;
    if ((infile = fopen(filename, "rb")) == nullptr) {
        fprintf(stderr, "Can't open %s\n", filename);
        return -1;
    }

    jpeg_decompress_struct cinfo{};
    auto dct_coefficients = vector<vector<double>>();

    read_jpeg_info(infile, cinfo, dct_coefficients);

    int Y = cinfo.image_height;
    int X = cinfo.image_width;

    auto Qm1 = vector(8, vector<int>(8));
    auto Qm2 = vector(8, vector<int>(8));

    // Generate quantization matrix for the given quality
    Qmatrix(Q1, Qm1);
    Qmatrix(Q2, Qm2);

    // Print the generated quantization matrix
    /*cout << "Quantization Matrix Qm1 for quality " << Q1 << endl;
    printMatrix(Qm1);
    cout << "Quantization Matrix Qm2 for quality " << Q2 << endl;
    printMatrix(Qm2);*/

    // Create a 3D vector with dimensions [121][121][3]
    auto E = vector(121, vector<int>(121));

    // Call the DC_pairs function
    DC_pairs(E);

    // Create the output vector D1
    auto D1 = vector<vector<vector<int>>>();

    // Call the PlaneToVec function
    PlaneToVec(dct_coefficients, Y, X, D1);

    // Decompress D1 to the spatial domain
    auto X1 = showD(D1, Qm1);

    int MB = floor(X1.size() / 8);
    int NB = floor(X1[0].size() / 8);

    //second compression, write the compressed image into cover directory
    save_jpeg("/home/mato/CLionProjects/StegoDisk/cover.jpg", X1, Q2);

    // dispose resources
    fclose(infile);
    jpeg_destroy_decompress(&cinfo);


    //-------------------------------
    // SECOND COMPRESSION STARTS HERE
    //-------------------------------


    FILE *infileTwo;
    if ((infileTwo = fopen("/home/mato/CLionProjects/StegoDisk/cover.jpg", "rb")) == nullptr) {
        fprintf(stderr, "Can't open %s\n", "/home/mato/CLionProjects/StegoDisk/cover.jpg");
        return -1;
    }

    jpeg_decompress_struct cinfoTwo{};
    auto dct_coefficients_two = vector<vector<double>>();

    read_jpeg_info(infileTwo, cinfoTwo, dct_coefficients_two);

    // Create the output vector D1
    auto D1_two = vector<vector<vector<int>>>();

    // Call the PlaneToVec function
    PlaneToVec(dct_coefficients_two, Y, X, D1_two);

    auto DCTs_DoubleCompressed = D1_two;

    // For storing absolute message length
    const double capacity = 0.1;
    size_t AbsoluteMsgLength = 0;
    const string nonzero_spec = "AC-DC";

    // 'DC-DC' condition
    if (nonzero_spec == "DC-DC") {
        // Include all terms (DC included) in capacity calculations
        int Nz_all = count_nonzeros(DCTs_DoubleCompressed, 0, 63); // Count non-zeros in entire matrix
        AbsoluteMsgLength = round(Nz_all * capacity);
    }

    // 'AC-DC' condition
    if (nonzero_spec == "AC-DC") {
        // Do not include DC terms in capacity calculations, but use them for embedding
        int Nz_ac = count_nonzeros(DCTs_DoubleCompressed, 1, 63); // Count non-zeros in AC terms (ignoring the first DC component)
        AbsoluteMsgLength = round(Nz_ac * capacity);
    }

    // 'AC-AC' condition
    if (nonzero_spec == "AC-AC") {
        // Do not include DC terms in capacity calculations and do not embed into them
        int Nz_ac = count_nonzeros(DCTs_DoubleCompressed, 1, 63); // Count non-zeros in AC terms
        AbsoluteMsgLength = round(Nz_ac * capacity);

        // Modify the quantization matrix to make DC terms non-contributing
        Qm1[0][0] = 6;
    }

    size_t Cap;
    int N0;

    // Calculate the capacity and the number of non-zero coefficients
    // embedding capacity calculation (after the first compression)
    DC_capacity_indiv(Qm1, Qm2, D1, E, Cap, N0);

    auto D2raw = vector<vector<vector<double>>>();

    // Second compression, without quantization
    DCTcutraw(X1, D2raw);

    // Generate random secret message 1's and -1's of length MessageLength
    vector<int> Message;
    size_t MessageLength = Cap > AbsoluteMsgLength ? AbsoluteMsgLength : Cap;
    generateMessage(Message, Cap, MessageLength);

    // Initialize Diff_vector and Diff_matrix
    const string method = "pqt";
    vector Diff_vector(Cap, 0.0);
    vector Diff_matrix(MB, vector(NB, vector(64, 0.0)));

    if (method == "midpoint" || method == "pq" || method == "PQ") {
        size_t c = 0;  // Message bit counter

        // Iterate over the 64 DCT modes
        for (int k = 0; k < 64; ++k) {
            int q1 = Qm1[k % 8][k / 8];  // q1 = first quantization step
            int q2 = Qm2[k % 8][k / 8];  // q2 = second quantization step

            //if (E[q1-1][q2-1][0] != 0) {  // If q1, q2 is a contributing pair
            if (E[q1-1][q2-1] != 0 && E[q1-1][q2-1] % 2 == 0) {
                // Loop over all blocks
                for (int i = 0; i < MB; ++i) {
                    for (int j = 0; j < NB; ++j) {
                        //if (std::fmod(D1[i][j][k] - E[q1-1][q2-1][1], E[q1-1][q2-1][2]) == 0) {
                        if (std::fmod(D1[i][j][k] - E[q1-1][q2-1] * 0.5, E[q1-1][q2-1]) == 0) {
                            // D1(i,j,k) is the contributing multiple
                            c++;  // Update message bit counter

                            // Count the diff_vector value
                            double D2_scaled = D2raw[i][j][k] / static_cast<double>(q2);
                            Diff_vector[c - 1] = std::fabs(std::fabs(D2_scaled - std::round(D2_scaled)) - 0.5);
                        }
                    }
                }
            }
        }
    }

    if (method == "texture" || method == "pqt" || method == "PQt") {
        int c = 0;  // Message bit counter

        vector X1_T(MB, std::vector(NB, std::vector(64, 0)));

        PlaneToVec(X1, Y, X, X1_T);

        for (int i = 0; i < MB; ++i) {
            for (int j = 0; j < NB; ++j) {
                // Calculate the texture function 'Energy'
                double Energy = 0.0;

                // Calculate energy based on 8x8 block of X1_T
                for (int k = 0; k < 8; k += 2) {
                    for (int l = 0; l < 8; l += 2) {
                        // Extract values from X1_T (assuming 64 coefficients per block)
                        std::vector<int> values = {
                            X1_T[i][j][k * 8 + l],
                            X1_T[i][j][k * 8 + l + 1],
                            X1_T[i][j][(k + 1) * 8 + l],
                            X1_T[i][j][(k + 1) * 8 + l + 1]
                        };

                        // Calculate the max and min difference and update Energy
                        Energy += *std::max_element(values.begin(), values.end()) -
                                  *std::min_element(values.begin(), values.end());
                    }
                }

                // Iterate over the 64 DCT modes
                for (int k = 0; k < 64; ++k) {
                    int ii = k / 8; // row in 8x8 block
                    int jj = k % 8; // column in 8x8 block

                    int q1 = Qm1[jj][ii];  // First quantization step
                    int q2 = Qm2[jj][ii];  // Second quantization step

                    // Check if q1, q2 is a contributing pair
                    if (E[q1-1][q2-1] != 0 && E[q1-1][q2-1] % 2 == 0) {
                        // Check for contributing multiple
                        if (static_cast<int>(D1[i][j][k] - E[q1 - 1][q2 - 1] * 0.5) % E[q1 - 1][q2 - 1] == 0) {
                            // Increment message bit counter
                            ++c;

                            // Update Diff_matrix
                            Diff_matrix[i][j][k] = -Energy - 1;
                        }
                    }
                }
            }
        }

        // Create the diff_vector
        int counter = 0;
        for (int k = 0; k < 64; ++k) {
            for (int i = 0; i < MB; ++i) {
                for (int j = 0; j < NB; ++j) {
                    if (Diff_matrix[i][j][k] < 0) {
                        Diff_vector[counter] = Diff_matrix[i][j][k];
                        ++counter;
                    }
                }
            }
        }
    }

    if (method == "-pqt") {
        int c = 0;  // Message bit counter

        vector X1_T(MB, std::vector(NB, std::vector(64, 0)));

        PlaneToVec(X1, Y, X, X1_T);

        for (int i = 0; i < MB; ++i) {
            for (int j = 0; j < NB; ++j) {
                // Calculate the texture function 'Energy'
                double Energy = 0.0;

                // Calculate energy based on 8x8 block of X1_T
                for (int k = 0; k < 8; k += 2) {
                    for (int l = 0; l < 8; l += 2) {
                        // Extract values from X1_T (assuming 64 coefficients per block)
                        std::vector<int> values = {
                            X1_T[i][j][k * 8 + l],
                            X1_T[i][j][k * 8 + l + 1],
                            X1_T[i][j][(k + 1) * 8 + l],
                            X1_T[i][j][(k + 1) * 8 + l + 1]
                        };

                        // Calculate the max and min difference and update Energy
                        Energy += *std::max_element(values.begin(), values.end()) -
                                  *std::min_element(values.begin(), values.end());
                    }
                }

                // Iterate over the 64 DCT modes
                for (int k = 0; k < 64; ++k) {
                    int ii = k / 8; // row in 8x8 block
                    int jj = k % 8; // column in 8x8 block

                    int q1 = Qm1[jj][ii];  // First quantization step
                    int q2 = Qm2[jj][ii];  // Second quantization step

                    // Check if q1, q2 is a contributing pair
                    if (E[q1-1][q2-1] != 0 && E[q1-1][q2-1] % 2 == 0) {
                        // Check for contributing multiple
                        if (static_cast<int>(D1[i][j][k] - E[q1 - 1][q2 - 1] * 0.5) % E[q1 - 1][q2 - 1] == 0) {
                            // Increment message bit counter
                            ++c;

                            // Update Diff_matrix
                            Diff_matrix[i][j][k] = Energy + 1;
                        }
                    }
                }
            }
        }

        // Create the diff_vector
        int counter = 0;
        for (int k = 0; k < 64; ++k) {
            for (int i = 0; i < MB; ++i) {
                for (int j = 0; j < NB; ++j) {
                    if (Diff_matrix[i][j][k] < 0) {
                        Diff_vector[counter] = Diff_matrix[i][j][k];
                        ++counter;
                    }
                }
            }
        }
    }

    if (method == "dct energy" || method == "pqe" || method == "PQe") {
        int c = 0;  // Message bit counter

        for (int i = 0; i < MB; ++i) {
            for (int j = 0; j < NB; ++j) {
                // Calculate the texture function 'Energy'
                double Energy = 0.0;

                for (int w = 0; w < 64; w++) {
                    Energy += D1[i][j][w] * D1[i][j][w];

                    // Iterate over the 64 DCT mode
                    // assign -energy to all contributing multiples
                    for (int k = 0; k < 64; ++k) {
                        int ii = k / 8; // row in 8x8 block
                        int jj = k % 8; // column in 8x8 block

                        int q1 = Qm1[jj][ii];  // First quantization step
                        int q2 = Qm2[jj][ii];  // Second quantization step

                        // Check if q1, q2 is a contributing pair
                        if (E[q1-1][q2-1] != 0 && E[q1-1][q2-1] % 2 == 0) {
                            // Check for contributing multiple
                            if (static_cast<int>(D1[i][j][k] - E[q1 - 1][q2 - 1] * 0.5) % E[q1 - 1][q2 - 1] == 0) {
                                // Increment message bit counter
                                ++c;

                                // Update Diff_matrix
                                Diff_matrix[i][j][k] = -Energy - 1;
                            }
                        }
                    }
                }
            }

            // Create the diff_vector
            int counter = 0;
            for (int x = 0; x < 64; ++x) {
                for (int y = 0; y < MB; ++y) {
                    for (int z = 0; z < NB; ++z) {
                        if (Diff_matrix[y][z][x] < 0) {
                            Diff_vector[counter] = Diff_matrix[y][z][x];
                            ++counter;
                        }
                    }
                }
            }
        }
    }

    std::vector<int> indices(Cap);

    sortRowsEmulation(Diff_vector, indices);

    std::fill(Diff_vector.begin(), Diff_vector.end(), 0);

    for (size_t i = 0; i < MessageLength; ++i) {
        Diff_vector[indices[i]] = 1.0; // Set the first 'MessageLength' indices to 1.0
    }

    auto D2n = embedMessage(Qm1, Qm2, E, D1, Diff_vector, Message, D2raw, nonzero_spec, MB, NB);

    vector<vector<double>> resultPlane;
    VecToPlane(resultPlane, MB, NB, D2n);

    // Save perturbed coefficients
    save_perturbed_jpeg(resultPlane);

    cout << "\n ---- Perturbed jpeg image saved successfully! ----" << endl;

    /* CHECK STORED COEFFICIENTS BY RELOADING STORED JPEG */

    auto filenameLast = "/home/mato/CLionProjects/StegoDisk/output.jpg";

    FILE *infileLast;
    if ((infileLast = fopen(filenameLast, "rb")) == nullptr) {
        fprintf(stderr, "Can't open %s\n", filenameLast);
        return -1;
    }

    jpeg_decompress_struct cinfoLast{};
    auto dct_coefficientsLast = vector<vector<double>>();

    read_jpeg_info(infileLast, cinfoLast, dct_coefficientsLast);

    writeVectorToFile(dct_coefficientsLast, "/home/mato/Desktop/result_plane.txt");

    return 0;
}
