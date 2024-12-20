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

#pragma region HELPER FUNCTIONS/MACROS

#define DCTSIZE 8
#define DCTSIZE2 64
#define CONST_BITS 13
#define PASS1_BITS 2
#define INT32 long
#define JSAMPLE unsigned char
#define DCTELEM	int
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

// Write 1D vector
template<typename T>
void WriteVectorToFile(const std::vector<T>& vec1D, const std::string& filename) {
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
void WriteVectorToFile(const std::vector<std::vector<T>>& vec2D, const std::string& filename) {
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
void WriteVectorToFile(const std::vector<std::vector<std::vector<T>>>& vec3D, const std::string& filename) {
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

void PrintQuantTable(jpeg_decompress_struct* cinfo) {
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

void PrintMatrix(const vector<vector<int>>& Q) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            cout << Q[i][j] << " "; // Dereference the shared_ptr and access elements
        }
        cout << endl;
    }
}

int RoundToInt(const double val) {
    return static_cast<int>(std::round(val));
}

int Sign(const int val) {
    return (val > 0) - (val < 0);
}

#pragma endregion HELPER FUNCTIONS/MACROS

/**
 * @brief Performs the forward discrete cosine transform (FDCT) on a block of spatial data.
 *
 * This function computes the forward discrete cosine transform (FDCT) for an 8x8 block of spatial data using a slow implementation.
 * It takes the spatial domain values (image block) and applies the DCT, storing the resulting frequency domain coefficients.
 *
 * @param[out] out Output array to store the 8x8 DCT coefficients after the transformation.
 *                 The resulting coefficients are stored in row-major order.
 * @param[in] in Input array containing the 8x8 spatial block to be transformed.
 *               These are typically pixel values of a single block from the image.
 *
 * @note The input `in` is expected to contain the spatial domain pixel values in a flattened form (64 elements).
 * @note The function uses a "slow" implementation of DCT, which may be less efficient but provides a simple method for the transformation.
 */
void PerformSlowFDCT (double out[], const double in[])
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

/**
 * @brief Performs the inverse discrete cosine transform (IDCT) on a block of DCT coefficients.
 *
 * This function computes the inverse DCT for an 8x8 block of DCT coefficients using a slow implementation.
 * It takes the quantized DCT coefficients, applies the inverse DCT transformation, and stores the resulting
 * spatial domain values in the output array.
 *
 * @param[out] out Output array to store the 8x8 spatial block after applying the inverse DCT.
 *                 The resulting values are stored in row-major order.
 * @param[in] coef_block Input array containing the quantized DCT coefficients for the block.
 *                      These are typically 64 DCT coefficients from a single 8x8 block.
 * @param[in] qt Quantization matrix used to scale the DCT coefficients back to their original values before IDCT.
 *              This matrix is used for the inverse quantization process before applying the IDCT.
 *
 * @note The input `coef_block` is expected to contain the quantized DCT coefficients in a flattened form (64 elements).
 * @note The function uses a "slow" implementation of IDCT, which may be less efficient but provides a simple method for inverse transformation.
 */
void PerformSlowIDCT (double out[], const double coef_block[], const double qt[])
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

/**
 * @brief Computes the quantization matrix for a given JPEG quality factor.
 *
 * This function generates an 8x8 quantization matrix `Q` based on the specified quality factor.
 * It scales the default 50% quality quantization matrix (`Q50`) proportionally to the quality value,
 * where higher quality results in lower quantization values (less compression).
 *
 * @param[in] quality JPEG quality factor (range: 1-100).
 *                    - Values closer to 1 result in higher compression and lower image quality.
 *                    - Values closer to 100 result in lower compression and higher image quality.
 * @param[out] Q 8x8 quantization matrix to store the computed values.
 *               - Higher values increase compression by reducing precision in the frequency domain.
 *               - Lower values preserve more detail at the cost of higher file size.
 *
 * @note The input quality value is clamped to the range [1, 100] if it exceeds the bounds.
 * @note The generated matrix is based on scaling the default JPEG 50% quality matrix.
 */
void ComputeQmatrix(int quality, vector<vector<int>>& Q) {
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

/**
 * @brief Extracts DCT coefficients from a JPEG image.
 *
 * This function reads the DCT coefficients from a JPEG image file using libjpeg.
 * The coefficients are organized into a 2D matrix representing the entire image.
 *
 * @param[in] cinfo JPEG decompression structure for reading the file.
 * @param[out] dct_coefficients 2D vector to store the DCT coefficients of the image.
 */
void ExtractDCTCoefficients(j_decompress_ptr cinfo, vector<vector<double>>& dct_coefficients) {
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

/**
 * @brief Reads JPEG file information and extracts its DCT coefficients.
 *
 * This function initializes a JPEG decompression structure, reads the header information
 * of the input JPEG file, and extracts its quantized DCT coefficients into a 2D matrix.
 *
 * @param[in] infile Pointer to the opened JPEG file to be read.
 * @param[out] cinfo JPEG decompression structure used to store the file's metadata and state.
 * @param[out] dct_coefficients 2D vector to store the quantized DCT coefficients from the JPEG image.
 *                              - Rows correspond to image height in pixels.
 *                              - Columns correspond to image width in pixels.
 *
 * @note This function assumes the input file is a valid JPEG image.
 * @note Calls `ExtractDCTCoefficients` to process the coefficient arrays.
 */
void ReadJPEGInfo(FILE *infile, jpeg_decompress_struct& cinfo, vector<vector<double>>& dct_coefficients) {
    struct jpeg_error_mgr jerr{};

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    jpeg_read_header(&cinfo, TRUE);

    ExtractDCTCoefficients(&cinfo, dct_coefficients);
}

/**
 * @brief Initializes contributing quantization pairs.
 *
 * Populates the matrix `E` with contributing quantization pairs (q1, q2)
 * based on the conditions that determine if a pair contributes to the embedding process.
 *
 * @param[out] E 3D matrix to store contributing pairs and their properties.
 * @param[in] MaxQ Maximum quantization value to consider.
 */
void InitializeContributingPairs(vector<vector<int>>& E) {
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

/**
 * @brief Converts a 2D matrix into a 3D DCT coefficient structure.
 *
 * Reshapes a 2D image plane into a 3D vector structure where the first two dimensions represent
 * the block indices, and the third dimension stores DCT coefficients in column-major order.
 *
 * @param[in] plane 2D matrix representing the image.
 * @param[in] Y Height of the image in pixels.
 * @param[in] X Width of the image in pixels.
 * @param[out] D1 3D vector to store the reshaped blocks of DCT coefficients.
 */
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

/**
 * @brief Converts a 3D DCT coefficient structure into a 2D image plane.
 *
 * Reconstructs a 2D matrix from a 3D structure where blocks are organized column-major.
 * Each block contributes an 8x8 region in the output 2D matrix.
 *
 * @param[out] plane 2D matrix to store the reconstructed image.
 * @param[in] MB Number of block rows.
 * @param[in] NB Number of block columns.
 * @param[in] cube 3D vector of DCT coefficients.
 */
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

/**
 * @brief Decompresses a DCT-transformed image back into the spatial domain.
 *
 * This function takes quantized DCT coefficients organized into 8x8 blocks (`D`), applies
 * inverse quantization and inverse discrete cosine transform (IDCT) using the quantization
 * matrix (`Q`), and reconstructs the decompressed image (`X`) in the spatial domain.
 *
 * @param[in] D 3D vector containing quantized DCT coefficients.
 *              - Dimensions: [MD][ND][64], where MD and ND are the number of vertical and horizontal blocks.
 *              - Each block contains 64 coefficients in a flattened 8x8 structure, arranged column-major.
 * @param[in] Q 2D vector representing the quantization matrix (8x8).
 * @param[out] X 2D vector to store the decompressed image in the spatial domain.
 *               - Dimensions: [MD*8][ND*8], where each block contributes an 8x8 region.
 *
 * @note This function assumes the DCT coefficients and quantization matrices are valid.
 * @note Uses the `jpeg_idct_islow` function to perform IDCT with input quantized coefficients and quantization steps.
 */
void DecompressImage(const vector<vector<vector<int>>>& D, const vector<vector<int>>& Q, vector<vector<double>>& X) {
    const int MD = D.size();       // Number of 8x8 blocks in the vertical direction
    const int ND = D[0].size();    // Number of 8x8 blocks in the horizontal direction

    // Temporary block to store the spatial representation of a DCT block
    vector Dblock(8, vector<int>(8));

    // Loop over each block
    for (int i = 0; i < MD; i++) {
        for (int j = 0; j < ND; j++) {
            int B = 8;
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
            PerformSlowIDCT(out_array, coefs_array, quant_array);

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
}

/**
 * @brief Saves a 2D image matrix as a JPEG file.
 *
 * Compresses and writes a 2D image matrix into a JPEG file with a specified quality.
 *
 * @param[in] filename Path to the output JPEG file.
 * @param[in] image 2D image matrix to save.
 * @param[in] quality JPEG compression quality (1-100).
 */
void SaveJPEG(const char* filename, const vector<vector<double>>& image, int quality) {
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

/**
 * @brief Counts the number of non-zero coefficients in a 3D matrix within a specified range.
 *
 * This function iterates through the given 3D matrix and counts the number of non-zero elements
 * within the specified range of the last dimension. The range is defined by the `start` and `end` indices
 * (inclusive).
 *
 * @param[in] matrix Input 3D matrix to search for non-zero coefficients.
 * @param[in] start Starting index for the range in the third dimension (inclusive).
 * @param[in] end Ending index for the range in the third dimension (inclusive).
 * @return The count of non-zero coefficients in the specified range.
 *
 * @note Assumes that `start` and `end` are valid indices within the bounds of the third dimension of the matrix.
 */
int CountNonZeroCoefficients(const vector<vector<vector<int>>>& matrix, const int start, const int end) {
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

/**
 * @brief Computes the capacity and number of non-zero coefficients.
 *
 * Determines the embedding capacity and the total number of non-zero DCT coefficients
 * in the image by evaluating the contributing quantization pairs.
 *
 * @param[in] Qm1 First quantization matrix.
 * @param[in] Qm2 Second quantization matrix.
 * @param[in] D 3D array of DCT coefficients.
 * @param[in] E Matrix of contributing quantization pairs.
 * @param[out] Cap Calculated embedding capacity.
 * @param[out] N0 Count of non-zero DCT coefficients.
 */
void ComputeCapacity(const vector<vector<int>>& Qm1, const vector<vector<int>>& Qm2,
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

/**
 * @brief Converts a spatial 8x8 block into quantized DCT coefficients.
 *
 * Performs a discrete cosine transform (DCT) on a spatial 8x8 block and quantizes the resulting coefficients using the provided quantization matrix.
 *
 * @param[in] Z Input 8x8 spatial block.
 * @param[in] Qf Quantization matrix (8x8).
 * @param[out] QD Output 8x8 block of quantized DCT coefficients.
 * @return Quantized DCT coefficients as an 8x8 block.
 */
void ComputeQuantizedDCT(const vector<vector<double>>& Z, const vector<vector<double>>& Qf, vector<vector<double>>& QD) {
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

    PerformSlowFDCT(out, in);

    // Column-first
    std::vector D(8, vector(8, 0.0));
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
            D[i][j] = out[ 8 * i + j];

    // Quantize the DCT coefficients
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            QD[i][j] = D[i][j] / (Qf[i][j] * 8);
        }
    }
}

/**
 * @brief Calculates the DCT coefficients for an image and organizes them into blocks.
 *
 * This function computes the discrete cosine transform (DCT) coefficients for the input 2D image matrix `X`.
 * The image is divided into 8x8 blocks, and the DCT coefficients for each block are stored in a 3D array `D_raw`.
 * If the dimensions of the image are not multiples of 8, the image is cropped to the nearest 8-pixel boundary.
 *
 * @param[in] X Input 2D matrix representing the image (grayscale intensity values).
 * @param[out] D_raw 3D array to store the calculated DCT coefficients.
 *                  - The first two dimensions represent the block indices.
 *                  - The third dimension contains the 64 coefficients of the corresponding 8x8 DCT block, arranged column-major.
 */
void ComputeDCTBlocks(const vector<vector<double>>& X, vector<vector<vector<double>>>& D_raw) {
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

            auto Qf = std::vector(B, std::vector<double>(B, 1));
            auto Dblock = std::vector(8, std::vector(8, 0.0));

            // Perform DCT (raw2jpg02raw function to be implemented later)
            ComputeQuantizedDCT(Block, Qf, Dblock);

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

/**
 * @brief Generates a random message for embedding.
 *
 * Creates a binary random message with 1's and -1's, generated with a specific bias.
 * The message length is determined by the minimum of the available capacity and the absolute message length.
 *
 * @param[out] Message Vector to store the generated message.
 * @param[in] Cap Maximum embedding capacity.
 * @param[in] AbsoluteMsgLength Desired length of the message to generate.
 * @param[in] seed Seed for random number generation.
 */
void GenerateMessage(vector<int>& Message, const size_t Cap, const size_t AbsoluteMsgLength) {
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

void SortRowsEmulation(vector<double>& diff_vector, vector<int>& indices) {
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

/**
 * @brief Embeds a message into the DCT coefficients.
 *
 * Performs the embedding of a random message into the double-compressed DCT coefficients (D2n)
 * while ensuring the coefficients adhere to specific quantization conditions.
 *
 * @param[in] D1 Original DCT coefficients.
 * @param[in] D2raw Double-compressed raw DCT coefficients.
 * @param[in] Qm1 First quantization matrix.
 * @param[in] Qm2 Second quantization matrix.
 * @param[in] E Matrix of contributing quantization pairs.
 * @param[in] Message The message to embed.
 * @param[out] D2n Resulting DCT coefficients after embedding.
 * @param[out] changes Counter for the number of coefficient changes.
 * @param[out] Zeros Counter for zero-valued coefficients in the result.
 * @param[out] NonZeros Counter for non-zero coefficients in the result.
 * @param[in] nonzero_spec Specification to determine which coefficients to include in calculations.
 */
void EmbedMessage(const vector<vector<int>>& Qm1, const vector<vector<int>>& Qm2,
                  const vector<vector<int>>& E, const vector<vector<vector<int>>>& D1,
                  const vector<double>& Diff_vector, const vector<int>& Message,
                  const vector<vector<vector<double>>>& D2raw, vector<vector<vector<double>>>& D2n,
                  const string& nonzero_spec, const int MB, const int NB) {
    int c = -1;                   // Message bit counter
    int d = -1;                   // Contributing multiples counter
    int changes = 0;             // Number of changes against cover
    int NonZeros = 0, Zeros = 0; // Zero and non-zero counters

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
                            D2n[i][j][k] = (D1[i][j][k] * q1 + Sign(D1[i][j][k]) * Message[c] * q2 / 2) / q2;

                            int cover_value = RoundToInt(D2raw[i][j][k] / q2);
                            if (std::abs(q2 * cover_value - q1 * RoundToInt(q2 * cover_value / q1)) == q2 / 2) {
                                cover_value += Sign(D2raw[i][j][k] / q2 - cover_value);
                            }

                            if (D2n[i][j][k] != cover_value) {
                                changes++;
                            }
                        } else {
                            // D1[i][j][k] is not suitable for embedding
                            D2n[i][j][k] = RoundToInt(D2raw[i][j][k] / q2);
                            if (std::abs(q2 * D2n[i][j][k] - q1 * RoundToInt(q2 * D2n[i][j][k] / q1)) == q2 / 2) {
                                D2n[i][j][k] += Sign(D2raw[i][j][k] / q2 - D2n[i][j][k]);
                            }
                        }
                    } else {
                        // D1[i][j][k] is not a contributing multiple
                        D2n[i][j][k] = RoundToInt(D2raw[i][j][k] / q2);

                        if (std::abs(q2 * D2n[i][j][k] - q1 * RoundToInt(q2 * D2n[i][j][k] / q1)) == q2 / 2) {
                            D2n[i][j][k] += Sign(D2raw[i][j][k] / q2 - D2n[i][j][k]);

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
                    D2n[i][j][k] = RoundToInt(D2raw[i][j][k] / q2);
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
}

/**
 * @brief Saves a perturbed image (after applying steganography or other transformations) as a JPEG file.
 *
 * This function takes a 2D matrix representing the perturbed image in the spatial domain (e.g., after applying DCT, quantization, or embedding)
 * and saves it as a JPEG file. The image is encoded into the JPEG format with a specified compression quality.
 *
 * @param[in] Plane 2D matrix representing the perturbed image in the spatial domain.
 *                  - This matrix is typically a grayscale image, where each element corresponds to a pixel value.
 * @return Returns `true` on successful saving, or `false` if an error occurs during the process.
 *
 * @note The compression quality of the JPEG file can be adjusted by modifying the quality parameter, which is typically set between 1 (lowest quality) and 100 (highest quality).
 */
int SavePerturbedJPEG(vector<vector<double>>& Plane) {
    FILE *infile = fopen("/home/mato/CLionProjects/StegoDisk/cover.jpg", "rb");
    if (!infile) {
        fprintf(stderr, "Can't open %s\n", "/home/mato/CLionProjects/StegoDisk/cover.jpg");
        return -1;
    }

    jpeg_decompress_struct cinfo{};
    jpeg_error_mgr jerr{};
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    jpeg_read_header(&cinfo, TRUE);

    jvirt_barray_ptr* coef_arrays = jpeg_read_coefficients(&cinfo);

    jpeg_compress_struct cinfo_out{};
    jpeg_error_mgr jerr_out{};

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

    cout << "\n ---- Perturbed jpeg image saved successfully! ----" << endl;

    return 0;
}

int main() {
    cout << "\n ---- Algorithm started successfully! ----" << endl;

    const string method = "pqt";
    const string nonzero_spec = "AC-DC";
    constexpr double capacity = 0.1;
    constexpr int Q1 = 85;
    constexpr int Q2 = 70;

    FILE *infile;
    jpeg_decompress_struct cinfo{};
    auto dct_coefficients = vector<vector<double>>();
    auto Qm1 = vector(8, vector<int>(8));
    auto Qm2 = vector(8, vector<int>(8));
    auto E = vector(121, vector<int>(121));
    auto D1 = vector<vector<vector<int>>>();

    //std::srand(12345);

    auto filename = "/home/mato/CLionProjects/StegoDisk/img_1200.jpg";
    if ((infile = fopen(filename, "rb")) == nullptr) {
        fprintf(stderr, "Can't open %s\n", filename);
        return -1;
    }

    ReadJPEGInfo(infile, cinfo, dct_coefficients);

    int image_height = cinfo.image_height;
    int image_width = cinfo.image_width;

    ComputeQmatrix(Q1, Qm1);
    ComputeQmatrix(Q2, Qm2);

    // Print the generated quantization matrix
    //printMatrix(Qm1);
    //printMatrix(Qm2);

    InitializeContributingPairs(E);

    PlaneToVec(dct_coefficients, image_height, image_width, D1);

    auto X1 = vector(image_height, vector(image_width, 0.0));
    DecompressImage(D1, Qm1, X1);

    int MB = floor(X1.size() / 8);
    int NB = floor(X1[0].size() / 8);

    // Write the compressed image into cover directory
    SaveJPEG("/home/mato/CLionProjects/StegoDisk/cover.jpg", X1, Q2);

    fclose(infile);
    jpeg_destroy_decompress(&cinfo);

    //-------------------------------
    // SECOND COMPRESSION STARTS HERE
    //-------------------------------

    size_t AbsoluteMsgLength = 0;
    size_t Cap;
    int N0;

    FILE *infileTwo;
    jpeg_decompress_struct cinfoTwo{};
    auto dct_coefficients_two = vector<vector<double>>();
    auto D1_two = vector<vector<vector<int>>>();
    auto D2raw = vector<vector<vector<double>>>();
    vector<int> Message;
    std::vector D2n(MB, std::vector(NB, std::vector<double>(64)));
    vector<vector<double>> resultPlane;

    if ((infileTwo = fopen("/home/mato/CLionProjects/StegoDisk/cover.jpg", "rb")) == nullptr) {
        fprintf(stderr, "Can't open %s\n", "/home/mato/CLionProjects/StegoDisk/cover.jpg");
        return -1;
    }

    ReadJPEGInfo(infileTwo, cinfoTwo, dct_coefficients_two);

    PlaneToVec(dct_coefficients_two, image_height, image_width, D1_two);

    if (nonzero_spec == "DC-DC") {
        // Include all terms (DC included) in capacity calculations
        int Nz_all = CountNonZeroCoefficients(D1_two, 0, 63); // Count non-zeros in entire matrix
        AbsoluteMsgLength = round(Nz_all * capacity);
    }

    if (nonzero_spec == "AC-DC") {
        // Do not include DC terms in capacity calculations, but use them for embedding
        int Nz_ac = CountNonZeroCoefficients(D1_two, 1, 63); // Count non-zeros in AC terms (ignoring the first DC component)
        AbsoluteMsgLength = round(Nz_ac * capacity);
    }

    if (nonzero_spec == "AC-AC") {
        // Do not include DC terms in capacity calculations and do not embed into them
        int Nz_ac = CountNonZeroCoefficients(D1_two, 1, 63); // Count non-zeros in AC terms
        AbsoluteMsgLength = round(Nz_ac * capacity);

        // Modify the quantization matrix to make DC terms non-contributing
        Qm1[0][0] = 6;
    }

    ComputeCapacity(Qm1, Qm2, D1, E, Cap, N0);

    // Second compression, without quantization
    ComputeDCTBlocks(X1, D2raw);

    // Generate random secret message 1's and -1's of length MessageLength
    size_t MessageLength = Cap > AbsoluteMsgLength ? AbsoluteMsgLength : Cap;
    GenerateMessage(Message, Cap, MessageLength);

    // Initialize Diff_vector and Diff_matrix
    vector Diff_vector(Cap, 0.0);
    vector Diff_matrix(MB, vector(NB, vector(64, 0.0)));

    if (method == "midpoint" || method == "pq" || method == "PQ") {
        size_t c = 0;  // Message bit counter

        // Iterate over the 64 DCT modes
        for (int k = 0; k < 64; ++k) {
            int q1 = Qm1[k % 8][k / 8];  // q1 = first quantization step
            int q2 = Qm2[k % 8][k / 8];  // q2 = second quantization step

            if (E[q1-1][q2-1] != 0 && E[q1-1][q2-1] % 2 == 0) {
                // Loop over all blocks
                for (int i = 0; i < MB; ++i) {
                    for (int j = 0; j < NB; ++j) {
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

        PlaneToVec(X1, image_height, image_width, X1_T);

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

        PlaneToVec(X1, image_height, image_width, X1_T);

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

    vector<int> indices(Cap);
    SortRowsEmulation(Diff_vector, indices);

    std::fill(Diff_vector.begin(), Diff_vector.end(), 0);
    for (size_t i = 0; i < MessageLength; ++i) {
        Diff_vector[indices[i]] = 1.0; // Set the first 'MessageLength' indices to 1.0
    }

    EmbedMessage(Qm1, Qm2, E, D1, Diff_vector, Message, D2raw, D2n, nonzero_spec, MB, NB);
    VecToPlane(resultPlane, MB, NB, D2n);
    SavePerturbedJPEG(resultPlane);

    //-------------------------------
    /* CHECK STORED JPEG */
    //-------------------------------

    FILE *infileLast;
    jpeg_decompress_struct cinfoLast{};
    auto dct_coefficientsLast = vector<vector<double>>();

    auto filenameLast = "/home/mato/CLionProjects/StegoDisk/output.jpg";
    if ((infileLast = fopen(filenameLast, "rb")) == nullptr) {
        fprintf(stderr, "Can't open %s\n", filenameLast);
        return -1;
    }

    ReadJPEGInfo(infileLast, cinfoLast, dct_coefficientsLast);
    WriteVectorToFile(dct_coefficientsLast, "/home/mato/Desktop/result_plane.txt");

    return 0;
}
