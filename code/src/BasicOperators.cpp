/*! \file	BasicOperators.cpp
 *
 *  \brief	Basic operators;
 *
 */


#include "Head.h"
#include "Data.h"
#include "DataOperators.h"
#include "functs.h"


/**
 * \brief	Maximal absolute value;
 */
double normRealInfty ( double *src, int n )
{
    double tmp;
    double rslt = 0.0;
    for (int i = 0; i < n; i++)
    {
        tmp = fabs(src[i]);
        rslt = (rslt > tmp ? rslt : tmp);
    }
    return rslt;
}


/**
 * \brief	Maximal absolute value of an integer matrix;
 */
int fn_max_abs ( int **src, int row, int col )
{
	int tmp, rslt = 0;
	for ( int i = 0; i < row; i++ )
	{
		for ( int j = 0; j < col; j++ )
		{
			tmp = fabs(src[i][j]);
			rslt = ( rslt > tmp ? rslt : tmp );
		}
	}
	return rslt;
}


/**
 * \brief	Define the sort rule;
 */
bool fn_compare (stu_sort a, stu_sort b)
{
	return a.val < b.val;	// increase;
}

/* --------------------		Value (complex)	--------------------*/

/**
 * \brief	Set zero for a complex value;
 */
void fn_complex_setZero ( fftw_complex rslt )
{
	rslt[0] = 0.0; rslt[1] = 0.0;
}


/**
 * \brief	Absolute value for a complex value;
 */
double fn_complex_abs ( fftw_complex src )
{
	double rslt;
	rslt = src[0]*src[0] + src[1]*src[1];
	return sqrt( rslt );
}


/**
 * \brief	complex F1 multiplies complex F2;
 */
void fn_complex_multiply	(	fftw_complex	F1, 
								fftw_complex	F2,
								fftw_complex	rslt )
{
	// (a+ib)(c+id) = (ac-bd) + i(ad+bc);
	rslt[0] = F1[0]*F2[0] - F1[1]*F2[1];
	rslt[1] = F1[0]*F2[1] + F1[1]*F2[0];
}


/**
 * \brief	complex F1 divides complex F2;
 */
void fn_complex_divide		(	fftw_complex	F1, 
								fftw_complex	F2,
								fftw_complex	rslt )
{
	// (a+ib)/(c+id) = (ac+bd)/(c^2+d^2) + i(bc-ad)/(c^2+d^2);
	double F2norm = F2[0]*F2[0] + F2[1]*F2[1];
	rslt[0] = ( F1[0]*F2[0] + F1[1]*F2[1] ) / F2norm;
	rslt[1] = ( F1[1]*F2[0] - F1[0]*F2[1] ) / F2norm;
}

/* --------------------		Value (complex)	--------------------*/

/* --------------------		tvec (complex)	--------------------*/

/**
 * \brief	Set zero of a 'tvec' vector;
 */
void fn_tvec_setZero_complex ( tvec<fftw_complex> src )
{
	for ( int i = 0; i < src.len; i++ )
		src.val[i][0] = 0, src.val[i][1] = 0;
}


/**
 * \brief	Transpose a matrix whose is straightened as a 'tvec' vector;
 *
 * \param	row:	the number of rows of 'src';
 *			col:	the number of columns of 'src';
 
 * \note!	rslt.len == src.len;
 */
void fn_tvec_trans_complex (	tvec<fftw_complex>		rslt,
								tvec<fftw_complex>		src,
										int				row,
										int				col )
{
	for ( int i = 0; i < row; i++ )
	{
		for ( int j = 0; j < col; j++ )
		{
			int rsltInd = i*col + j;
			int srcInd	= j*row + i;
			rslt.val[rsltInd][0] = src.val[srcInd][0];
			rslt.val[rsltInd][1] = src.val[srcInd][1];
		}
	}
}


/**
 * \brief	Maximal absolute value of a complex 'tvec' vector;
 */
double fn_tvec_maxAbs_complex ( tvec<fftw_complex> src )
{
	double tmp;
	double rslt = 0.0;
	for ( int i = 0; i < src.len; i++ )
	{
		tmp = fn_complex_abs ( src.val[i] );
		rslt = (rslt > tmp ? rslt : tmp);
	}
	return rslt;
}


/**
 * \brief	Norm 2 of a complex 'tvec' vector;
 */
double fn_tvec_norm_complex ( tvec<fftw_complex> src )
{
	double rslt = 0.0;
	for ( int i = 0; i < src.len; i++ )
	{
		rslt += pow(src.val[i][0],2) + pow(src.val[i][1],2);
	}
	return sqrt(rslt);
}


/**
 * \brief	The adding result of a complex 'tvec' vector and a double value;
 */
void fn_tvec_constAdd_complex			( tvec<fftw_complex>		rslt,
												double				a )
{
	for ( int i = 0; i < rslt.len; i++ )
	{
		rslt.val[i][0] += a;
		rslt.val[i][1] += a;
	}
}


/**
 * \brief	The multiplying result of a complex 'tvec' vector and a double value;
 */
void fn_tvec_constMultiply_complex		(	tvec<fftw_complex>		src,
												double				a )
{
	for ( int i = 0; i < src.len; i++ )
	{
		src.val[i][0] *= a;
		src.val[i][1] *= a;
	}
}


/**
 * \brief	Result of adding two complex 'tvec' vectors;
 */
void fn_tvec_add_complex (	tvec<fftw_complex>		dst, 
							tvec<fftw_complex>		src, 
								double				dstCoeff, 
								double				srcCoeff )
{
	for ( int i = 0; i < dst.len; i++ )
	{
		dst.val[i][0] = dstCoeff*dst.val[i][0] + srcCoeff*src.val[i][0];
		dst.val[i][1] = dstCoeff*dst.val[i][1] + srcCoeff*src.val[i][1];
	}
}


/**
 * \brief	The dot multiplying result of two complex 'tvec' vectors;
 *
 * \note!	F1.len == F2.len == rslt.len;
 */
void fn_tvec_dotMultiply_complex (	tvec<fftw_complex>		F1,
									tvec<fftw_complex>		F2,
									tvec<fftw_complex>		rslt )
{
	for (int i = 0; i < rslt.len; i++)
	{
		// (a+ib)(c+id) = (ac-bd) + i(ad+bc);
		rslt.val[i][0] = F1.val[i][0]*F2.val[i][0] - F1.val[i][1]*F2.val[i][1];
		rslt.val[i][1] = F1.val[i][0]*F2.val[i][1] + F1.val[i][1]*F2.val[i][0];
	}
}


/**
 * \brief	The summation of the dot multiplying result of two complex 'tvec' vectors;
 * 
 * \note!	F1.len == F2.len;
 */
void fn_tvec_dotMultiplySum_complex	(	tvec<fftw_complex>		F1,
										tvec<fftw_complex>		F2,
											fftw_complex		rslt )
{
	fn_complex_setZero ( rslt );
	for ( int i = 0; i < F1.len; i++ )
	{
		// (a+ib)(c+id) = (ac-bd) + i(ad+bc);
		rslt[0] += F1.val[i][0]*F2.val[i][0] - F1.val[i][1]*F2.val[i][1];
		rslt[1] += F1.val[i][0]*F2.val[i][1] + F1.val[i][1]*F2.val[i][0];
	}
}


/**
 * \brief	Save 'tvec' vector;
 */
void fn_tvec_save_complex ( tvec<fftw_complex> src, char filename[] )
{
	/* open file; */
	FILE *fwPath = fopen(filename, "w");
	if ( fwPath == NULL )
		printf("%s error.\n", filename);

	/* save data; */
	if ( src.row == 1 && src.col == 1 )
		fprintf(fwPath, "%d\n", src.len);
	else
		fprintf(fwPath, "%d\t%d\t%d\n", src.row, src.col, src.len);
	for ( int i = 0; i < src.len; i++ )
		fprintf(fwPath, "%+.15E\t%+.15E\n", src.val[i][0], src.val[i][1]);
	fclose(fwPath);
}


/**
 * \brief	Print 'tvec' vector;
 */
void fn_tvec_print_complex ( tvec<fftw_complex> src )
{
	/* print data; */
	if ( src.row == 1 && src.col == 1 )
		printf("len = %d\n", src.len);
	else
		printf("nrow = %d, ncol = %d, nnz = %d\n", src.row, src.col, src.len);
	for ( int i = 0; i < src.len; i++ )
		printf("%+.10E\t%+.10E\n", src.val[i][0], src.val[i][1]);
}

/* --------------------		tvec (complex)	--------------------*/

/* --------------------		tmat (complex)	--------------------*/

/**
 * \brief	Set zero of 'tmat' matrix;
 */
void fn_tmat_setZero_complex ( tmat<fftw_complex> src )
{
	for ( int i = 0; i < src.row; i++ )
		for ( int j = 0; j < src.col; j++ )
			src.val[i][j][0] = 0, src.val[i][j][1] = 0;
}


/**
 * \brief	Save 'tmat' matrix;
 */
void fn_tmat_save_complex ( tmat<fftw_complex> src, char filename[] )
{
	/* open file; */
	FILE *fwPath = fopen(filename, "w");
	if ( fwPath == NULL )
		printf("%s error.\n", filename);

	/* save data; */
	fprintf(fwPath, "%d\t%d\t%d\n", src.row, src.col, src.ele);
	for ( int i = 0; i < src.row; i++ )
		for ( int j = 0; j < src.col; j++ )
			fprintf(fwPath, "%+.15E\t%+.15E\n", 
					src.val[i][j][0], src.val[i][j][1]);
	fclose(fwPath);
}


/**
 * \brief	Print 'tmat' matrix;
 */
void fn_tmat_print_complex ( tmat<fftw_complex> src )
{
	/* print data; */
	printf("nrow = %d, ncol = %d, nnz = %d\n", src.row, src.col, src.ele);
	for ( int i = 0; i < src.row; i++ )
		for ( int j = 0; j < src.col; j++ )
			printf("%+.10E\t%+.10E\n", src.val[i][j][0], src.val[i][j][1]);
}

/* --------------------		tmat (complex)	--------------------*/
