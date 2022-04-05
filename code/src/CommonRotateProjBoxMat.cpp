/*! \file	CommonRotateProjBoxMat.cpp
 *
 * \brief	Generate the common 'rotateProjBoxMat';
 */

#include "Data.h"
#include "Head.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"


/**
 * \brief	Calculate the common 'rotateProjBoxMat';
 *
 * \param	sbulkparam1:		structure body for the left   'rotateProjBoxMat';
 * \param	sbulkparam2:		structure body for the right  'rotateProjBoxMat';
 * \param	ssysparam:			structure body for the common 'rotateProjBoxMat';
 *
 * \note!	here all 'rotateProjBoxMat' not consider the first row;
 *			the first row is corresponding to x (GJP);
 */
int fn_obt_commom_rotateProjBoxMat  (	stu_bulk_param		*sbulkparam1,
										stu_bulk_param		*sbulkparam2,
										stu_system_param	*ssysparam )
{
	if ( strcmp ( ssysparam->com_projmat_way, "direct" ) == 0 )
	{
		/* parameter file; */
		char paraFile[FILELEN];
		sprintf(paraFile, "./para/%s/input%s.dat", paraDir, para_flag);
		printf("parameter file: %s.\n\n", paraFile);
		fn_obt_com_projmat ( ssysparam, paraFile );
	}
	else
	{
		printf(" <========== Calculation of the common projection matrix ==========> \n\n");
		mytimer_t timer;
		timer.reset();
		timer.start();

		double tol  = 1e-6;
		int int_reg = ssysparam->searchReg;

		if ( sbulkparam1->dimPhy != sbulkparam2->dimPhy )
		{
			printf("Error using 'fn_obt_commom_rotateProjBoxMat'\n");
			printf("'dimPhy' of the left/right bulk phases must be equal.\n");
			return 0;
		}
		int nr = sbulkparam1->dimPhy - 1;

		/* matrix of the left bulk phase; */
		tmat<double> xmat;
		xmat.row = nr;
		xmat.col = sbulkparam1->dimCpt;
		xmat.val = (double **) malloc( sizeof(double *) * xmat.col );
		for ( int i = 0; i < xmat.col; i++ )
			xmat.val[i] = (double *) malloc( sizeof(double) * xmat.row );
		for ( int i = 0; i < xmat.col; i++ )
			for ( int j = 0; j < xmat.row; j++ )
				xmat.val[i][j] = sbulkparam1->rotateProjBoxMat.val[j+1][i];

		/* matrix of the right bulk phase; */
		tmat<double> ymat;
		ymat.row = nr;
		ymat.col = sbulkparam2->dimCpt;
		ymat.val = (double **) malloc( sizeof(double *) * ymat.col );
		for ( int i = 0; i < ymat.col; i++ )
			ymat.val[i] = (double *) malloc( sizeof(double) * ymat.row );
		for ( int i = 0; i < ymat.col; i++ )
			for ( int j = 0; j < ymat.row; j++ )
				ymat.val[i][j] = sbulkparam2->rotateProjBoxMat.val[j+1][i];


		/* connect xmat and ymat; */
		tmat<double> xymat = fn_matrix_connect ( xmat, ymat );
		printf("xymat:\n");
		for ( int i = 0; i < xymat.row; i++ )
		{
			for ( int j = 0; j < xymat.col; j++ )
				printf("% .15f\t", xymat.val[j][i]);
			printf("\n");
		}
		printf("\n");

		/* calculate the common 'rotateProjBoxMat'; */
		tmat<double> rsltmat = fn_matrix_rational_dependent ( xymat, tol, int_reg );

		/* save results; */
		ssysparam->dimRePhy	= rsltmat.row;
		ssysparam->dimReCpt	= rsltmat.col;
		ssysparam->rotateProjBoxMat = fn_tmat_init<double> ( 
									ssysparam->dimRePhy, ssysparam->dimReCpt );
		for ( int i = 0; i < ssysparam->dimRePhy; i++ )	// transpose;
			for ( int j = 0; j < ssysparam->dimReCpt; j++ )
				ssysparam->rotateProjBoxMat.val[i][j] = rsltmat.val[j][i];

		/* release memory; */
		for ( int i = 0; i < xmat.col; i++ )
			free(xmat.val[i]);
		free(xmat.val);
		for ( int i = 0; i < ymat.col; i++ )
			free(ymat.val[i]);
		free(ymat.val);
		for ( int i = 0; i < xymat.col; i++ )
			free(xymat.val[i]);
		free(xymat.val);
		for ( int i = 0; i < rsltmat.col; i++ )
			free(rsltmat.val[i]);
		free(rsltmat.val);

		timer.pause();
		printf("\t ***** time cost of calculation of the common projection matrix: ");
		printf("%f seconds *****\n\n", timer.get_current_time());
	}

	printf("\t ---> The common projection matrix (%s) is \n", ssysparam->com_projmat_way);
	fn_matrix_print ( ssysparam->rotateProjBoxMat );
	printf("\n");

	/* create fold; */
	/*
	mkdir("./result/", 0755);
	mkdir(rsltDir, 0755);
	char fname[FILELEN];
	sprintf(fname, "%s/parameter_opt.dat", rsltDir);	// file for saving parameters;
	*/

	/* save parameters; */
	/*
	FILE *fp = fopen(fname, "a");
	fprintf(fp, "\n\n#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n");
	fprintf(fp, "#		the way to get the common projection matrix		#\n");
	fprintf(fp, "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#\n\n");
	fprintf(fp, "# com_projmat_way = calculate: calculate the common projection matrix;\n");
	fprintf(fp, "# com_projmat_way = direct:    directly input the common projection matrix;\n\n");
	fprintf(fp, "com_projmat_way\t\t= %s\n", ssysparam->com_projmat_way);
	fprintf(fp, "com_projmat_size\t= %d\t%d\n", ssysparam->dimRePhy, ssysparam->dimReCpt);
	fprintf(fp, "com_projmat_mat\t\t= \n");
	fn_matrix_save ( ssysparam->rotateProjBoxMat, fp );
	fclose(fp);
	*/

	return 1;
}


/**
 * \brief	Connect 'rotateProjBoxMat' of two bulk phases;
 *			delete the zero vector;
 */
tmat<double> fn_matrix_connect	(	tmat<double>	xmat, 
									tmat<double>	ymat )
{
	int nr = xmat.row; // or ymat.row;

	/* initialize result; */
	tmat<double> xymat;
	xymat.row = nr;
	xymat.col = xmat.col + ymat.col;
	xymat.val = (double **) malloc( sizeof(double *) * xymat.col );
	for ( int i = 0; i < xymat.col; i++ )
		xymat.val[i] = (double *) malloc( sizeof(double) * nr );

	/* connect; */
	int ind0 = 0;
	for ( int i = 0; i < xmat.col; i++ )	// xmat part;
	{
		if ( normRealInfty(xmat.val[i], nr) > ZEROTOL )	// start position;
		{
			for ( int j = 0; j < nr; j++ )
				xymat.val[ind0][j] = xmat.val[i][j];
			ind0 ++ ;
		}
	}
	for ( int i = 0; i < ymat.col; i++ )	// ymat part;
	{
		if ( normRealInfty(ymat.val[i], nr) > ZEROTOL )	// start position;
		{
			for ( int j = 0; j < nr; j++ )
				xymat.val[ind0][j] = ymat.val[i][j];
			ind0 ++ ;
		}
	}
	for ( int i = ind0; i < xymat.col; i++ ) // release the extra memory;
		free(xymat.val[i]);
	xymat.col = ind0;						 // update the number of columns;

	return xymat;
}


/**
 * \brief	Numerically check whether two vectors are linear dependent
 *				over the rational number domain;
 *
 * \param	xmat, ymat:		two real matrices; they have the same number of rows;
 *							have been transposed;
 *
 * \return	rsltmat:		return a real matrix;
 */
tmat<double> fn_matrix_rational_dependent (	tmat<double>	xymat, 
												double		tol, 
												int			int_reg )
{
	int nr = xymat.row;	// invariant; or ymat.row;
//	bool isPrint = true;
	bool isPrint = false;

	/* memory allocation for result; */
	tmat<double> rsltmat;
	rsltmat.row = nr;
	rsltmat.col = xymat.col;
	rsltmat.val = (double **) malloc( sizeof(double *) * rsltmat.col );
	for ( int i = 0; i < rsltmat.col; i++ )
		rsltmat.val[i] = (double *) malloc( sizeof(double) * nr );
	for ( int i = 0; i < rsltmat.col; i++ )	// initialization;
		for ( int j = 0; j < nr; j++ )
			rsltmat.val[i][j] = 0.0;

	/* memory allocation for the integer coefficients; */
	tmat<int> coeffmat;
	coeffmat.row = xymat.col;
	coeffmat.col = rsltmat.col;
	coeffmat.val = (int **) malloc( sizeof(int *) * coeffmat.col );
	for ( int i = 0; i < coeffmat.col; i++ )
		coeffmat.val[i] = (int *) malloc( sizeof(int) * coeffmat.row );
	for ( int i = 0; i < coeffmat.col; i++ )
		for ( int j = 0; j < coeffmat.row; j++ )
			coeffmat.val[i][j] = 0.0;
	if ( isPrint )
	{
		printf("Memory allocation of coeffmat: \n");
		for ( int i = 0; i < coeffmat.col; i++ )
		{
			for ( int j = 0; j < coeffmat.row; j++ )
				printf("% d\t", coeffmat.val[i][j]);
			printf("\n");
		}
		printf("\n");
	}
	
	/* -------------------------------------------------------------------- */
	/**
	 * find the maximal linearly independent group over the rational number domain;
	 */
	int xyIntCoeff, dpdflag;
	double xyTmp, err, errTmp;
	int intSpaceFlag = 0;					// 0: create 'intSpace';
	int rsltLen = 1;
	for ( int j = 0; j < nr; j++ )			// select the first column vector;
		rsltmat.val[0][j] = xymat.val[0][j];
	coeffmat.val[0][0] = 1;

	/**
	 * check each column vector of 'xymat';
	 *	here is row vector since 'xymat' has been transposed;
	 *	transposing matrix is more convenient;
	 */
	int cplxDofs, comDiv;
	int **intSpace;
	for ( int xyInd = 1; xyInd < xymat.col; xyInd++ )	// the first one has been used;
	{
		if ( intSpaceFlag == 0 )
		{
			cplxDofs = pow(int_reg, rsltLen);
			intSpace = (int **) malloc( sizeof(int *) * cplxDofs );
			for ( int i = 0; i < cplxDofs; i++ )
				intSpace[i] = (int *) malloc( sizeof(int) * rsltLen );
			/* generate all cases of integer coefficients; */
			fn_obt_int_coeff (intSpace, rsltLen, int_reg, cplxDofs);
		}

		/**
		 * find a nonzero element in the 'xyInd'-th column vector of 'xymat';
		 * it should exist since we have deleted all zero column vectors;
		 */
		int nonzeroInd = 0;
		for ( int j = 0; j < nr; j++ )
		{
			if ( fabs(xymat.val[xyInd][j]) > tol )
			{
				nonzeroInd = j;
				break;
			}
		}

		/* search representation; */
		dpdflag = 0;	// the flag to check whether they are linearly independent;
		for ( int j0 = 0; j0 < cplxDofs; j0++ )
		{
			/* calculate the coefficient by one element; */
			xyTmp = 0.0;
			for ( int i = 0; i < rsltLen; i++ )	// 'rsltLen' is the dimensionality of 'intSpace';
				xyTmp += intSpace[j0][i] * rsltmat.val[i][nonzeroInd];
			xyIntCoeff = round( xyTmp / xymat.val[xyInd][nonzeroInd] );
			/* xyIntCoeff cannot be 0 and cannot be greater than int_reg; */
			if ( xyIntCoeff == 0 || fabs(xyIntCoeff) > int_reg ) continue;

			/* calculate the representation error; */
			err = 0.0;
			for ( int j = 0; j < nr; j++ )
			{
				errTmp = 0.0;
				for ( int i = 0; i < rsltLen; i++ )
					errTmp += intSpace[j0][i] * rsltmat.val[i][j];
				errTmp = fabs( errTmp - xyIntCoeff*xymat.val[xyInd][j] );
				err = ( err > errTmp ? err : errTmp );
			}

			if ( err < tol )
			{
				/* linearly dependent; */
				for ( int i = 0; i < rsltLen; i++ )
				{
					coeffmat.val[i][xyInd] = intSpace[j0][i];		 // save coefficients;
					if ( intSpace[j0][i] != 0 )
					{
						comDiv = __gcd(xyIntCoeff, intSpace[j0][i]); // common divisor;
						for ( int j = 0; j < rsltLen; j++ )
							coeffmat.val[i][j] *= (xyIntCoeff/comDiv); // update the calculate 'rsltmat';
						coeffmat.val[i][xyInd] /= comDiv;			 // update coefficients;
						for ( int j = 0; j < nr; j++ )
						{
							rsltmat.val[i][j] = rsltmat.val[i][j] / (xyIntCoeff/comDiv);
						}
					}
				}
				dpdflag = 1;
				break;
			}
		}

		if ( isPrint )
		{
			printf("xyInd = %d \t coeffmat: \n", xyInd);
			for ( int i = 0; i < coeffmat.col; i++ )
			{
				for ( int j = 0; j < coeffmat.row; j++ )
					printf("% d\t", coeffmat.val[i][j]);
				printf("\n");
			}
			printf("\n");
		}
		
		/* add the linearly independent column vector; */
		if ( dpdflag == 0 )
		{
			coeffmat.val[rsltLen][xyInd] = 1;
			for ( int j = 0; j < nr; j++ )
				rsltmat.val[rsltLen][j] = xymat.val[xyInd][j];
			rsltLen ++ ;
		}

		/* release memory of 'intSpace'; */
		if ( dpdflag == 0 || xyInd == xymat.col-1 )
		{
			for ( int i = 0; i < cplxDofs; i++ )
				free(intSpace[i]);
			free(intSpace);
			intSpaceFlag = 0;
		}
		else
			intSpaceFlag = 1;
	}

	/* release the extra memory; */
	for ( int i = rsltLen; i < rsltmat.col; i++ )
		free(rsltmat.val[i]);
	rsltmat.col = rsltLen;
	for ( int i = rsltLen; i < coeffmat.col; i++ )
		free(coeffmat.val[i]);
	coeffmat.col = rsltLen;

	/* the maximal absolute value as the check criterion; */
	int check = fn_max_abs ( coeffmat.val, coeffmat.col, coeffmat.row );

	if ( isPrint )
	{
		printf("before check: \n");
		for ( int i = 0; i < coeffmat.col; i++ )
		{
			for ( int j = 0; j < coeffmat.row; j++ )
				printf("% d\t", coeffmat.val[i][j]);
			printf("\n");
		}
		printf("\n");
		printf("-- value: \n");
		for ( int i = 0; i < coeffmat.row; i++ )
		{
			for ( int j = 0; j < rsltmat.row; j++ )
			{
				double chetmp = 0.0;
				for ( int j0 = 0; j0 < coeffmat.col; j0++ )
				{
					chetmp += coeffmat.val[j0][i] * rsltmat.val[j0][j];
				}
				printf("[%d,%d] % .15f\t", i, j, chetmp);
			}
			printf("\n");
		}
		printf("-- rsltmat: \n");
		for ( int i = 0; i < rsltmat.row; i++ )
		{
			for ( int j = 0; j < rsltmat.col; j++ )
				printf("[%d,%d] % .15f\t", i, j, rsltmat.val[j][i]);
			printf("\n");
		}
		printf("\n\n");
	}

	if ( check > 2 && rsltLen > 1 )
	{
		/* optimize the integer coefficients; */
		int opt_int_reg = check;
		int dim = coeffmat.col - 1;
		cplxDofs = pow(opt_int_reg, dim);
		intSpace = (int **) malloc( sizeof(int *) * cplxDofs );
		for ( int i = 0; i < cplxDofs; i++ )
			intSpace[i] = (int *) malloc( sizeof(int) * dim );
		/* generate all cases of integer coefficients; */
		fn_opt_int_coeff (intSpace, dim, opt_int_reg, cplxDofs);

		/* try different cases; */
		int checkTmp;
		int **coeffmatTmp = (int **) malloc( sizeof(int *) * coeffmat.col );
		for ( int i = 0; i < coeffmat.col; i++ )
			coeffmatTmp[i] = (int *) malloc( sizeof(int) * coeffmat.row );
		int **coeffOptMat = (int **) malloc( sizeof(int *) * coeffmat.col );
		for ( int i = 0; i < coeffmat.col; i++ )
			coeffOptMat[i] = (int *) malloc( sizeof(int) * coeffmat.row );
		for ( int i = 0; i < coeffmat.col; i++ )
			for ( int j = 0; j < coeffmat.row; j++ )
				coeffOptMat[i][j] = coeffmat.val[i][j];
		double **rsltOptMat = (double **) malloc( sizeof(double *) * rsltmat.col );
		for ( int i = 0; i < rsltmat.col; i++ )
			rsltOptMat[i] = (double *) malloc( sizeof(double) * rsltmat.row );
		for ( int i = 0; i < rsltmat.col; i++ )
			for ( int j = 0; j < rsltmat.row; j++ )
				rsltOptMat[i][j] = rsltmat.val[i][j];
		bool isbreakloop = true;

		/** 
		 * optimize each column of 'coeffmat';
		 * 'coeffmat.col' == 'rsltmat.col';
		 */
		for ( int ci = 0; ci < coeffmat.col && isbreakloop; ci++ )
		{
			/* try; */
			for ( int i = 0; i < cplxDofs && isbreakloop; i++ )
			{
				/* reset; */
				for ( int i = 0; i < coeffmat.col; i++ )
					for ( int j = 0; j < coeffmat.row; j++ )
						coeffmatTmp[i][j] = coeffmat.val[i][j];
				comDiv = 1;
				for ( int j = 0; j < coeffmat.row; j++ ) // check all rows of 'xymat';
				{
					coeffmatTmp[ci][j] = coeffmat.val[ci][j];
					/* except ci; */
					for ( int j0 = 0; j0 < ci; j0++ )
						coeffmatTmp[ci][j] += intSpace[i][j0] * coeffmat.val[j0][j];
					for ( int j0 = ci+1; j0 < coeffmat.col; j0++ )
						coeffmatTmp[ci][j] += intSpace[i][j0-1] * coeffmat.val[j0][j];

					if ( j == 0 )
						comDiv = coeffmatTmp[ci][0];
					else
						comDiv = __gcd(comDiv, coeffmatTmp[ci][j]); // common divisor;
				}
				for ( int j = 0; j < coeffmat.row; j++ )
					coeffmatTmp[ci][j] /= comDiv;
				checkTmp = fn_max_abs ( coeffmatTmp, coeffmat.col, coeffmat.row );
				if ( checkTmp < check )		// better;
				{
					check = checkTmp;
					/* update integer coefficients; */
					for ( int j = 0; j < coeffmat.row; j++ )
						coeffOptMat[ci][j] = coeffmatTmp[ci][j];

					/** 
					 * update 'rsltmat.val';
					 *	'coeffmat.val' borrow, but 'rsltmat.val' return;
					 */
					for ( int j = 0; j < rsltmat.row; j++ )
					{
						for ( int jj = 0; jj < rsltmat.col; jj++ ) // reset;
							rsltOptMat[jj][j] = rsltmat.val[jj][j];
						for ( int j0 = 0; j0 < ci; j0++ )
							rsltOptMat[j0][j] = intSpace[i][j0] * rsltOptMat[ci][j];
						for ( int j0 = ci+1; j0 < rsltmat.col; j0++ )
							rsltOptMat[j0][j] -= intSpace[i][j0-1] * rsltOptMat[ci][j];
						rsltOptMat[ci][j] *= comDiv;
					}

					if ( isPrint )
					{
						printf("check: %d\n", check);
						printf("ci = %d\n", ci);
						for ( int jj = 0; jj < rsltmat.col-1; jj++ )
							printf("% d\t", intSpace[i][jj]);
						printf("\n");
						printf("comDiv = %d\n", comDiv);
						printf("-- coeffmat: \n");
						for ( int i = 0; i < coeffmat.col; i++ )
						{
							for ( int j = 0; j < coeffmat.row; j++ )
								printf("% d\t", coeffmatTmp[i][j]);
							printf("\n");
						}
						printf("\n");
						printf("-- value: \n");
						for ( int i = 0; i < coeffmat.row; i++ )
						{
							for ( int j = 0; j < rsltmat.row; j++ )
							{
								double chetmp = 0.0;
								for ( int j0 = 0; j0 < coeffmat.col; j0++ )
								{
									chetmp += coeffmatTmp[j0][i] * rsltOptMat[j0][j];
								}
								printf("[%d,%d] % .15f\t", i, j, chetmp);
							}
							printf("\n");
						}
						printf("-- rsltmat: \n");
						for ( int i = 0; i < rsltmat.row; i++ )
						{
							for ( int j = 0; j < rsltmat.col; j++ )
								printf("[%d,%d] % .15f\t", i, j, rsltOptMat[j][i]);
							printf("\n");
						}
						printf("\n\n");
					}
				}
				if ( check <= 2 )			// 2 is very good;
					isbreakloop = false;
			}
		}

		/* update 'coeffmat.val' and 'rsltmat.val'; */
		for ( int i = 0; i < coeffmat.col; i++ )
			for ( int j = 0; j < coeffmat.row; j++ )
				coeffmat.val[i][j] = coeffOptMat[i][j];
		for ( int i = 0; i < rsltmat.col; i++ )
			for ( int j = 0; j < rsltmat.row; j++ )
				rsltmat.val[i][j] = rsltOptMat[i][j];

		/* release memory; */
		for ( int i = 0; i < coeffmat.col; i++ )
			free(coeffmatTmp[i]);
		free(coeffmatTmp);
		for ( int i = 0; i < coeffmat.col; i++ )
			free(coeffOptMat[i]);
		free(coeffOptMat);
		for ( int i = 0; i < rsltmat.col; i++ )
			free(rsltOptMat[i]);
		free(rsltOptMat);
		for ( int i = 0; i < cplxDofs; i++ )
			free(intSpace[i]);
		free(intSpace);

		if ( isPrint )
		{
			printf("after check: \n");
			printf("-- coeffmat: \n");
			for ( int i = 0; i < coeffmat.col; i++ )
			{
				for ( int j = 0; j < coeffmat.row; j++ )
					printf("% d\t", coeffmat.val[i][j]);
				printf("\n");
			}
			printf("\n");
			printf("-- value: \n");
			for ( int i = 0; i < coeffmat.row; i++ )
			{
				for ( int j = 0; j < rsltmat.row; j++ )
				{
					double chetmp = 0.0;
					for ( int j0 = 0; j0 < coeffmat.col; j0++ )
					{
						chetmp += coeffmat.val[j0][i] * rsltmat.val[j0][j];
					}
					printf("[%d,%d] % .15f\t", i, j, chetmp);
				}
				printf("\n");
			}
			printf("-- rsltmat: \n");
			for ( int i = 0; i < rsltmat.row; i++ )
			{
				for ( int j = 0; j < rsltmat.col; j++ )
					printf("[%d,%d] % .15f\t", i, j, rsltmat.val[j][i]);
				printf("\n");
			}
			printf("\n\n");
		}
	}

	/* print result; */
	err = 0.0;
	for ( int i = 0; i < coeffmat.row; i++ )
	{
		for ( int j = 0; j < rsltmat.row; j++ )
		{
			double chetmp = 0.0;
			for ( int j0 = 0; j0 < coeffmat.col; j0++ )
			{
				chetmp += coeffmat.val[j0][i] * rsltmat.val[j0][j];
			}
			chetmp -= xymat.val[i][j];
			chetmp  = fabs(chetmp);
			err = ( err > chetmp ? err : chetmp );
		}
	}
	printf("\t Common projection matrix: \n");
	printf("\t ---> the maximal absolute coefficient is %d\n", check);
	printf("\t ---> the error between the original matrix and ");
	printf("the new representation is % .5e\n", err);

	/* release memory; */
	for ( int i = 0; i < coeffmat.col; i++ )
		free(coeffmat.val[i]);
	free(coeffmat.val);

	return rsltmat;
}


/**
 * \brief	Generate all possible cases of the integer coefficients;
 *
 * \param	intSpace:		return value; all possible cases;
 * \param	dim:		dimensionality;
 * \param	N:			the value range; 1,2,...,N;
 */
void fn_obt_int_coeff (int **intSpace, int dim, int N, int cplxDofs)
{
	int *k = (int *)malloc( sizeof(int) * dim);
	for ( int i = 0; i < dim; i++ ) k[i] = -N/2+1;

	for ( int i = 0; i < cplxDofs; i++ )
	{
		for ( int j = 0; j < dim; j++ )
			intSpace[i][j] = k[j];
		k[dim-1] ++;
		for ( int jj = dim-1; jj > 0; jj-- )
		{
			if (k[jj] > N/2)
			{
				k[jj] = -N/2+1;
				k[jj-1] ++;
			}
		}
	}
	free(k);

	/* sort 'intSpace' by modules; */
	double tmp;
	vector <stu_sort> sort_array(cplxDofs);
	for ( int i = 0; i < cplxDofs; ++i )
	{
		tmp = 0.0;
		for ( int j = 0; j < dim; j++ )
			tmp += pow(intSpace[i][j], 2);
		sort_array[i].ind = i;
		sort_array[i].val = tmp;
	}
	sort(sort_array.begin(), sort_array.end(), fn_compare);

	/* obtain the sorted 'intSpace'; */
	int **rslt = (int **) malloc( sizeof(int *) * cplxDofs );
	for ( int i = 0; i < cplxDofs; i++ )
		rslt[i] = (int *) malloc( sizeof(int) * dim );
	for ( int i = 0; i < cplxDofs; i++ )
		for ( int j = 0; j < dim; j++ )
			rslt[i][j] = intSpace[sort_array[i].ind][j];
	for ( int i = 0; i < cplxDofs; i++ ) // copy;
		for ( int j = 0; j < dim; j++ )
			intSpace[i][j] = rslt[i][j];

	/* release memory; */
	for ( int i = 0; i < cplxDofs; i++ )
		free(rslt[i]);
	free(rslt);
}


/**
 * \brief	Generate all possible cases to optimize the integer coefficients;
 *
 * \param	intSpace:		return value; all possible cases;
 * \param	dim:		dimensionality;
 * \param	N:			the value range; 1,2,...,N;
 */
void fn_opt_int_coeff (int **intSpace, int dim, int N, int cplxDofs)
{
	int *k = (int *)malloc( sizeof(int) * dim);
	for ( int i = 0; i < dim; i++ ) k[i] = -N/2;

	for ( int i = 0; i < cplxDofs; i++ )
	{
		for ( int j = 0; j < dim; j++ )
		{
			if ( fabs(k[j]) == 0 )
				k[j] ++;			// skip 0;
			intSpace[i][j] = k[j];
		}
		k[dim-1] ++;
		for ( int jj = dim-1; jj > 0; jj-- )
		{
			if (k[jj] > N/2)
			{
				k[jj] = -N/2;
				k[jj-1] ++;
			}
		}
	}
	free(k);
}


/**
 * \brief	Numerically check whether two vectors are linear dependent
 *				over the rational number domain;
 * \param	x, y:			two real vectors (nonzero);
 * \param	z:				return a real vector;
 * \param	n:				the length of x,y;
 *
 * \return	an integer:		0: x,y are nonzero and independent;
 *							1: x,y are nonzero and dependent, select x;
 *							2: x,y are nonzero and dependent, select y;
 *							-1: |y| = 0, x is nonzero;
 *							-2: |x| = 0, y is nonzero;
 *							-3: |x| = |y| = 0;
 */
int fn_check_vector_rational_dependent (double *x,	double *y,	double *z,	
										int		n,	double tol,	int	   int_reg)
{
	if ( normRealInfty(x, n) < tol )
	{
		if ( normRealInfty(y, n) < tol )	// |x| = |y| = 0;
			return -3;
		else
			return -2;						// |x| = 0, y is nonzero;
	}
	else if ( normRealInfty(y, n) < tol )
		return -1;							// |y| = 0, x is nonzero;

	/* check zero; */
	for ( int i = 0; i < n; i++ )
	{
		if ( ( fabs(x[i]) < tol && fabs(y[i]) > tol ) ||
		   ( fabs(x[i]) > tol && fabs(y[i]) < tol ) )
			   return 0;					// x[i], y[i]: only one zero;
	}

	/* check that x[i] and y[i] both are nonzero; */
	int *ind = (int *) malloc( sizeof(int) * n );
	int len = 0;
	for ( int i = 0; i < n; i++ )
	{
		if ( fabs(x[i]) > tol && fabs(y[i]) > tol ) // x[i], y[i]: both nonzero;
		{
			ind[len] = i;
			len++;
		}
	}

	/** 
	 * check whether x,y are linearly dependent over the rational number domain;
	 *		x, y are nonzero;
	 *	 dependent: ax=by, a,b are integers belonging to [1, int_reg);
	 * independent: cannot find two integers a,b belonging to [1, int_reg) 
	 *				to satisfy ax=by;
	 */
	double b, err, tmp;
	int	*intb = (int *) malloc( sizeof(int) * len );
	int  intb0;
	for ( int i = 0; i < n; i++ )
		z[i] = 0.0;
	for ( int a = 1; a < int_reg; a++ )
	{
		err = 0.0;
		for ( int i = 0; i < len; i++ )	// search b;
		{
			b		= a * x[ind[i]] / y[ind[i]];
			intb[i] = round(b);
			tmp		= fabs( b - intb[i] );
			err		= (err > tmp ? err : tmp);
		}
		intb0 = intb[0];
		for ( int i = 1; i < len; i++ )	// all elements of b should be equal;
		{
			if ( fabs(intb0) < tol || fabs(intb[i] - intb0) > tol )
			{
				err = 1.0;
				break;
			}
		}
		if ( err < tol )
		{
			if ( fabs(a) > fabs(intb0) )
			{
				for ( int i = 0; i < len; i++ ) // return x;
					z[ind[i]] = x[ind[i]] / intb0;
				free(ind);	free(intb);
				return 1;		// x,y are nonzero and dependent, select x;
			}
			else
			{
				for ( int i = 0; i < len; i++ ) // return y;
					z[ind[i]] = y[ind[i]] / a;
				free(ind);	free(intb);
				return 2;		// x,y are nonzero and dependent, select y;
			}
		}
	}
	free(ind);	free(intb);
	return 0;					// x,y are nonzero and independent;
}


/**
 * \brief	Numerically check whether two real numbers are linearly dependent 
 *				over the rational number domain;
 *
 * \param	x, y:			two real numbers;
 *
 * \return	an integer:		0: x,y are nonzero and independent;
 *							1: x,y are nonzero and dependent, select x;
 *							2: x,y are nonzero and dependent, select y;
 *							-1: y = 0, x is nonzero;
 *							-2: x = 0, y is nonzero;
 *							-3: x = y = 0;
 *							
 */
int fn_check_number_rational_dependent (double x,	double y,
									    double tol,	int	   int_reg)
{
	if ( fabs(x) < tol )
	{
		if ( fabs(y) < tol )	// x = y = 0;
			return -3;
		else
			return -2;			// x = 0, y is nonzero;
	}
	else if ( fabs(y) < tol )
		return -1;				// y = 0, x is nonzero;

	/** 
	 * check whether x,y are linearly dependent over the rational number domain;
	 *		x, y are nonzero;
	 *	 dependent: ax=by, a,b are integers belonging to [1, int_reg);
	 * independent: cannot find two integers a,b belonging to [1, int_reg) 
	 *				to satisfy ax=by;
	 */
	double b, err;
	int intb;
	for ( int a = 1; a < int_reg; a++ )
	{
		b	 = a * x / y;
		intb = round(b);
		err  = fabs( b - round(b) );
		if ( err < tol )
		{
			if ( fabs(a) > fabs(intb) )
				return 1;		// x,y are nonzero and dependent, select x;
			else
				return 2;		// x,y are nonzero and dependent, select y;
		}
	}
	return 0;					// x,y are nonzero and independent;
}
