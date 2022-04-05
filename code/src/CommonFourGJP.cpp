/*! \file	CommonFourGJP.cpp
 *
 * \brief	Transform rhoJCplx (in the space of Fourier and GJPs) to 
 *			the common space (Fourier and GJPs);
 *
 * \note!	Please note 'break' in 'fn_re_represent_common';
 *			deleting 'break' will bring the same result as the MATLAB code;
 */

#include "Data.h"
#include "DataOperators.h"
#include "Head.h"
#include "Mytimer.h"
#include "functs.h"


/**
 * \brief	Initialization and preparation;
 *			Transform rhoJCplx of the two bulk phases in a common space;	
 *
 * \param	sbulkparam:		the structure body for bulk phases;
 * \param	ssysparam:		the structure body for the interface system;
 *
 */
void fn_common_Fourier_GJP (	stu_bulk_param		*sbulkparam1,
								stu_bulk_param		*sbulkparam2,
								stu_system_param	*ssysparam )
{
	printf(" <========== Transform to the common space ==========> \n\n");
	mytimer_t timer;
	timer.reset();
	timer.start();

	/* obtain the projection plane about the common 'rotateProjBoxMat'; */
	fn_get_system_projection_plane ( ssysparam );

//	sbulkparam1->srebndv->isSave  = true;
//	sbulkparam2->srebndv->isSave  = true;
//	sbulkparam1->srebndv->isTest  = true;
//	sbulkparam2->srebndv->isTest  = true;

	/**
	 * re-represent by the common 'rotateProjBoxMat';
	 *	including 'rhoJCplx', 'd0bnd', ...;
	 */
	fn_total_re_represent_common ( sbulkparam1, ssysparam );
	fn_total_re_represent_common ( sbulkparam2, ssysparam );

	/* release memory about 'sbndv'; */
	fn_bnd_memory_free ( sbulkparam1->sbndv, "total" );
	fn_bnd_memory_free ( sbulkparam2->sbndv, "total" );

	timer.pause();
	printf("\n\t ***** time cost of transformation to the common space: ");
	printf("%f seconds *****\n\n", timer.get_current_time());
}


/**
 * \brief	Get projection plane in the interface system; i.e. R'PBk;
 *			R is rotateMat;
 *			P is projMat;
 *			B is rcpBox;
 *			k is Fourier index;
 *			R'PB: [dimRePhy, dimReCpt];		k: [dimReCpt, 1];
 */
void fn_get_system_projection_plane ( stu_system_param		*ssysparam )
{
	/* parameters; */
	int dimRePhy	  = ssysparam->dimRePhy;
	int dimReCpt	  = ssysparam->dimReCpt;
	int Four_num	  = ssysparam->Four_num;
	int cplxReDofs    = 1;
	for ( int i = 0; i < dimReCpt; i++ )
		cplxReDofs *= Four_num;
	ssysparam->cplxReDofs = cplxReDofs;

	/* memory allocation; */
	ssysparam->NCpt = fn_tvec_init<int> ( dimReCpt );
	for ( int i = 0; i < dimReCpt; i++ )
		ssysparam->NCpt.val[i] = Four_num;
	ssysparam->sfftv->indKspace = fn_tmat_init<int>		( cplxReDofs, dimReCpt );
	ssysparam->sfftv->projPlane	= fn_tmat_init<double>	( cplxReDofs, dimRePhy );

	/* obtain k index; */
	fn_obt_kIndex ( ssysparam->sfftv->indKspace, ssysparam->NCpt );

	/* obtain projection wave; */
	double RPBKtmp;
	for (int i = 0; i < cplxReDofs; i++)
	{
		for (int kk = 0; kk < dimRePhy; kk++)
		{
			RPBKtmp = 0.0;
			for (int jj = 0; jj < dimReCpt; jj++) // R'PB * k;
				RPBKtmp += ssysparam->rotateProjBoxMat.val[kk][jj] * 
							ssysparam->sfftv->indKspace.val[i][jj];
			ssysparam->sfftv->projPlane.val[i][kk] = RPBKtmp;
		}
	}
}


/**
 * \brief	summarization of 'fn_re_represent_common';
 */
void fn_total_re_represent_common (stu_bulk_param		*sbulkparam,
								   stu_system_param		*ssysparam)
{
	/* memory allocation; */
	fn_common_memory_alloc (sbulkparam->sbndv, sbulkparam->srebndv, 
							ssysparam->cplxReDofs, "total");

	/* re-represent 'rhoJCplx' by the common 'rotateProjBoxMat'; */
	fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->rhoJCplx, 
							sbulkparam->srebndv->rhoJCplx, "rhoReJCplx");

	/* re-represent 'd0bnd', ..., by the common 'rotateProjBoxMat'; */
	fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d0bnd, 
							sbulkparam->srebndv->d0bnd, "d0rebnd");
	fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d1bnd, 
							sbulkparam->srebndv->d1bnd, "d1rebnd");
	fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d2bnd, 
							sbulkparam->srebndv->d2bnd, "d2rebnd");
	if ( strcmp(model_type, "LP") == 0 )
	{
		fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d3bnd, 
								sbulkparam->srebndv->d3bnd, "d3rebnd");
		fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d4bnd, 
								sbulkparam->srebndv->d4bnd, "d4rebnd");
		fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d5bnd, 
								sbulkparam->srebndv->d5bnd, "d5rebnd");
		fn_re_represent_common (sbulkparam, ssysparam, sbulkparam->sbndv->d6bnd, 
								sbulkparam->srebndv->d6bnd, "d6rebnd");
	}

	/* save density calculating by re-represented rhoJCplx; */
	if ( sbulkparam->srebndv->isSave )
	{
		fn_disp_bulk_common_density (sbulkparam->srebndv->rhoJCplx, sbulkparam, ssysparam);
		/* calculate error; */
		fn_proj_common_error ( sbulkparam, ssysparam );
	}
}


/**
 * \brief	memory allocation;
 *
 * \param	sbndv:		provide parameters;
 * \param	sobjbndv:	the objective structure body;
 * \param	cplxReDofs:	a parameter provided by 'ssysparam';
 */
void fn_common_memory_alloc (stu_bnd_var		*sbndv,
							 stu_bnd_var		*sobjbndv,
								int				cplxReDofs,
							const char			*memType)
{
	sobjbndv->polyDegree = sbndv->polyDegree;
	sobjbndv->nd		 = sbndv->nd;
	sobjbndv->xlen		 = sbndv->xlen;
	sobjbndv->cplxDofs	 = cplxReDofs;
	sobjbndv->x_range	 = sbndv->x_range;

	int rhoLen = sobjbndv->nd   * sobjbndv->cplxDofs;
	sobjbndv->rhoJCplx	= fn_tvec_init<fftw_complex> ( rhoLen );

	int bndLen = sobjbndv->xlen * sobjbndv->cplxDofs;
	sobjbndv->d0bnd		= fn_tvec_init<fftw_complex> (bndLen);
	sobjbndv->d2bnd		= fn_tvec_init<fftw_complex> (bndLen);
	if ( strcmp(model_type, "LP") == 0 )
	{
		sobjbndv->d4bnd		= fn_tvec_init<fftw_complex> (bndLen);
		sobjbndv->d6bnd		= fn_tvec_init<fftw_complex> (bndLen);
	}

	if ( strcmp(memType, "total") == 0 )
	{
		sobjbndv->d1bnd		= fn_tvec_init<fftw_complex> (bndLen);
		if ( strcmp(model_type, "LP") == 0 )
		{
			sobjbndv->d3bnd		= fn_tvec_init<fftw_complex> (bndLen);
			sobjbndv->d5bnd		= fn_tvec_init<fftw_complex> (bndLen);
		}
	}
}


/**
 * \brief	Re-represent 'rhoJCplx' in the common 'rotateProjBoxMat';
 *			Re-represent 'd0bnd', ... in the common 'rotateProjBoxMat';
 *
 * \param	sbulkparam:		provide projection plane about bulk phases;
 * \param	ssysparam:		provide projection plane about the interface system;
 * \param	origBulk:		the original bulk phase (projected in the system);
 * \param	reComBulk:		the bulk phase transformed by the common 'rotateProjBoxMat';
 * \param	saveStr:		a string for saving data;
 *
 * \note!	The essential spectra are consistent;
 */
void fn_re_represent_common	(	stu_bulk_param		*sbulkparam,
								stu_system_param	*ssysparam,
								tvec<fftw_complex>	origBulk,
								tvec<fftw_complex>	reComBulk,
									const char		*saveStr )
{
	/* load parameters; */
	double tol = 1e-10;
	int nd		   = sbulkparam->sbndv->nd; 
	int xlen	   = sbulkparam->sbndv->xlen;
	int cplxDofs   = sbulkparam->sbndv->cplxDofs; // the number of spectra in bulk phase;
	int dimRePhy   = ssysparam->dimRePhy;
	int dimReCpt   = ssysparam->dimReCpt;
	int cplxReDofs = ssysparam->cplxReDofs;		 // the number of spectra for re-representation;

	/* considering two cases: 'rhoJCplx'; 'd0bnd',...; */
	int nj;
	if ( strcmp(saveStr, "rhoReJCplx") == 0 ) nj = nd;
	else nj = xlen;

	/* initialization of results; */
	int lenRe = nj * cplxReDofs;
	for ( int i = 0; i < lenRe; i++ )
		fn_complex_setZero ( reComBulk.val[i] );

	/* re-representation of 'rhoJCplx'; */
	int ind0, ind1;
	double res, tmp0, tmp1;
	for ( int i = 0; i < cplxDofs; i++ )
	{
		for ( int j = 0; j < cplxReDofs; j++ )
		{
			res = 0.0;
			for ( int j0 = 0; j0 < dimRePhy; j0++ )	// calculate residual; norm 2;
			{
				tmp0 = sbulkparam->sfftv->projPlane.val[i][j0+1];
				tmp1 = ssysparam->sfftv->projPlane.val[j][j0];
				res += pow(tmp0-tmp1, 2);
			}
			res = sqrt(res);
			if ( res < tol )
			{
				for ( int j1 = 0; j1 < nj; j1++ )	// same spectra;
				{
					ind0 = j1 * cplxReDofs + j;
					ind1 = j1 * cplxDofs + i;
					reComBulk.val[ind0][0] += origBulk.val[ind1][0];
					reComBulk.val[ind0][1] += origBulk.val[ind1][1];
				}
				/*
				if ( strcmp(saveStr, "rhoReJCplx") == 0 && j == 20 )
				{
					printf("reComBulk[%d] [%d,%d] : %.4e,%.4e\n", sbulkparam->sflag,
						i, j, reComBulk.val[0][0], reComBulk.val[0][1]);
					printf("\t --> %.4e,%.4e,%.4e \t %.4e,%.4e\n",
						sbulkparam->sfftv->projPlane.val[i][0],
						sbulkparam->sfftv->projPlane.val[i][1],
						sbulkparam->sfftv->projPlane.val[i][2],
						ssysparam->sfftv->projPlane.val[j][0],
						ssysparam->sfftv->projPlane.val[j][1]);
					printf("\n");
				}
				*/
				break;
			}
		}
	}

	/* save results; */
	if ( sbulkparam->srebndv->isTest )
	{
		char fwFile[FILELEN];
		sprintf(fwFile, "%s/bulk%d_%s.dat", rsltDir, sbulkparam->sflag, saveStr);
		FILE *fwPath = fopen(fwFile, "w");
		fprintf(fwPath, "%d\t%d\t%d\n", nj, cplxReDofs, lenRe);
		for ( int j = 0; j < lenRe; j++ )
		{
			fprintf(fwPath, "%+.15E\t%+.15E\n", reComBulk.val[j][0], reComBulk.val[j][1]);
		}
		fclose(fwPath);
		/* save 'projPlane'; */
		sprintf(fwFile, "%s/bulk%d_projPlane.dat", rsltDir, sbulkparam->sflag);
		fwPath = fopen(fwFile, "w");
		fprintf(fwPath, "%d\t%d\t%d\n", cplxDofs, dimRePhy+1, cplxDofs*(dimRePhy+1));
		for ( int i = 0; i < cplxDofs; i++ )
		{
			for ( int j = 0; j < dimRePhy+1; j++ )
				fprintf(fwPath, "%+.15E\t", sbulkparam->sfftv->projPlane.val[i][j]);
			fprintf(fwPath, "\n");
		}
		fclose(fwPath);
		sprintf(fwFile, "%s/sys_projPlane.dat", rsltDir);
		fwPath = fopen(fwFile, "w");
		fprintf(fwPath, "%d\t%d\t%d\n", cplxReDofs, dimRePhy, cplxReDofs*(dimRePhy+1));
		for ( int i = 0; i < cplxReDofs; i++ )
		{
			for ( int j = 0; j < dimRePhy; j++ )
				fprintf(fwPath, "%+.15E\t", ssysparam->sfftv->projPlane.val[i][j]);
			fprintf(fwPath, "\n");
		}
		fclose(fwPath);
	}
}


/**
 * \brief	Calculate the error between 'bulk1_FGJP_density' and 'bulk1_FGJPre_density';
 *			also for 'bulk2_FGJP_density' and 'bulk2_FGJPre_density';
 */
void fn_proj_common_error ( stu_bulk_param		*sbulkparam,
							stu_system_param	*ssysparam )
{
	char fname0[FILELEN];
	char fname1[FILELEN];
	sprintf(fname0, "%s/bulk%d_FGJP_density.dat", rsltDir, sbulkparam->sflag);
	sprintf(fname1, "%s/bulk%d_FGJPre_density.dat", rsltDir, sbulkparam->sflag);

	int		xlen0	= sbulkparam->sbndv->xlen;
	int		xlen1	= sbulkparam->srebndv->xlen;
	int		y_num	= ssysparam->y_num;
	int		z_num	= ssysparam->z_num;

	int status = 1;
	if ( sbulkparam->dimPhy == 2 )
	{
		/* read data from 'fname0'; */
		tvec<double> x0		= fn_tvec_init<double> ( xlen0 );
		tvec<double> y0		= fn_tvec_init<double> ( y_num );
		tvec<double> data0	= fn_tvec_init<double> ( xlen0 * y_num );
		status = fn_read_data_2d ( fname0, x0, y0, data0 );

 		/* read data from 'fname1'; */
		tvec<double> x1		= fn_tvec_init<double> ( xlen1 );
		tvec<double> y1		= fn_tvec_init<double> ( y_num );
		tvec<double> data1	= fn_tvec_init<double> ( xlen1 * y_num );
		status = fn_read_data_2d ( fname1, x1, y1, data1 );

		/* calculate error; */
		fn_tvec_add<double> ( x0,	 x1,	1.0, -1.0 ); // x0 - x1;
		fn_tvec_add<double> ( y0,	 y1,	1.0, -1.0 ); // y0 - y1;
		fn_tvec_add<double> ( data0, data1, 1.0, -1.0 ); // data0 - data1;
		double err_x	= fn_tvec_maxAbs<double> ( x0 );
		double err_y	= fn_tvec_maxAbs<double> ( y0 );
		double err_data = fn_tvec_maxAbs<double> ( data0 );
		printf("\t Error between the 'phi' before and after changing space.\n");
		printf("\t ---> x error = %.6e\n", err_x);
		printf("\t ---> y error = %.6e\n", err_y);
		printf("\t ---> density error = %.6e\n", err_data);

		/* release memory; */
		fn_tvec_free<double> ( x0	 );
		fn_tvec_free<double> ( x1	 );
		fn_tvec_free<double> ( y0	 );
		fn_tvec_free<double> ( y1	 );
		fn_tvec_free<double> ( data0 );
		fn_tvec_free<double> ( data1 );
	} 
	else if ( sbulkparam->dimPhy == 3 )
	{
		/* read data from 'fname0'; */
		tvec<double> x0		= fn_tvec_init<double> ( xlen0 );
		tvec<double> y0		= fn_tvec_init<double> ( y_num );
		tvec<double> z0		= fn_tvec_init<double> ( z_num );
		tvec<double> data0	= fn_tvec_init<double> ( xlen0 * y_num * z_num );
		fn_read_data_3d ( fname0, x0, y0, z0, data0 );

 		/* read data from 'fname1'; */
		tvec<double> x1		= fn_tvec_init<double> ( xlen1 );
		tvec<double> y1		= fn_tvec_init<double> ( y_num );
		tvec<double> z1		= fn_tvec_init<double> ( z_num );
		tvec<double> data1	= fn_tvec_init<double> ( xlen1 * y_num * z_num );
		fn_read_data_3d ( fname1, x1, y1, z1, data1 );

		/* calculate error; */
		fn_tvec_add<double> ( x0,	 x1,	1.0, -1.0 ); // x0 - x1;
		fn_tvec_add<double> ( y0,	 y1,	1.0, -1.0 ); // y0 - y1;
		fn_tvec_add<double> ( z0,	 z1,	1.0, -1.0 ); // z0 - z1;
		fn_tvec_add<double> ( data0, data1, 1.0, -1.0 ); // data0 - data1;
		double err_x	= fn_tvec_maxAbs<double> ( x0 );
		double err_y	= fn_tvec_maxAbs<double> ( y0 );
		double err_z	= fn_tvec_maxAbs<double> ( z0 );
		double err_data = fn_tvec_maxAbs<double> ( data0 );
		printf("\t Error between the 'phi' before and after changing space.\n");
		printf("\t ---> x error = %.6e\n", err_x);
		printf("\t ---> y error = %.6e\n", err_y);
		printf("\t ---> z error = %.6e\n", err_z);
		printf("\t ---> density error = %.6e\n", err_data);

		/* release memory; */
		fn_tvec_free<double> ( x0	 );
		fn_tvec_free<double> ( x1	 );
		fn_tvec_free<double> ( y0	 );
		fn_tvec_free<double> ( y1	 );
		fn_tvec_free<double> ( z0	 );
		fn_tvec_free<double> ( z1	 );
		fn_tvec_free<double> ( data0 );
		fn_tvec_free<double> ( data1 );
	}
}


/**
 * \brief	Read file for 2D data;
 */
int fn_read_data_2d		(		char		fname[FILELEN],		
							tvec<double>	x,
							tvec<double>	y,
							tvec<double>	data )
{
	FILE	*fp = fopen(fname, "r");
	int		status;
    if ( fp == NULL )	// invalid file;
	{
		printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 0;
	}
    
	int nx, ny, len;
    fscanf(fp, "%d", &(nx));
	fscanf(fp, "%d", &(ny));
	fscanf(fp, "%d", &(len));

	for ( int i = 0; i < nx; i++ )
		fscanf(fp, "%lf", &(x.val[i]));
	for ( int i = 0; i < ny; i++ )
		fscanf(fp, "%lf", &(y.val[i]));
	for ( int i = 0; i < len; i++ )
		fscanf(fp, "%lf", &(data.val[i]));
	fclose(fp);
	return 1;
}


/**
 * \brief	Read file for 3D data;
 */
int fn_read_data_3d		(		char		fname[FILELEN],		
							tvec<double>	x,
							tvec<double>	y,
							tvec<double>	z,
							tvec<double>	data )
{
	FILE	*fp = fopen(fname, "r");
	int		status;
    if ( fp == NULL )	// invalid file;
	{
		printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 0;
	}
    
	int nx, ny, nz, len;
    
    fscanf(fp, "%d", &(nx));
	fscanf(fp, "%d", &(ny));
	fscanf(fp, "%d", &(nz));
	fscanf(fp, "%d", &(len));

	for ( int i = 0; i < nx; i++ )
		fscanf(fp, "%lf", &(x.val[i]));
	for ( int i = 0; i < ny; i++ )
		fscanf(fp, "%lf", &(y.val[i]));
	for ( int i = 0; i < nz; i++ )
		fscanf(fp, "%lf", &(z.val[i]));
	for ( int i = 0; i < len; i++ )
		fscanf(fp, "%lf", &(data.val[i]));
	fclose(fp);
	return 1;
}
