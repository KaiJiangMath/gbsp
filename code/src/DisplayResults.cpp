/*! \file	DisplayResults.cpp
 *
 *  \brief	Display results;
 *			calculate the corresponding drawing data for displaying results;
 *
 */

#include "Data.h"
#include "Head.h"
#include "DataOperators.h"
#include "umfpack.h"
#include "Mytimer.h"
#include "functs.h"


/**
 * \brief	display densities of bulk phase;
 */
void fn_disp_bulk_density	(	tvec<fftw_complex>		src, 
								stu_bulk_param			*sbulkparam, 
									int					opt_iter,
									int					step )
{
	mytimer_t timer;
	timer.reset();
	timer.start();

	FILE *densityfile;
	char densityName[FILELEN];
	sprintf(densityName, "%s/bulk%d_opt%d_density%d.dat", 
			rsltDir, sbulkparam->sflag, opt_iter, step);
	densityfile = fopen(densityName, "w");

	if ( sbulkparam->dimPhy == 2 )
	{
		double enlarge = sbulkparam->enlarge;
		int ny = (int) sbulkparam->plot_num, nz = ny; // the number of discrete points in real space;
		double sizey = enlarge*sbulkparam->dirBox.val[0][0];
		double sizez = enlarge*sbulkparam->dirBox.val[0][0];
		double dy = sizey / (double)(ny-1);
		double dz = sizez / (double)(nz-1);

		// save the number of discrete points;
		fprintf(densityfile, "%d\t%d\t%d\n", ny, nz, ny*nz);

		// save discrete grids along each direction;
		for (int ky = 0; ky < ny; ky++) 
			fprintf(densityfile, "%+.15E\n", dy * (ky-ny/2.0));
		for (int kz = 0; kz < nz; kz++)
			fprintf(densityfile, "%+.15E\t", dz * (kz-nz/2.0));

		// save density; projection;
		for (int ky = 0; ky < ny; ky++)
		{
			double yVal = dy * (ky - ny/2.0);
			for (int kz = 0; kz < nz; kz++)
			{
				double zVal = dz * (kz - nz/2.0);
				double rho = 0.0;
				double tmpphase;

				// rhoCplx * exp(i (k1*y+k2*z));
				for (int i = 0; i < sbulkparam->cplxDofs; i++)
				{
					double elm = fn_complex_abs ( src.val[i] );
					if ( elm > ZEROTOL )
					{
						tmpphase = sbulkparam->sfftv->projPlane.val[i][0]*yVal + 
								sbulkparam->sfftv->projPlane.val[i][1]*zVal;
						rho += (src.val[i][0]*cos(tmpphase) - src.val[i][1]*sin(tmpphase));
					}
				}
				fprintf(densityfile, "%+.15E\n", rho);
			}
		}
	}
	else if ( sbulkparam->dimPhy == 3 )
	{
		double enlarge = sbulkparam->enlarge;
		int nx = (int) sbulkparam->plot_num, ny = nx, nz = nx; // the number of discrete points in real space;
		double sizex = enlarge*sbulkparam->dirBox.val[0][0];
		double sizey = enlarge*sbulkparam->dirBox.val[0][0];
		double sizez = enlarge*sbulkparam->dirBox.val[0][0];
		double dx = sizex / (double)(nx-1);
		double dy = sizey / (double)(ny-1);
		double dz = sizez / (double)(nz-1);

		// save the number of discrete points;
		fprintf(densityfile, "%d\t%d\t%d\t%d\n", nx, ny, nz, nx*ny*nz);

		// save discrete grids along each direction;
		for (int kx = 0; kx < nx; kx++)
			fprintf(densityfile, "%+.15E\n", dx * (kx-nx/2.0));
		for (int ky = 0; ky < ny; ky++)
			fprintf(densityfile, "%+.15E\n", dy * (ky-ny/2.0));
		for (int kz = 0; kz < nz; kz++)
			fprintf(densityfile, "%+.15E\n", dz * (kz-nz/2.0));

		// save density; projection;
		for (int kx = 0; kx < nx; kx++)
		{
			double xVal = dx * (kx - nx/2.0);
			for (int ky = 0; ky < ny; ky++)
			{
				double yVal = dy * (ky - ny/2.0);
				for (int kz = 0; kz < nz; kz++)
				{
					double zVal = dz * (kz - nz/2.0);
					double rho = 0.0;
					double tmpphase;

					// rhoCplx * exp(i (k1*y+k2*z));
					for (int i = 0; i < sbulkparam->cplxDofs; i++)
					{
						double elm = fn_complex_abs ( src.val[i] );
						if ( elm > ZEROTOL )
						{
							tmpphase = sbulkparam->sfftv->projPlane.val[i][0]*xVal + 
								sbulkparam->sfftv->projPlane.val[i][1]*yVal + 
								sbulkparam->sfftv->projPlane.val[i][2]*zVal;
							rho += (src.val[i][0]*cos(tmpphase) - src.val[i][1]*sin(tmpphase));
						}
					}
					fprintf(densityfile, "%+.15E\n", rho);
				}
			}
		}
	}
	timer.pause();
	if ( sbulkparam->sflag == 1 )
		printf("\n===> Output plotted data (density of left bulk phase), ");
	else if ( sbulkparam->sflag == 2 )
		printf("\n===> Output plotted data (density of right bulk phase), ");
	printf("Step %d: %f seconds\n\n", step, timer.get_current_time());
	fclose(densityfile);
}


bool myComp(mySortVec a, mySortVec b)
{
	return (a.Data[0]*a.Data[0]+a.Data[1]*a.Data[1] > b.Data[0]*b.Data[0]+b.Data[1]*b.Data[1]);
}


bool scaleComp(mySortWave a, mySortWave b)
{
	return (a.Data > b.Data);
}


/**
 * \brief	display Fourier coefficients of bulk phase;
 */
void fn_disp_Fourier_coeff	(	tvec<fftw_complex>		src, 
								stu_bulk_param			*sbulkparam, 
									int					opt_iter,
									int					step )
{
	FILE *fFourCoeff;
	char fname[FILELEN];
	sprintf(fname, "%s/bulk%d_opt%d_field%d.dat", 
			rsltDir, sbulkparam->sflag, opt_iter, step);
	fFourCoeff = fopen(fname, "w");

	mySortVec *myVector = (mySortVec*)malloc(sizeof(mySortVec)*sbulkparam->cplxDofs);
	for (int k = 0; k < sbulkparam->cplxDofs; k++)
	{
		myVector[k].Data  = (double *)malloc(sizeof(double)*2);
		myVector[k].Index = (int *)malloc(sizeof(int)*sbulkparam->dimCpt);
		for (int i = 0; i < sbulkparam->dimCpt; i++)
			myVector[k].Index[i] = sbulkparam->sfftv->indKspace.val[k][i];
		for (int j = 0; j < 2; j++)
			myVector[k].Data[j] = src.val[k][j];
	}

	// save the Fourier coefficients whose values are greater than 1.0e-16;
	std::sort(myVector, myVector+sbulkparam->cplxDofs, myComp);
	for (int k = 0; k < sbulkparam->cplxDofs; k++)
	{
		double tmp = pow(myVector[k].Data[0], 2) + pow(myVector[k].Data[1], 2);
		tmp = sqrt(tmp);
		if (tmp > 1.0e-16)
		{
			for (int i = 0; i < sbulkparam->dimCpt; i++)
			{
				fprintf(fFourCoeff, "% d\t", myVector[k].Index[i]);
			}
			fprintf(fFourCoeff, "% e\t% e\n", myVector[k].Data[0], myVector[k].Data[1]);
		}
	}

	for (int i = 0; i < sbulkparam->cplxDofs; i++)
	{
		free(myVector[i].Index);
		free(myVector[i].Data);
	}
	free(myVector);
	fclose(fFourCoeff);
}


/**
 * \brief	display plane waves of bulk phase;
 */
void fn_disp_bulk_plane_wave	(	tvec<fftw_complex>		src, 
									stu_bulk_param			*sbulkparam, 
										int					opt_iter,
										int					step)
{
	FILE *fprojPlane;
	char projPlaneName[FILELEN];
	sprintf(projPlaneName, "%s/bulk%d_opt%d_planeWave%d.dat", 
			rsltDir, sbulkparam->sflag, opt_iter, step);
	fprojPlane = fopen(projPlaneName, "w");

	mySortWave *myPlaneWave = (mySortWave*)malloc(sizeof(mySortWave)*sbulkparam->cplxDofs);
	for (int k = 0; k < sbulkparam->cplxDofs; k++)
	{
		myPlaneWave[k].Wave = (double *)malloc(sizeof(double) * sbulkparam->dimPhy);
		for (int i = 0; i < sbulkparam->dimPhy; i++)
			myPlaneWave[k].Wave[i] = sbulkparam->sfftv->projPlane.val[k][i];
		myPlaneWave[k].Data = sqrt(src.val[k][0]*src.val[k][0]+src.val[k][1]*src.val[k][1]);
	}

	std::sort(myPlaneWave, myPlaneWave+sbulkparam->cplxDofs, scaleComp);

	// save the Fourier coefficients whose values are greater than 1.0e-16;
	for (int i = 0; i < sbulkparam->cplxDofs; i ++)
	{
		if (myPlaneWave[i].Data > 1.0e-16)
		{
			for (int j = 0; j < sbulkparam->dimPhy; j++)
				fprintf(fprojPlane, "% f\t ", myPlaneWave[i].Wave[j]);
			fprintf(fprojPlane, "% e\n", myPlaneWave[i].Data);
		}
	}
	
	for (int i = 0; i < sbulkparam->cplxDofs; i++)
	{
		free(myPlaneWave[i].Wave);
	}
	free(myPlaneWave);
	fclose(fprojPlane);
}


/*
 * \brief	display density of bulk phases after projecting to the space of GJPs;
 *
 * \param	rhoJCplx:		coefficients in Fourier and GJP space at the same time;
 * \param	sbulkparam:		structure body for bulk phase;
 * \param	ssysparam:		structure body for the interface system;
 *							including some plotting parameters and sGJPv;
 */
void fn_disp_bulk_proj_density (	tvec<fftw_complex>	rhoJCplx, 
									stu_bulk_param		*sbulkparam, 
									stu_system_param	*ssysparam )
{
	/* parameters; */
	int		nd		= sbulkparam->sbndv->nd;
	int		xlen	= sbulkparam->sbndv->xlen;
	int	cplxDofs	= sbulkparam->sbndv->cplxDofs;
	double	y_range = ssysparam->y_range;
	double	z_range = ssysparam->z_range;
	int		y_num	= ssysparam->y_num;
	int		z_num	= ssysparam->z_num;

	double dy = 2.0*PI*y_range / y_num;		// step along y-direction;
	double dz = 2.0*PI*z_range / z_num;		// step along z-direction;

	mytimer_t timer;
	timer.reset();
	timer.start();

	/* memory allocation for temporary variables; */
	int ind0, ind1;
	tvec<fftw_complex> tmp = fn_tvec_init<fftw_complex> ( xlen * cplxDofs );
	tmp.row = xlen;
	tmp.col = cplxDofs;
	fn_tvec_setZero_complex ( tmp );

	/** 
	 * project to Fourier space;
	 *		rhoJCplx * d0JJ + d0bndlr;
	 */
	for ( int i = 0; i < xlen; i++ )
	{
		for ( int j0 = 0; j0 < cplxDofs; j0++ )
		{
			ind1 = i*cplxDofs + j0;
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				ind0 = j1*cplxDofs + j0;
				tmp.val[ind1][0] += rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
				tmp.val[ind1][1] += rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
			}
			tmp.val[ind1][0] += sbulkparam->sbndv->d0bnd.val[ind1][0];
			tmp.val[ind1][1] += sbulkparam->sbndv->d0bnd.val[ind1][1];
		}
	}

	/* obtain global indexes of special positions (N/2); */
	tvec<int> gIndex = fn_obt_half_gIndex ( sbulkparam->NCpt, cplxDofs );

	/* open file; */
	FILE *densityfile;
	char densityName[FILELEN];
	sprintf(densityName, "%s/bulk%d_FGJP_density.dat", rsltDir, sbulkparam->sflag);
	densityfile = fopen(densityName, "w");

	/**
	 * calculate density by 
	 *		sum_{k} sum_{j} rhoJCplx d0JJ exp(i(tilde{R}'PBk)'tilde{r});
	 */
	if ( sbulkparam->dimPhy == 2 )
	{
		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\n", xlen, y_num, xlen*y_num);

		/* save discrete grids along each direction; */
		for ( int i = 0; i < xlen; i++ )
			fprintf(densityfile, "%+.15E\n", ssysparam->sGJPv->x.val[i]);
		for ( int i = 0; i < y_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", dy*i);

		/**
		 * save density rho;
		 *		rho:	size is xlen*y_num;
		 */
		double tmpPhase, rho;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j = 0; j < y_num; j++ )
			{
				rho = 0.0;
				int ind = 0;
				for ( int j0 = 0; j0 < cplxDofs; j0++ )
				{
					if ( j0 == gIndex.val[ind] ) // skip the special index;
						ind ++ ;
					else
					{
						ind1 = i*cplxDofs + j0;
						/* save elements whose modulus are greater than ZEROTOL; */
						double elm = fn_complex_abs ( tmp.val[ind1] );
						if ( elm > ZEROTOL )
						{
							tmpPhase = sbulkparam->sfftv->projPlane.val[j0][1]*dy*j;
							rho += tmp.val[ind1][0]*cos(tmpPhase) - tmp.val[ind1][1]*sin(tmpPhase);
						}
					}
				}
				fprintf(densityfile, "%+.15E\n", rho);
			}
		}
	}
	else if ( sbulkparam->dimPhy == 3 )
	{
		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\t%d\n", xlen, y_num, z_num, xlen*y_num*z_num);

		/* save discrete grids along each direction; */
		for ( int i = 0; i < xlen; i++ )
			fprintf(densityfile, "%+.15E\n", ssysparam->sGJPv->x.val[i]);
		for ( int i = 0; i < y_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", dy*i);
		for ( int i = 0; i < z_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", dz*i);

		/**
		 * save density rho;
		 *		rho:	size is xlen*y_num;
		 */
		double tmpPhase, rho;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j = 0; j < y_num; j++ )
			{
				for ( int k = 0; k < z_num; k++ )
				{
					rho = 0.0;
					int ind = 0;
					for ( int j0 = 0; j0 < cplxDofs; j0++ )
					{
						if ( j0 == gIndex.val[ind] ) // skip the special index;
							ind ++ ;
						else
						{
							ind1 = i*cplxDofs + j0;
							/* save elements whose modulus are greater than ZEROTOL; */
							double elm = fn_complex_abs ( tmp.val[ind1] );
							if ( elm > ZEROTOL )
							{
								tmpPhase = sbulkparam->sfftv->projPlane.val[j0][1]*dy*j +
										   sbulkparam->sfftv->projPlane.val[j0][2]*dz*k;
								rho += tmp.val[ind1][0]*cos(tmpPhase) - tmp.val[ind1][1]*sin(tmpPhase);
							}
						}
					}
					fprintf(densityfile, "%+.15E\n", rho);
				}
			}
		}
	}
	timer.pause();
	if ( sbulkparam->sflag == 1 )
		printf("\n===> Output plotted data (density of left bulk phase after projecting): ");
	else if ( sbulkparam->sflag == 2 )
		printf("\n===> Output plotted data (density of right bulk phase after projecting): ");
	printf("%f seconds\n\n", timer.get_current_time());
	fclose(densityfile);

	/* release memory; */
	fn_tvec_free<fftw_complex> ( tmp );
	fn_tvec_free<int> ( gIndex );
}


/*
 * \brief	display density of bulk phases after projecting to the space of GJPs;
 *
 * \param	rhoReJCplx:		coefficients in the common space;
 * \param	sbulkparam:		structure body for bulk phase;
 * \param	ssysparam:		structure body for the interface system;
 *							including some plotting parameters and sGJPv;
 */
void fn_disp_bulk_common_density (	tvec<fftw_complex>	rhoReJCplx, 
									stu_bulk_param		*sbulkparam, 
									stu_system_param	*ssysparam )
{
	/* parameters; */
	int		nd		= sbulkparam->srebndv->nd;
	int		xlen	= sbulkparam->srebndv->xlen;
	int cplxReDofs  = sbulkparam->srebndv->cplxDofs;
	double	y_range = ssysparam->y_range;
	double	z_range = ssysparam->z_range;
	int		y_num	= ssysparam->y_num;
	int		z_num	= ssysparam->z_num;

	double dy = 2.0*PI*y_range / y_num;		// step along y-direction;
	double dz = 2.0*PI*z_range / z_num;		// step along z-direction;

	mytimer_t timer;
	timer.reset();
	timer.start();

	/* memory allocation for temporary variables; */
	int ind0, ind1;
	tvec<fftw_complex> tmp = fn_tvec_init<fftw_complex> ( xlen * cplxReDofs );
	tmp.row = xlen;
	tmp.col = cplxReDofs;
	fn_tvec_setZero_complex ( tmp );

	/** 
	 * project to Fourier space;
	 *		rhoReJCplx * d0JJ + d0rebndlr;
	 */
	for ( int i = 0; i < xlen; i++ )
	{
		for ( int j0 = 0; j0 < cplxReDofs; j0++ )
		{
			ind1 = i*cplxReDofs + j0;
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				ind0 = j1*cplxReDofs + j0;
				tmp.val[ind1][0] += rhoReJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
				tmp.val[ind1][1] += rhoReJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
			}
			tmp.val[ind1][0] += sbulkparam->srebndv->d0bnd.val[ind1][0];
			tmp.val[ind1][1] += sbulkparam->srebndv->d0bnd.val[ind1][1];
		}
	}

	/* obtain global indexes of special positions (N/2); */
	tvec<int> gIndex = fn_obt_half_gIndex ( ssysparam->NCpt, cplxReDofs );

	/* open file; */
	FILE *densityfile;
	char densityName[FILELEN];
	sprintf(densityName, "%s/bulk%d_FGJPre_density.dat", rsltDir, sbulkparam->sflag);
	densityfile = fopen(densityName, "w");

	/**
	 * calculate density by 
	 *		sum_{k} sum_{j} rhoReJCplx d0JJ exp(i(common rotateProjBoxMat)'tilde{r});
	 */
	if ( sbulkparam->dimPhy == 2 )
	{
		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\n", xlen, y_num, xlen*y_num);

		/* save discrete grids along each direction; */
		for ( int i = 0; i < xlen; i++ )
			fprintf(densityfile, "%+.15E\n", ssysparam->sGJPv->x.val[i]);
		for ( int i = 0; i < y_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", dy*i);

		/**
		 * save density rho;
		 *		rho:	size is xlen*y_num;
		 */
		double tmpPhase, rho;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j = 0; j < y_num; j++ )
			{
				rho = 0.0;
				int ind = 0;
				for ( int j0 = 0; j0 < cplxReDofs; j0++ )
				{
					if ( j0 == gIndex.val[ind] ) // skip the special index;
						ind ++ ;
					else
					{
						ind1 = i*cplxReDofs + j0;
						/* save elements whose modulus are greater than ZEROTOL; */
						double elm = fn_complex_abs ( tmp.val[ind1] );
						if ( elm > ZEROTOL )
						{
							tmpPhase = ssysparam->sfftv->projPlane.val[j0][0]*dy*j;
							rho += tmp.val[ind1][0]*cos(tmpPhase) - tmp.val[ind1][1]*sin(tmpPhase);
						}
					}
				}
				fprintf(densityfile, "%+.15E\n", rho);
			}
		}
	}
	else if ( sbulkparam->dimPhy == 3 )
	{
		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\t%d\n", xlen, y_num, z_num, xlen*y_num*z_num);

		/* save discrete grids along each direction; */
		for ( int i = 0; i < xlen; i++ )
			fprintf(densityfile, "%+.15E\n", ssysparam->sGJPv->x.val[i]);
		for ( int i = 0; i < y_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", dy*i);
		for ( int i = 0; i < z_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", dz*i);
	
		/**
		 * save density rho;
		 *		rho:	size is xlen*y_num;
		 */
		double tmpPhase, rho;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j = 0; j < y_num; j++ )
			{
				for ( int k = 0; k < z_num; k++ )
				{
					rho = 0.0;
					int ind = 0;
					for ( int j0 = 0; j0 < cplxReDofs; j0++ )
					{
						if ( j0 == gIndex.val[ind] ) // skip the special index;
							ind ++ ;
						else
						{
							ind1 = i*cplxReDofs + j0;
							/* save elements whose modulus are greater than ZEROTOL; */
							double elm = fn_complex_abs ( tmp.val[ind1] );
							if ( elm > ZEROTOL )
							{
								tmpPhase = ssysparam->sfftv->projPlane.val[j0][0]*dy*j +
										   ssysparam->sfftv->projPlane.val[j0][1]*dz*k;
								rho += tmp.val[ind1][0]*cos(tmpPhase) - tmp.val[ind1][1]*sin(tmpPhase);
							}
						}
					}
					fprintf(densityfile, "%+.15E\n", rho);
				}
			}
		}
	}
	timer.pause();
	if ( sbulkparam->sflag == 1 )
		printf("\n===> Output plotted data (density of left bulk phase int the common space): ");
	else if ( sbulkparam->sflag == 2 )
		printf("\n===> Output plotted data (density of right bulk phase in the common space): ");
	printf("%f seconds\n\n", timer.get_current_time());
	fclose(densityfile);

	/* release memory; */
	fn_tvec_free<fftw_complex> ( tmp );
	fn_tvec_free<int> ( gIndex );
}


/*
 * \brief	display density of the interfce system;
 *
 * \param	rhoJCplx:		coefficients of the interface system;
 * \param	ssysparam:		structure body for the interface system;
 *							including some plotting parameters and sGJPv;
 */
void fn_disp_system_density		(	tvec<fftw_complex>		rhoJCplx, 
									stu_system_param		*ssysparam,
										int					step )
{
	/* parameters; */
	int		nd		= ssysparam->scbndv->nd;
	int		xlen	= ssysparam->scbndv->xlen;
	int cplxReDofs  = ssysparam->scbndv->cplxDofs;
	double	y_range = ssysparam->y_range;
	double	z_range = ssysparam->z_range;
	int		y_num	= ssysparam->y_num;
	int		z_num	= ssysparam->z_num;

	double dy = 2.0*PI*y_range / y_num;		// step along y-direction;
	double dz = 2.0*PI*z_range / z_num;		// step along z-direction;

	mytimer_t timer;
	timer.reset();
	timer.start();

	/* memory allocation for temporary variables; */
	int ind0, ind1;
	tvec<fftw_complex> tmp = fn_tvec_init<fftw_complex> ( xlen * cplxReDofs );
	tmp.row = xlen;
	tmp.col = cplxReDofs;
	fn_tvec_setZero_complex ( tmp );

	/** 
	 * project to Fourier space;
	 *		rhoJCplx * d0JJ + d0rebndlr;
	 */
	for ( int i = 0; i < xlen; i++ )
	{
		for ( int j0 = 0; j0 < cplxReDofs; j0++ )
		{
			ind1 = i*cplxReDofs + j0;
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				ind0 = j1*cplxReDofs + j0;
				tmp.val[ind1][0] += rhoJCplx.val[ind0][0] * ssysparam->sGJPv->d0JJ.val[i][j1];
				tmp.val[ind1][1] += rhoJCplx.val[ind0][1] * ssysparam->sGJPv->d0JJ.val[i][j1];
			}
			tmp.val[ind1][0] += ssysparam->scbndv->d0bnd.val[ind1][0];
			tmp.val[ind1][1] += ssysparam->scbndv->d0bnd.val[ind1][1];
		}
	}

	/* obtain global indexes of special positions (N/2); */
	tvec<int> gIndex = fn_obt_half_gIndex ( ssysparam->NCpt, cplxReDofs );
//	fn_tvec_print<int> ( gIndex );

	/* open file; */
	FILE *densityfile;
	char densityName[FILELEN];
	sprintf(densityName, "%s/sys_density%d.dat", rsltDir, step);
	densityfile = fopen(densityName, "w");

	/**
	 * calculate density by 
	 *		sum_{k} sum_{j} rhoJCplx d0JJ exp(i(common rotateProjBoxMat)'tilde{r});
	 */
	if ( ssysparam->dimRePhy == 1 )
	{
		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\n", xlen, y_num, xlen*y_num);

		/* save discrete grids along each direction; */
		for ( int i = 0; i < xlen; i++ )
			fprintf(densityfile, "%+.15E\n", ssysparam->sGJPv->x.val[i]);
		for ( int i = 0; i < y_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", dy*i);

		/**
		 * save density rho;
		 *		rho:	size is xlen*y_num;
		 */
		double tmpPhase, rho;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j = 0; j < y_num; j++ )
			{
				rho = 0.0;
				int ind = 0;
				for ( int j0 = 0; j0 < cplxReDofs; j0++ )
				{
					if ( j0 == gIndex.val[ind] ) // skip the special index;
						ind ++ ;
					else
					{
						ind1 = i*cplxReDofs + j0;
						/* save elements whose modulus are greater than ZEROTOL; */
						double elm = fn_complex_abs ( tmp.val[ind1] );
						if ( elm > ZEROTOL )
						{
							tmpPhase = ssysparam->sfftv->projPlane.val[j0][0]*dy*j;
							rho += tmp.val[ind1][0]*cos(tmpPhase) - tmp.val[ind1][1]*sin(tmpPhase);
						}
					}
				}
				fprintf(densityfile, "%+.15E\n", rho);
			}
		}
	}
	else if ( ssysparam->dimRePhy == 2 )
	{
		/* save the number of discrete points; */
		fprintf(densityfile, "%d\t%d\t%d\t%d\n", xlen, y_num, z_num, xlen*y_num*z_num);

		/* save discrete grids along each direction; */
		for ( int i = 0; i < xlen; i++ )
			fprintf(densityfile, "%+.15E\n", ssysparam->sGJPv->x.val[i]);
		for ( int i = 0; i < y_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", dy*i);
		for ( int i = 0; i < z_num; i++ ) 
			fprintf(densityfile, "%+.15E\n", dz*i);

		/**
		 * save density rho;
		 *		rho:	size is xlen*y_num;
		 */
		double tmpPhase, rho;
		for ( int i = 0; i < xlen; i++ )
		{
			for ( int j = 0; j < y_num; j++ )
			{
				for ( int k = 0; k < z_num; k++ )
				{
					rho = 0.0;
					int ind = 0;
					for ( int j0 = 0; j0 < cplxReDofs; j0++ )
					{
						if ( j0 == gIndex.val[ind] ) // skip the special index;
							ind ++ ;
						else
						{
							ind1 = i*cplxReDofs + j0;
							/* save elements whose modulus are greater than ZEROTOL; */
							double elm = fn_complex_abs ( tmp.val[ind1] );
							if ( elm > ZEROTOL )
							{
								tmpPhase = ssysparam->sfftv->projPlane.val[j0][0]*dy*j +
										   ssysparam->sfftv->projPlane.val[j0][1]*dz*k;
								rho += tmp.val[ind1][0]*cos(tmpPhase) - tmp.val[ind1][1]*sin(tmpPhase);
							}
						}
					}
					fprintf(densityfile, "%+.15E\n", rho);
				}
			}
		}
	}
	timer.pause();
	printf("\n===> Output plotted data (density of interface system): ");
	printf("%f seconds\n\n", timer.get_current_time());
	fclose(densityfile);

	/* release memory; */
	fn_tvec_free<fftw_complex> ( tmp );
	fn_tvec_free<int> ( gIndex );
}
