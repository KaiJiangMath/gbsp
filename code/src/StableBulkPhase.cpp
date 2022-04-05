/*! \file	StableBulkPhase.cpp
 *
 * \brief	obtain the stable state of the bulk phase;
 *
 */

#include "Head.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"


/**
 * \brief	Initialization and preparation;
 *			update bulk phase to obtain the stable state;
 */
double fn_obt_stable_bulk_phase (	stu_bulk_param	*sbulkparam )
{
	printf(" <========== Fixed box, generate equilibrium state of ");
	if ( sbulkparam->sflag == 1 )
		printf("left bulk phase ==========> \n\n");
	else if ( sbulkparam->sflag == 2 )
		printf("right bulk phase ==========> \n\n");
	mytimer_t timer;
	timer.reset();
	timer.start();

	/* initialization and preparation; */
	fn_bulk_memory_allocation	( sbulkparam );			// memory allocation;
	fn_obt_kIndex ( sbulkparam->sfftv->indKspace, sbulkparam->NCpt ); // obtain Fourier k;
	fn_get_projection_plane	( sbulkparam );				// project plane; R'PBk;
	fn_get_Gsquare			( sbulkparam );				// Gsquare; |R'PBk|^2;
	fn_get_init_bulk_phase	( sbulkparam );				// obtain the initial rho in Fourier space;
	fn_save_bulk_phase		( sbulkparam, 0, 0 );		// save initial value;

	/* optimize computational box; */
	int opt_iter_max = sbulkparam->opt_iter_max;
	double opt_tol	 = sbulkparam->opt_tol;
	int opt_iter = 0;
	double opt_err = 1.0;
	double hamilton, hamiltonOld;

	/* update bulk phase to obtain the stable state; */
	hamilton = fn_update_bulk_phase	( sbulkparam, opt_iter );
	hamiltonOld = hamilton;

	while ( opt_err > opt_tol && opt_iter < opt_iter_max )
	{
		fn_opt_rcpBox ( sbulkparam );	// optimize computational box;
		hamilton = fn_update_bulk_phase	( sbulkparam, opt_iter ); // calculate stable state;

		/* calculate error; */
		opt_err = fabs( hamilton - hamiltonOld );
		printf("\t ***** Iter %d: hamilton = % .20e\terr = % .10e.\n\n", 
				opt_iter, hamilton, opt_err);
		
		/* update energy; */
		hamiltonOld = hamilton;
		opt_iter ++ ;
	}

	timer.pause();
	printf("\n\t ***** time cost of the calculation of stable ");
	if ( sbulkparam->sflag == 1 )
		printf("left bulk phase: %f seconds *****\n\n", timer.get_current_time());
	else if ( sbulkparam->sflag == 2 )
		printf("right bulk phase: %f seconds *****\n\n", timer.get_current_time());

	/* update 'dirBox' and 'translProjBoxVec'; */
	obtDualBox(sbulkparam->rcpBox.val, sbulkparam->dirBox.val, sbulkparam->dimCpt);
	fn_obt_translProjBox_matrix ( sbulkparam );
	printf("rcpBox: \n");
	fn_tmat_print<double> ( sbulkparam->rcpBox );
	printf("\n");

	/* release memory; only save rhoCplx; */
	fn_bulk_memory_free	(sbulkparam, "stable bulk");

	return hamilton;
}


/**
 * \brief	Allocation memory for calculating stable bulk phase;
 */
void fn_bulk_memory_allocation (stu_bulk_param *sbulkparam)
{
	/* load parameters; */
	int cplxDofs = sbulkparam->cplxDofs;
	int dimPhy	 = sbulkparam->dimPhy;
	int dimCpt	 = sbulkparam->dimCpt;

	sbulkparam->sfftv->indKspace = fn_tmat_init<int>			( cplxDofs, dimCpt );
	sbulkparam->sfftv->projPlane = fn_tmat_init<double>			( cplxDofs, dimPhy );
	sbulkparam->sfftv->Gsquare   = fn_tvec_init<double>			( cplxDofs );

	sbulkparam->sfftv->rhoCplx   = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->rhoReal	 = fn_tvec_init<fftw_complex>	( cplxDofs );

	sbulkparam->sfftv->fftw_Ctmp = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->fftw_Rtmp = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->gradient  = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->cplxTmp   = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->quadTerm  = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->cubTerm   = fn_tvec_init<fftw_complex>	( cplxDofs );
	sbulkparam->sfftv->quarTerm  = fn_tvec_init<fftw_complex>	( cplxDofs );

	/* initialization; */
	fn_tmat_setZero<int>	( sbulkparam->sfftv->indKspace );
	fn_tmat_setZero<double> ( sbulkparam->sfftv->projPlane );
	fn_tvec_setZero<double> ( sbulkparam->sfftv->Gsquare  );
	fn_tvec_setZero_complex ( sbulkparam->sfftv->rhoCplx  );
	fn_tvec_setZero_complex ( sbulkparam->sfftv->rhoReal  );
	fn_tvec_setZero_complex ( sbulkparam->sfftv->gradient );
	fn_tvec_setZero_complex ( sbulkparam->sfftv->cplxTmp  );

	/* prepare FFTW plan; */
	sbulkparam->sfftv->Planc2cFord = fftw_plan_dft ( dimCpt, 
			sbulkparam->NCpt.val, sbulkparam->sfftv->fftw_Rtmp.val, 
			sbulkparam->sfftv->fftw_Ctmp.val, FFTW_FORWARD,  FFTW_MEASURE );  // real to cplx
	sbulkparam->sfftv->Planc2cBack = fftw_plan_dft ( dimCpt, 
			sbulkparam->NCpt.val, sbulkparam->sfftv->fftw_Ctmp.val, 
			sbulkparam->sfftv->fftw_Rtmp.val, FFTW_BACKWARD, FFTW_MEASURE );  // cplx to real 
}


/**
 * \brief	Get projection plane; i.e. R'PBk;
 *			R is rotateMat;
 *			P is projMat;
 *			B is rcpBox;
 *			k is Fourier index;
 *			R'PB: [dimPhy, dimCpt];		k: [dimCpt, 1];
 */
void fn_get_projection_plane(stu_bulk_param *sbulkparam)
{
	double RPBKtmp;
	for (int i = 0; i < sbulkparam->cplxDofs; i++)
	{
		for (int kk = 0; kk < sbulkparam->dimPhy; kk++)
		{
			RPBKtmp = 0.0;
			for (int jj = 0; jj < sbulkparam->dimCpt; jj++)
				RPBKtmp += sbulkparam->rotateProjBoxMat.val[kk][jj] * 
							sbulkparam->sfftv->indKspace.val[i][jj]; // R'PB * k;
			sbulkparam->sfftv->projPlane.val[i][kk] = RPBKtmp;
		}
	}
}


/**
 * \brief	Get Gsquare; i.e. |R'PBk|^2;
 *			based on 'fn_get_projection_plane';
 */
void fn_get_Gsquare(stu_bulk_param		*sbulkparam)
{
	for (int i = 0; i < sbulkparam->cplxDofs; i++)
	{
		sbulkparam->sfftv->Gsquare.val[i] = 0.0;
		for (int kk = 0; kk < sbulkparam->dimPhy; kk++)
			sbulkparam->sfftv->Gsquare.val[i] += pow(sbulkparam->sfftv->projPlane.val[i][kk], 2);
	}
}


/**
 * \brief	Get the initial state of bulk phase;
 */
void fn_get_init_bulk_phase	(	stu_bulk_param	*sbulkparam )
{
	int dimCpt = sbulkparam->dimCpt;
	int kk;
	double translVal;
	tvec<int> kIndex = fn_tvec_init<int> ( dimCpt );
	for (int i = 0; i < sbulkparam->initNum; i++)
	{
		/* k to global index; */
		kk = fn_kIndex_to_gIndex ( sbulkparam->initIndex.val[i], sbulkparam->NCpt );
		/* consider translation; calculate (PBk)'t = k'B'P't; superscript ' means transpose; */
		/* t'PB * k; */
		translVal = 0.0;
		for (int j = 0; j < dimCpt; j++)
		{
			translVal += sbulkparam->translProjBoxVec.val[j] * 
						sbulkparam->initIndex.val[i][j]; // t'PB * k;
		}
		/* rhoCplx * exp(i(PBk)'t); t is translation matrix; */
		double tmpReal = sbulkparam->initCoeff.val[i][0] * cos(translVal) -
			sbulkparam->initCoeff.val[i][1] * sin(translVal);
		double tmpCplx = sbulkparam->initCoeff.val[i][0] * sin(translVal) + 
			sbulkparam->initCoeff.val[i][1] * cos(translVal);
		sbulkparam->sfftv->rhoCplx.val[kk][0] = tmpReal;
		sbulkparam->sfftv->rhoCplx.val[kk][1] = tmpCplx;

		/* symmetry; */
		if ( sbulkparam->initIndex.val[i][dimCpt-1] > 0 )
		{
			for ( int j = 0; j < dimCpt; j++ )
				kIndex.val[j] = -1 * sbulkparam->initIndex.val[i][j];
			/* k to global index; */
			kk = fn_kIndex_to_gIndex ( kIndex.val, sbulkparam->NCpt );
			/* consider translation; calculate (PBk)'t = k'B'P't; superscript ' means transpose; */
			/* t'PB * k; */
			translVal = 0.0;
			for (int j = 0; j < dimCpt; j++)
			{
				translVal += sbulkparam->translProjBoxVec.val[j] * 
							sbulkparam->initIndex.val[i][j]; // t'PB * k;
			}
			/* rhoCplx * exp(i(PBk)'t); t is translation matrix; */
			sbulkparam->sfftv->rhoCplx.val[kk][0] = tmpReal;
			sbulkparam->sfftv->rhoCplx.val[kk][1] = tmpCplx;
		}
	}
	fn_complex_setZero ( sbulkparam->sfftv->rhoCplx.val[0] );	// mass conservation;
	fn_tvec_free<int>  ( kIndex );
}


/**
 * \brief	Update bulk phase and obtain the stable bulk phase;
 */
double fn_update_bulk_phase (	stu_bulk_param		*sbulkparam,
									int				opt_iter )
{
	/* load iteration parameters; */
	double tol		= sbulkparam->tol;
	double stepSize = sbulkparam->step_size;
	int itMax		= sbulkparam->iter_max;
	int stepSave	= sbulkparam->save_step;

	/* initialization; */
	int iterator = 0;
	double res = 1.0;
	double tmp, temp;
	clock_t start, finish;
	double duration, hamilton, oldHamilton, diffham;
	char dataName[300];
	sprintf(dataName, "%s/bulk%d_energy_error.dat", rsltDir, sbulkparam->sflag);
	FILE *fdata = fopen(dataName, "w");

	oldHamilton = 100.0;
	diffham		= 0.0;
	tvec<double> energy = fn_tvec_init<double> ( 3 );
	fn_tvec_setZero<double> ( energy );

	while (res > tol && iterator < itMax)
	{
		/* calculate hamilton energy; */
		fn_calc_bulk_energy ( sbulkparam, energy );
		hamilton = energy.val[2];

		/* save data; */
		fprintf(fdata, "%+.15E\t%+.15E\t%+.20E\t%+.15E\n", 
				energy.val[0], energy.val[1], energy.val[2], res);
		if ( iterator%sbulkparam->print_step == 0 )
		{
			printf("\nIter %d: stepSize = %.5e \t res = %.10e \t diffham = %.10e\n", 
					iterator, stepSize, res, diffham);
			printf("\t Laplace term: % .10e \t Nonlinear term: % .10e,", energy.val[0], energy.val[1]);
			printf("\t hamilton = % .10e\n", energy.val[2]);
		}

		if ( (iterator > 0 ) && (iterator % stepSave) == 0 )
			fn_save_bulk_phase (sbulkparam, opt_iter, iterator);	// save data;

		diffham = hamilton - oldHamilton;
		oldHamilton = hamilton;

		for (int i = 0; i < sbulkparam->cplxDofs; i++)
		{
			/* [(q0^2-|G|^2)*(q1^2-|G|^2)*...]^2; */
			temp = 1.0;
			for (int j = 0; j < sbulkparam->scale_num; j++)
			{
				tmp = pow(sbulkparam->scale_val[j], 2) - sbulkparam->sfftv->Gsquare.val[i];
				temp *= pow(tmp, 2);
			}

			/* calculate gradient part; */
			sbulkparam->sfftv->gradient.val[i][0] = 
				-sbulkparam->model_tau	* sbulkparam->sfftv->rhoCplx.val[i][0] - 
				sbulkparam->model_gamma	* sbulkparam->sfftv->quadTerm.val[i][0] +
				sbulkparam->model_kappa * sbulkparam->sfftv->cubTerm.val[i][0];
			sbulkparam->sfftv->gradient.val[i][1] = 
				-sbulkparam->model_tau	* sbulkparam->sfftv->rhoCplx.val[i][1] - 
				sbulkparam->model_gamma	* sbulkparam->sfftv->quadTerm.val[i][1] + 
				sbulkparam->model_kappa	* sbulkparam->sfftv->cubTerm.val[i][1];

			/* semi-implicit scheme to update rhoCplx; */
			sbulkparam->sfftv->rhoCplx.val[i][0] -= stepSize * sbulkparam->sfftv->gradient.val[i][0];
			sbulkparam->sfftv->rhoCplx.val[i][1] -= stepSize * sbulkparam->sfftv->gradient.val[i][1];
			sbulkparam->sfftv->rhoCplx.val[i][0] /= (1.0 + stepSize*temp*sbulkparam->model_xi);
			sbulkparam->sfftv->rhoCplx.val[i][1] /= (1.0 + stepSize*temp*sbulkparam->model_xi);

//			printf("%d : temp = %.4e \t gradient = %.4e, %.4e \t rhoCplx = %.4e, %.4e.\n", i, temp,
//					sbulkparam->sfftv->gradient.val[i][0], sbulkparam->sfftv->gradient.val[i][1],
//					sbulkparam->sfftv->rhoCplx.val[i][0],  sbulkparam->sfftv->rhoCplx.val[i][1]);

			/* update gradient error; */
			sbulkparam->sfftv->gradient.val[i][0] += 
				sbulkparam->model_xi * temp * sbulkparam->sfftv->rhoCplx.val[i][0];
			sbulkparam->sfftv->gradient.val[i][1] += 
				sbulkparam->model_xi * temp * sbulkparam->sfftv->rhoCplx.val[i][1];
		}
		fn_complex_setZero ( sbulkparam->sfftv->rhoCplx.val[0]  );
		fn_complex_setZero ( sbulkparam->sfftv->gradient.val[0] );
		res = fn_tvec_maxAbs_complex ( sbulkparam->sfftv->gradient );

		iterator ++;
	}

	if ( itMax > 0 )
		fn_save_bulk_phase (sbulkparam, opt_iter, -1);		// save stationary value;

	// output data about the stationary state;
	fn_calc_bulk_energy ( sbulkparam, energy );
	hamilton = energy.val[2];
	fprintf(fdata, "%+.15E\t%+.15E\t%+.20E\t%+.15E\n", 
			energy.val[0], energy.val[1], energy.val[2], res);
	printf("\n\t ***** the bulk phase is %s\n", sbulkparam->phase);
	printf("\t ***** iterator = %d \t step size = %.3e\n", iterator, stepSize);
	printf("\t ***** res = % .10e \t diffham = % .10e\n", res, diffham);
	printf("\t ***** Laplace term: % .15e,\t Nonlinear term: % .15e\n", 
			energy.val[0], energy.val[1]);
	printf("\t ***** hamilton = % .20e\n\n", energy.val[2]);

	fclose(fdata);
	fn_tvec_free<double> ( energy );
	return hamilton;
}


/**
 * \brief	Calculate the free energy of bulk phase;
 */
void fn_calc_bulk_energy (	stu_bulk_param		*sbulkparam, 
							tvec<double>		energy )
{
	int cplxDofs = sbulkparam->cplxDofs;

	/* initialization; */
	fn_tvec_setZero<double> ( energy );
	tvec<fftw_complex> DiffTerm		= fn_tvec_init<fftw_complex> ( cplxDofs );
	tvec<fftw_complex> DiffTermTmp	= fn_tvec_init<fftw_complex> ( cplxDofs );

	/* (q0^2-|G|^2) * (q1^2-|G|^2) * ... * rhoCplx; */
	double Difftmp;
	for (int i = 0; i < sbulkparam->cplxDofs; i++)
	{
		Difftmp = 1.0;
		for (int j=0; j < sbulkparam->scale_num; j++)
		{
			Difftmp *= pow(sbulkparam->scale_val[j], 2) - sbulkparam->sfftv->Gsquare.val[i];  
		}
		DiffTermTmp.val[i][0] = Difftmp*sbulkparam->sfftv->rhoCplx.val[i][0];
		DiffTermTmp.val[i][1] = Difftmp*sbulkparam->sfftv->rhoCplx.val[i][1];
	}

	/* convolution calculation; */
	fn_convolution(DiffTermTmp, DiffTerm, sbulkparam->sfftv, 2); // [(nabla^2+q0^2)*(nabla^2+q1^2)*...*rho]^2;
	fn_convolution(sbulkparam->sfftv->rhoCplx, sbulkparam->sfftv->quadTerm, 
					sbulkparam->sfftv, 2); // rho^2;
	fn_convolution(sbulkparam->sfftv->rhoCplx, sbulkparam->sfftv->cubTerm, 
					sbulkparam->sfftv, 3); // rho^3;
	fn_convolution(sbulkparam->sfftv->rhoCplx, sbulkparam->sfftv->quarTerm, 
					sbulkparam->sfftv, 4); // rho^4;

	/* interaction energy; operator part; */
	energy.val[0] = DiffTerm.val[0][0]*sbulkparam->model_xi/2.0;
	/* entropy energy; nonlinear part; */
	energy.val[1] = 
		-sbulkparam->sfftv->quadTerm.val[0][0]	* sbulkparam->model_tau/2.0 - 
		sbulkparam->sfftv->cubTerm.val[0][0]	* sbulkparam->model_gamma/3.0 + 
		sbulkparam->sfftv->quarTerm.val[0][0]	* sbulkparam->model_kappa/4.0; 
	/* hamilton energy; */
	energy.val[2] = energy.val[0] + energy.val[1];

	fn_tvec_free<fftw_complex> ( DiffTerm );
	fn_tvec_free<fftw_complex> ( DiffTermTmp );
}


/**
 * \brief	The integration of functions about saving data;
 *			save Fourier coefficients, densities, plane wave;
 */
void fn_save_bulk_phase (	stu_bulk_param		*sbulkparam,
								int				opt_iter,
								int				iterator )
{	
	if ( sbulkparam->save_type[0] == 'y' || sbulkparam->save_type[0] == 'Y' )
		fn_disp_Fourier_coeff  ( sbulkparam->sfftv->rhoCplx, sbulkparam, opt_iter, iterator );
	if ( sbulkparam->save_type[1] == 'y' || sbulkparam->save_type[1] == 'Y' )
		fn_disp_bulk_density   ( sbulkparam->sfftv->rhoCplx, sbulkparam, opt_iter, iterator );
	if ( sbulkparam->save_type[2] == 'y' || sbulkparam->save_type[2] == 'Y' )
		fn_disp_bulk_plane_wave( sbulkparam->sfftv->rhoCplx, sbulkparam, opt_iter, iterator );
}


/**
 * \brief	Optimize the computational box;
 */
void fn_opt_rcpBox		( stu_bulk_param		*sbulkparam )
{
	int bbType	 = sbulkparam->box_bbType;
	int	iter_max = sbulkparam->box_iter_max;
	double step0 = sbulkparam->box_step_size;
	double  tol	 = sbulkparam->box_tol;

	/* memory allocation; */
	int row = sbulkparam->rcpBox.row;
	int col = sbulkparam->rcpBox.col;
	tmat<double> box0		= fn_tmat_init<double> ( row, col );
	tmat<double> box1		= fn_tmat_init<double> ( row, col );
	tmat<double> box		= fn_tmat_init<double> ( row, col );
	tmat<double> boxDiff	= fn_tmat_init<double> ( row, col );
	tmat<double> gradBOld	= fn_tmat_init<double> ( row, col );
	tmat<double> gradB		= fn_tmat_init<double> ( row, col );
	tmat<double> gradBDiff	= fn_tmat_init<double> ( row, col );
	fn_tmat_copy<double> ( box0,	sbulkparam->rcpBox );
	fn_tmat_copy<double> ( box,		sbulkparam->rcpBox );

	/* the step length for computing difference quotient; */
	double dh = 1.0e-4;

	/* initialization; */
	int		iter = 0;
	double  err	 = 1.0;
	double	step_size = step0;
	double lham, rham;
	tvec<double> energy = fn_tvec_init<double> ( 3 );
	fn_tvec_setZero<double> ( energy );

	/* recurrence for updating the computational box; */
	while ( err > tol && iter < iter_max )
	{
		/* compute the gradient box; */
		fn_tmat_setZero<double> ( gradB );
		if ( strcmp ( sbulkparam->boxType, "cube" ) == 0 )
		{
			/* calculate the energy of boxOld + dh; */
			fn_tmat_copy<double> ( sbulkparam->rcpBox, box );
			for ( int j = 0; j < col; j++ )
				sbulkparam->rcpBox.val[j][j] += dh;
			fn_obt_rotateProjBox_matrix ( sbulkparam );
			fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
			fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
			fn_calc_bulk_energy ( sbulkparam, energy );
			rham = energy.val[2];

			/* calculate the energy of boxOld - dh; */
			for ( int j = 0; j < col; j++ )
				sbulkparam->rcpBox.val[j][j] -= 2.0 * dh;
			fn_obt_rotateProjBox_matrix ( sbulkparam );
			fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
			fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
			fn_calc_bulk_energy ( sbulkparam, energy );
			lham = energy.val[2];

			/* gradient box; */
			for ( int j = 0; j < col; j++ )
				gradB.val[j][j] = ( rham - lham ) / ( 2.0*dh );
		}
		else if ( strcmp ( sbulkparam->boxType, "cuboid" ) == 0 )
		{
			for ( int j = 0; j < col; j++ )
			{
				/* calculate the energy of boxOld + dh; */
				fn_tmat_copy<double> ( sbulkparam->rcpBox, box );
				sbulkparam->rcpBox.val[j][j] += dh;
				fn_obt_rotateProjBox_matrix ( sbulkparam );
				fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
				fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
				fn_calc_bulk_energy ( sbulkparam, energy );
				rham = energy.val[2];

				/* calculate the energy of boxOld - dh; */
				sbulkparam->rcpBox.val[j][j] -= 2.0 * dh;
				fn_obt_rotateProjBox_matrix ( sbulkparam );
				fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
				fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
				fn_calc_bulk_energy ( sbulkparam, energy );
				lham = energy.val[2];

				/* gradient box; */
				gradB.val[j][j] = ( rham - lham ) / ( 2.0*dh );
			}
		}
		else if ( strcmp ( sbulkparam->boxType, "hex" ) == 0 )
		{
			for ( int j = 0; j < col; j++ )
			{
				for ( int i = 0; i < j; i++ )
				{
					/* calculate the energy of boxOld + dh; */
					fn_tmat_copy<double> ( sbulkparam->rcpBox, box );
					sbulkparam->rcpBox.val[i][j] += dh;
					fn_obt_rotateProjBox_matrix ( sbulkparam );
					fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
					fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
					fn_calc_bulk_energy ( sbulkparam, energy );
					rham = energy.val[2];

					/* calculate the energy of boxOld + dh; */
					sbulkparam->rcpBox.val[i][j] -= 2.0 * dh;
					fn_obt_rotateProjBox_matrix ( sbulkparam );
					fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
					fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;
					fn_calc_bulk_energy ( sbulkparam, energy );
					lham = energy.val[2];

					/* gradient box; */
					gradB.val[i][j] = ( rham - lham ) / ( 2.0*dh );
				}
			}
		}

		/* stepping size; */
		if ( iter == 0 )
		{
			step_size = step0;	// the initial stepping size;
		}
		else
		{
			/* 'box' - 'boxOld'; */
			fn_tmat_copy<double> ( boxDiff, box1 );
			fn_tmat_add<double>	 ( boxDiff, box0, 1.0, -1.0 );

			/* 'gradB' - 'gradBOld'; */
			fn_tmat_copy<double> ( gradBDiff, gradB );
			fn_tmat_add<double>  ( gradBDiff, gradBOld, 1.0, -1.0 );

			/* BB stepping size; */
			double bbTmp0, bbTmp1, bbStep;
			if ( bbType == 1 )
			{
				/* BB stepping size (type 1); */
				bbTmp0 = fn_tmat_inner ( boxDiff,   boxDiff );
				bbTmp1 = fn_tmat_inner ( boxDiff, gradBDiff );
			}
			else
			{
				/* BB stepping size (type 2); */
				bbTmp0 = fn_tmat_inner (   boxDiff, gradBDiff );
				bbTmp1 = fn_tmat_inner ( gradBDiff, gradBDiff );
			}
			step_size  = bbTmp0 / bbTmp1;
		}

		/* update computational box; */
		fn_tmat_add<double>  ( box, gradB, 1.0, -step_size );

		/* update; */
		if ( iter > 0 )
			fn_tmat_copy<double> ( box0, box1 );
		fn_tmat_copy<double> ( box1, box );
		fn_tmat_copy<double> ( gradBOld, gradB );

		/* calculate gradient error; */
		err = fn_tmat_maxAbs ( gradB );
		printf("\t\t === > iter %d: step_size = %.5e \t res = % .15e\n", iter, step_size, err );
		iter ++ ;
	}

	printf("rcpBox: \n");
	fn_tmat_print<double> ( box );
	printf("\n");

	/* update the optimized result; */
	fn_tmat_copy<double> ( sbulkparam->rcpBox, box );
	fn_obt_rotateProjBox_matrix ( sbulkparam );
	fn_get_projection_plane	( sbulkparam );	// project plane; R'PBk;
	fn_get_Gsquare			( sbulkparam );	// Gsquare; |R'PBk|^2;

	/* release memory; */
	fn_tmat_free<double> ( box0 );
	fn_tmat_free<double> ( box1 );
	fn_tmat_free<double> ( box );
	fn_tmat_free<double> ( boxDiff );
	fn_tmat_free<double> ( gradBOld );
	fn_tmat_free<double> ( gradB );
	fn_tmat_free<double> ( gradBDiff );
	fn_tvec_free<double> ( energy );
}
