/*! \file	IterMethod.cpp
 *
 * \brief	Iteration methods for the calculation of stable interface structure;
 */

#include "Data.h"
#include "Head.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"

/**
 * \brief	Choose iteration method;
 */
double fn_choose_method		(	stu_system_param		*ssysparam,
									bool				isTest )
{
	printf(" <========== Iteration for the calculation of ");
	printf("stable interface structure ==========> \n\n");
	mytimer_t timer;
	timer.reset();
	timer.start();
	tvec<double> energy = fn_tvec_init<double> ( 3 );
	fn_tvec_setZero<double> ( energy );
	int iterator = 0;

	/* choose iteration method; */
	printf("\t Iteration method is %s.\n", ssysparam->iter_method);
	if ( strcmp(ssysparam->iter_method, "SIS") == 0 )
		iterator = fn_sis_method ( ssysparam, energy, isTest );
	else if ( strcmp(ssysparam->iter_method, "APG") == 0 )
		iterator = fn_apg_method ( ssysparam, energy, ssysparam->tol, isTest );
	else if ( strcmp(ssysparam->iter_method, "nAPG") == 0 )
	{
		printf("\n===> APG step: \n");
		iterator = fn_apg_method	( ssysparam, energy, ssysparam->newton_tol, isTest );
		printf("\n===> Newton step: \n");
		iterator = fn_newton_method ( ssysparam, energy, iterator, isTest );
	}
	else
	{
		printf("Error use 'fn_choose_method'\n");
		printf("Iteration method %s is out of consideration.\n", 
				ssysparam->iter_method);
	}

	/* save stationary value; */
	if ( ssysparam->iter_max > 0 )
		fn_save_system_phase (ssysparam, -1);

	double hamilton = energy.val[2];
	fn_tvec_free<double>		( energy );
	timer.pause();
	printf("\n\t ***** time cost of iteration: %f seconds *****\n\n", timer.get_current_time());

	return hamilton;
}


/**
 * \brief	Semi-implicit method;
 */
int fn_sis_method		(	stu_system_param		*ssysparam,
							tvec<double>			energy,
								bool				isTest )
{
	double	tol			= ssysparam->tol;
	double	tolham		= ssysparam->tolham;
	double	step_size	= ssysparam->step_size;
	int		iter_max	= ssysparam->iter_max;
	int		print_step	= ssysparam->print_step;
	int		save_step	= ssysparam->save_step;

	double	model_xi	= ssysparam->model_xi;
	double	model_tau   = ssysparam->model_tau;
	double	model_gamma = ssysparam->model_gamma;
	double	model_kappa = ssysparam->model_kappa;

	/* transpose 'interact_grad'; */
	mytimer_t timer;
	timer.reset();
	timer.start();
	tCCSmat<double> interact_grad_trans = fn_fast_trans_dCCSmat ( ssysparam->interact_grad );
	timer.pause();
	printf("\n\t ***** time cost of transpose dCCSmat: %f seconds *****\n\n", 
			timer.get_current_time());

	/* memory allocation; */
	int		nd	   = ssysparam->sGJPv->nd;
	int cplxReDofs = ssysparam->cplxReDofs;	
	int	rhsLen	   = nd * cplxReDofs;
	tvec<fftw_complex> rhoJCplxNew	 = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxTrans = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> iter_rhs		 = fn_tvec_init<fftw_complex> ( rhsLen );
	rhoJCplxNew.row   = nd;		rhoJCplxNew.col   = cplxReDofs;
	rhoJCplxTrans.col = nd;		rhoJCplxTrans.row = cplxReDofs;

	/* preparation for mass conservation; */
	tvec<fftw_complex> massVec = fn_tvec_init<fftw_complex> ( nd );
	for ( int j1 = 0; j1 < nd; j1++ )
	{
		int ind0 = j1*cplxReDofs;	// the first element along Fourier direction;
		massVec.val[j1][0] = ssysparam->scbndv->rhoJCplx.val[ind0][0];
		massVec.val[j1][1] = ssysparam->scbndv->rhoJCplx.val[ind0][1];
	}	

	/* obtain the iteration matrix; */
	fn_obt_iter_matrix ( ssysparam, step_size, -1, isTest );

	/* initialization; */
	int iterator = 0;
	int	status = 1;
	double res = 1.0;
	double hamilton, hamiltonOld, diffham;
	char dataName[FILELEN];
	sprintf(dataName, "%s/sys_energy_error.dat", rsltDir);
	FILE *fdata = fopen(dataName, "w");

	int isDecay = 1;
	hamiltonOld = 100.0;
	diffham		= 1.0;

	while ( iterator < iter_max )
	{
		/* calculate free energy and entropy gradient; */
		fn_calc_system_energy ( ssysparam, ssysparam->scbndv->rhoJCplx, energy, true, isTest );
		hamilton = energy.val[2];

		/* save data; */
		fprintf(fdata, "%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", 
						energy.val[0], energy.val[1], energy.val[2], res, isDecay);
		if ( iterator%print_step == 0 )
		{
			double mass = fn_obt_mass ( ssysparam, ssysparam->scbndv->rhoJCplx );
			printf("\nIterator %d: step_size = %.5e \t res = %.10e", iterator, step_size, res);
			printf(" \t diffham = %.10e \t isDecay = %d\n", diffham, isDecay);
			printf("\t Laplace term: % .10e \t Nonlinear term: % .10e,", 
					energy.val[0], energy.val[1]);
			printf("\t hamilton = % .10e\n", energy.val[2]);
			printf("\t Mass: %.10e\n", mass);
		}

		if ( (iterator > 0 ) && (iterator % save_step) == 0 )
			fn_save_system_phase (ssysparam, iterator);		// save data;

		diffham = hamilton - hamiltonOld;
		hamiltonOld = hamilton;

		/* set flag to check if the hamilton energy is decay; */
		if ( diffham < 0 )
			isDecay = 1;
		else
			isDecay = 0;
	
		/* generate the right term for iteration; */
		fn_obt_iter_rhs ( ssysparam, ssysparam->scbndv->rhoJCplx,
						ssysparam->sfftv->gradient, iterator, isTest );
		memcpy ( iter_rhs.val, ssysparam->iter_rho_rhs.val, sizeof(fftw_complex) * rhsLen );
		fn_tvec_add_complex ( iter_rhs, ssysparam->iter_entropy_rhs, 1.0/step_size, -1.0 );

		/* calculate the new 'rhoJCplx'; */
		status = fn_umfpack_complex_solver ( ssysparam->iter_matrix, iter_rhs, rhoJCplxNew );

		/* mass conservation; */
		for ( int j1 = 0; j1 < nd; j1++ )
		{
			int ind0 = j1;	// the first element along Fourier direction;
			rhoJCplxNew.val[ind0][0] = massVec.val[j1][0];
			rhoJCplxNew.val[ind0][1] = massVec.val[j1][1];
		}

		/* transpose; */
		fn_tvec_trans_complex ( rhoJCplxTrans, rhoJCplxNew, nd, cplxReDofs );

		/* calculate gradient error; */
		fn_cvec_multiply_dCCSmat ( interact_grad_trans, rhoJCplxNew, ssysparam->grad_err );
		fn_tvec_add_complex ( ssysparam->grad_err, ssysparam->iter_entropy_rhs, model_xi, 1.0 );
		for ( int i = 0; i < nd; i++ )
			fn_complex_setZero ( ssysparam->grad_err.val[i] );
		res = fn_tvec_maxAbs_complex ( ssysparam->grad_err );

		/* update 'rhoJCplx'; * transpose for consistency; */
		memcpy ( ssysparam->scbndv->rhoJCplx.val, rhoJCplxTrans.val,
					sizeof(fftw_complex) * rhoJCplxTrans.len );

		iterator ++;
//		if ( (res < tol && iterator > 2) || fabs(diffham) < 1.0e-15 )
//		if ( res < tol && iterator > 2 )
//		if ( res < tol )
		if ( res < tol || fabs(diffham) < tolham )
			break;
	}

	/* output data about the stationary state; */
	fn_calc_system_energy ( ssysparam, ssysparam->scbndv->rhoJCplx, energy, false, isTest );
	hamilton = energy.val[2];
	fprintf(fdata, "%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", 
			energy.val[0], energy.val[1], energy.val[2], res, isDecay);
	printf("\n\t ***** the interface system\n");
	printf("\t ***** iterator = %d \t step size = %.3e\n", iterator, step_size);
	printf("\t ***** res = % .10e \t diffham = % .10e\n", res, diffham);
	printf("\t ***** Laplace term: % .15e,\t Nonlinear term: % .15e\n", 
			energy.val[0], energy.val[1]);
	printf("\t ***** hamilton = % .20e\n\n", energy.val[2]);
	fclose(fdata);

	/* Releases memory; */
	fn_tCCSmat_free<double>		( interact_grad_trans );
	fn_tvec_free<fftw_complex>	( rhoJCplxNew );
	fn_tvec_free<fftw_complex>	( rhoJCplxTrans );
	fn_tvec_free<fftw_complex>	( iter_rhs );
	fn_tvec_free<fftw_complex>  ( massVec );

	return iterator;
}


/**
 * \brief	Adaptive accelerated Bregman proximal gradient method (P2 case);
 *			Unpublished name is accelerated proximal gradient method (APG);
 */
int fn_apg_method		(	stu_system_param		*ssysparam,
							tvec<double>			energy,
								double				tol,
								bool				isTest )
{
	double	tolham		= ssysparam->tolham;
	int		bbType	    = 1;						// the type of BB stepping size;
	double	step0	    = ssysparam->step_size;	// initial step size;
	double	step_min    = ssysparam->step_min;	// the lower bound of step size;
	double	step_max    = ssysparam->step_max;	// the upper bound of step size;
	int		iter_max    = ssysparam->iter_max;
	int		print_step  = ssysparam->print_step;
	int		save_step   = ssysparam->save_step;

	double	model_xi	= ssysparam->model_xi;
	double	model_tau   = ssysparam->model_tau;
	double	model_gamma = ssysparam->model_gamma;
	double	model_kappa = ssysparam->model_kappa;

	/* transpose 'interact_grad'; */
	mytimer_t timer;
	timer.reset();
	timer.start();
	tCCSmat<double> interact_grad_trans = fn_fast_trans_dCCSmat ( ssysparam->interact_grad );
	timer.pause();
	printf("\n\t ***** time cost of transpose dCCSmat: %f seconds *****\n\n", 
			timer.get_current_time());

	/* memory allocation; */
	int		nd	   = ssysparam->sGJPv->nd;
	int		xlen   = ssysparam->sGJPv->xlen;
	int cplxReDofs = ssysparam->cplxReDofs;	
	int	rhsLen	   = nd * cplxReDofs;
	tvec<fftw_complex> rhoJCplxNew   = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxTrans = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> iter_rhs		 = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> gradientOld	 = fn_tvec_init<fftw_complex> ( xlen * cplxReDofs );
	rhoJCplxNew.row   = nd;		rhoJCplxNew.col   = cplxReDofs;
	rhoJCplxTrans.col = nd;		rhoJCplxTrans.row = cplxReDofs;
	gradientOld.row	  = xlen;	gradientOld.col	  = cplxReDofs;

	/* temporary variables for APG method; */
	tvec<fftw_complex> rhoJCplxTmp0  = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxTmp1  = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxInp	 = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> rhoJCplxDiff  = fn_tvec_init<fftw_complex> ( cplxReDofs );
	tvec<fftw_complex> gradientDiff  = fn_tvec_init<fftw_complex> ( cplxReDofs );
	memcpy( rhoJCplxTmp0.val, ssysparam->scbndv->rhoJCplx.val, sizeof(fftw_complex) * rhsLen );
	memcpy( rhoJCplxInp.val,  ssysparam->scbndv->rhoJCplx.val, sizeof(fftw_complex) * rhsLen );
	rhoJCplxTmp0.row   = nd;		rhoJCplxTmp0.col   = cplxReDofs;
	rhoJCplxTmp1.row   = nd;		rhoJCplxTmp1.col   = cplxReDofs;
	rhoJCplxInp.row    = nd;		rhoJCplxInp.col    = cplxReDofs;

	/* preparation for mass conservation; */
	tvec<fftw_complex> massVec = fn_tvec_init<fftw_complex> ( nd );
	for ( int j1 = 0; j1 < nd; j1++ )
	{
		int ind0 = j1*cplxReDofs;	// the first element along Fourier direction;
		massVec.val[j1][0] = ssysparam->scbndv->rhoJCplx.val[ind0][0];
		massVec.val[j1][1] = ssysparam->scbndv->rhoJCplx.val[ind0][1];
	}	

	/* initialization; */
	int iterator = 0;
	double step_size = step0;
	int	status = 1;
	double res = 1.0;
	double hamilton, hamiltonOld, diffham;
	char dataName[FILELEN];
	sprintf(dataName, "%s/sys_energy_error.dat", rsltDir);
	FILE *fdata = fopen(dataName, "w");

	int reCount = 0;
	int isDecay = 1;
	hamiltonOld = 100.0;
	diffham		= 1.0;
	tvec<double> energyNew	= fn_tvec_init<double> ( 3 );
	fn_tvec_setZero<double> ( energyNew );

	/* innSMatd0JJ; */
	tCCSmat<double> G0innd0JJ = fn_tensor_diag_dCCSmat ( 
			ssysparam->sfftv->Gsquare, ssysparam->sGJPv->innSMatd0JJ, 0 );
	ssysparam->iter_matrix	  = fn_tCCSmat_init<double> (
			G0innd0JJ.row, G0innd0JJ.col, ssysparam->interact_grad.nnz + G0innd0JJ.nnz );

	/* parameters for Nesterov accelerated method; */
	double	theta = 1.0;
	double	q	  = 0.0;

	/* parameters for inexact linear search; */
	bool	isBreak = true;
	double	rho		= 0.5 * (sqrt(5) - 1.0);
	double	delta	= 1.0e-14;

	while ( iterator < iter_max )
	{
		/* calculate free energy and entropy gradient; */
		fn_calc_system_energy ( ssysparam, ssysparam->scbndv->rhoJCplx, energy, true, isTest );
		hamilton = energy.val[2];
	
		/* estimate the time stepping size by BB method; */
		if ( iterator > 0 )
		{
			/* difference of order parameters; */
			fn_diff_weight_rhoJCplx ( ssysparam, rhoJCplxTmp1, 
										rhoJCplxTmp0, rhoJCplxDiff );
			fn_diff_weight_gradient ( ssysparam, ssysparam->sfftv->gradient, 
										gradientOld, gradientDiff );
			if ( isTest && iterator == 1 )
			{
				char fwFile[FILELEN];
				sprintf(fwFile, "%s/rhoJCplxTmp0.dat", rsltDir);
				fn_tvec_save_complex ( rhoJCplxTmp0, fwFile );
				sprintf(fwFile, "%s/rhoJCplxTmp1.dat", rsltDir);
				fn_tvec_save_complex ( rhoJCplxTmp1, fwFile );
				sprintf(fwFile, "%s/rhoJCplxDiff.dat", rsltDir);
				fn_tvec_save_complex ( rhoJCplxDiff, fwFile );

				sprintf(fwFile, "%s/gradientOld.dat", rsltDir);
				fn_tvec_save_complex ( gradientOld, fwFile );
				sprintf(fwFile, "%s/gradient.dat", rsltDir);
				fn_tvec_save_complex ( ssysparam->sfftv->gradient, fwFile );
				sprintf(fwFile, "%s/gradientDiff.dat", rsltDir);
				fn_tvec_save_complex ( gradientDiff, fwFile );
			}
	
			/* calculate BB stepping size; */
			fftw_complex bbTmp0, bbTmp1, bbStep;
			if ( bbType == 1 )
			{
				/* BB stepping size (type 1); */
				fn_tvec_dotMultiplySum_complex ( rhoJCplxDiff, rhoJCplxDiff, bbTmp0 );
				fn_tvec_dotMultiplySum_complex ( rhoJCplxDiff, gradientDiff, bbTmp1 );
			}
			else
			{
				/* BB stepping size (type 2); */
				fn_tvec_dotMultiplySum_complex ( rhoJCplxDiff, gradientDiff, bbTmp0 );
				fn_tvec_dotMultiplySum_complex ( gradientDiff, gradientDiff, bbTmp1 );
			}
			fn_complex_divide ( bbTmp0, bbTmp1, bbStep ); // bbStep = bbTmp0 / bbTmp1;
			step_size  = fabs ( bbStep[0] );

			if ( isTest )
			{
				printf("\nBB step : [%.5e, %.5e] / [%.5e, %.5e] = [%.5e, %.5e].\n", 
						bbTmp0[0], bbTmp0[1], bbTmp1[0], bbTmp1[1], bbStep[0], bbStep[1]);
			}

			/* 'step_size' must belong to the region [step_min, step_max]; */
			step_size  = ( step_size > step_min ? step_size : step_min );
			step_size  = ( step_size < step_max ? step_size : step_max );
			isBreak = false;
		}

		/* generate the right term for iteration; */
		fn_obt_iter_rhs ( ssysparam, ssysparam->scbndv->rhoJCplx,
						ssysparam->sfftv->gradient, iterator, isTest );

		/* search time stepping size for energy decay; */
		while ( true )
		{
			/* adopt the minimal stepping size and then break recurrence; */
			if ( step_size <= step_min )
			{
				step_size  = step_min;
				isBreak = true;
			}

			/* obtain the iteration matrix; */
			//fn_obt_iter_matrix ( ssysparam, step_size, iterator, isTest );
			fn_tCCSmat_free<double> ( ssysparam->iter_matrix );
			ssysparam->iter_matrix	  = fn_add_dCCSmat ( 
					ssysparam->interact_grad, G0innd0JJ, model_xi, 1.0/step_size );

			/* generate the right term for iteration; */
			memcpy ( iter_rhs.val, ssysparam->iter_rho_rhs.val, sizeof(fftw_complex) * rhsLen );
			fn_tvec_add_complex ( iter_rhs, ssysparam->iter_entropy_rhs, 1.0/step_size, -1.0 );

			/* calculate the new 'rhoJCplx'; */
			status = fn_umfpack_complex_solver ( 
						ssysparam->iter_matrix, iter_rhs, rhoJCplxNew );

			/* mass conservation; */
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				int ind0 = j1;	// the first element along Fourier direction;
				rhoJCplxNew.val[ind0][0] = massVec.val[j1][0];
				rhoJCplxNew.val[ind0][1] = massVec.val[j1][1];
			}

			/* transpose; */
			fn_tvec_trans_complex ( rhoJCplxTrans, rhoJCplxNew, nd, cplxReDofs );

			/* compute energy about the temporal order parameter; */
			fn_calc_system_energy ( ssysparam, rhoJCplxTrans, energyNew, false, isTest );
			
			/* ensure enough energy decay; */
			fn_diff_weight_rhoJCplx ( ssysparam, rhoJCplxTrans, 
								ssysparam->scbndv->rhoJCplx, rhoJCplxDiff );
			double energy_decay = fn_tvec_norm_complex ( rhoJCplxDiff );
			energy_decay = energy.val[2] - energyNew.val[2] - delta * energy_decay;

			if ( isTest )
			{
				printf("step size = %.6e.\n", step_size);
				printf("hamilton = %.6e \t hamilton temp = %.6e\n", 
						energy.val[2], energyNew.val[2]);
				printf("energy decay = %.6e.\n\n", energy_decay);
			}
			if ( energy_decay <= 0 && !isBreak )
				step_size *= rho;
			else
				break;
		}

		if ( (iterator > 0 ) && (iterator % save_step) == 0 )
			fn_save_system_phase (ssysparam, iterator);		// save data;

		/* copy gradient values; must be before gradient error calculation; */
		memcpy ( gradientOld.val, ssysparam->sfftv->gradient.val, 
					sizeof(fftw_complex) * gradientOld.len );

		/* calculate gradient error; */
		fn_cvec_multiply_dCCSmat ( interact_grad_trans, rhoJCplxNew, ssysparam->grad_err );
		fn_tvec_add_complex ( ssysparam->grad_err, ssysparam->iter_entropy_rhs, model_xi, 1.0 );
		for ( int i = 0; i < nd; i++ )
			fn_complex_setZero ( ssysparam->grad_err.val[i] );
		res = fn_tvec_maxAbs_complex ( ssysparam->grad_err );

		if ( isTest )
		{
			char fwFile[FILELEN];
			sprintf(fwFile, "%s/grad_err%d.dat", rsltDir, iterator);
			fn_tvec_save_complex ( ssysparam->grad_err, fwFile );
		}
		
		/* compute the difference between the hamilton energies of two adjoin steps; */
		diffham = hamilton - hamiltonOld;
		hamiltonOld = hamilton;

		/* save data; */
		fprintf(fdata, "%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", 
						energy.val[0], energy.val[1], energy.val[2], res, isDecay);
		if ( iterator%print_step == 0 )
		{
			double mass = fn_obt_mass ( ssysparam, ssysparam->scbndv->rhoJCplx );
			printf("\nIterator %d: step_size = %.5e \t res = %.10e", iterator, step_size, res);
			printf(" \t diffham = %.10e \t isDecay = %d\n", diffham, isDecay);
			printf("\t Laplace term: % .10e \t Nonlinear term: % .10e,", 
					energy.val[0], energy.val[1]);
			printf("\t hamilton = % .10e\n", energy.val[2]);
			printf("\t Mass: %.10e\n", mass);
		}

		/* restart if energy is not decay; */
		if ( diffham < 0 )
		{
			/* set flag to check if the hamilton energy is decay; */
			isDecay = 1;
            /* the parameter for Lagrange extrapolation; */
			double theta2 = pow(theta, 2);
            double thetaTmp = -0.5*(theta2-q) + sqrt(0.25*pow(theta2-q,2) + theta2);
            double beta = theta*(1.0-theta) / (theta2+thetaTmp);
            theta = thetaTmp;
            /* Interpolation; */
			memcpy ( rhoJCplxInp.val, rhoJCplxTrans.val, sizeof(fftw_complex) * rhoJCplxTrans.len );
			fn_tvec_add_complex ( rhoJCplxInp, ssysparam->scbndv->rhoJCplx, 1.0+beta, -beta ); 
		}
		else
		{
			/* set flag to check if the hamilton energy is decay; */
			isDecay = 0;
			/* restart; */
            theta = 1.0;
			memcpy ( rhoJCplxInp.val, rhoJCplxTrans.val, sizeof(fftw_complex) * rhoJCplxTrans.len );
			printf("--> restart\n");
			reCount ++ ;
		}

		/* update 'rhoJCplx'; * transpose for consistency; */
		memcpy ( ssysparam->scbndv->rhoJCplx.val, rhoJCplxTrans.val,
					sizeof(fftw_complex) * rhoJCplxTrans.len );

		/* update 'rhoJCplxTmp0', 'rhoJCplxTmp1'; */
		if ( iterator > 0 )
			memcpy ( rhoJCplxTmp0.val, rhoJCplxTmp1.val, sizeof(fftw_complex) * rhoJCplxTmp0.len );
		memcpy ( rhoJCplxTmp1.val, rhoJCplxInp.val, sizeof(fftw_complex) * rhoJCplxTmp1.len );

		iterator ++;
//		if ( (res < tol && iterator > 2) || fabs(diffham) < 1.0e-15 )
//		if ( res < tol && iterator > 2 )
		if ( res < tol || fabs(diffham) < tolham )
			break;
		if ( reCount > 20 )
			break;
	}

	/* output data about the stationary state; */
	fn_calc_system_energy ( ssysparam, ssysparam->scbndv->rhoJCplx, energy, false, isTest );
	hamilton = energy.val[2];
	fprintf(fdata, "%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", 
			energy.val[0], energy.val[1], energy.val[2], res, isDecay);
	printf("\n\t ***** the interface system\n");
	printf("\t ***** iterator = %d\n", iterator);
	printf("\t ***** res = % .10e \t diffham = % .10e\n", res, diffham);
	printf("\t ***** Laplace term: % .15e,\t Nonlinear term: % .15e\n", 
			energy.val[0], energy.val[1]);
	printf("\t ***** hamilton = % .20e\n\n", energy.val[2]);
	fclose(fdata);

	/* Releases memory; */
	fn_tvec_free<double>		( energyNew );
	fn_tCCSmat_free<double>		( interact_grad_trans );
	fn_tvec_free<fftw_complex>	( rhoJCplxNew );
	fn_tvec_free<fftw_complex>	( rhoJCplxTrans );
	fn_tvec_free<fftw_complex>	( iter_rhs );
	fn_tvec_free<fftw_complex>	( gradientOld );
	fn_tCCSmat_free<double>		( G0innd0JJ );
	/* temporary variables for APG method; */
	fn_tvec_free<fftw_complex>	( rhoJCplxTmp0 );
	fn_tvec_free<fftw_complex>	( rhoJCplxTmp1 );
	fn_tvec_free<fftw_complex>	( rhoJCplxInp  );
	fn_tvec_free<fftw_complex>	( rhoJCplxDiff );
	fn_tvec_free<fftw_complex>	( gradientDiff );
	fn_tvec_free<fftw_complex>	( massVec );

	return iterator;
}


/**
 * \brief	Newton PCG method;
 */
int fn_newton_method	(	stu_system_param		*ssysparam,
							tvec<double>			energy,
								int					iterator,
								bool				isTest )
{
	double	tol			= ssysparam->tol;
	double	tolham		= ssysparam->tolham;
	double	step0		= ssysparam->newton_step_size;
	int		iter_max	= ssysparam->iter_max;
	int		print_step	= ssysparam->print_step;
	int		save_step	= ssysparam->save_step;
	double	mu_para		= ssysparam->pcg_mu_para;

	double	model_xi	= ssysparam->model_xi;
	double	model_tau   = ssysparam->model_tau;
	double	model_gamma = ssysparam->model_gamma;
	double	model_kappa = ssysparam->model_kappa;

	/* transpose 'interact_grad'; */
	mytimer_t timer;
	timer.reset();
	timer.start();
	tCCSmat<double> interact_grad_trans = fn_fast_trans_dCCSmat ( ssysparam->interact_grad );
	timer.pause();
	printf("\n\t ***** time cost of transpose dCCSmat: %f seconds *****\n\n", 
			timer.get_current_time());

	/* memory allocation; */
	int		nd	   = ssysparam->sGJPv->nd;
	int cplxReDofs = ssysparam->cplxReDofs;	
	int	rhsLen	   = nd * cplxReDofs;
	tvec<fftw_complex> rhoJCplxNew		  = fn_tvec_init<fftw_complex> ( rhsLen );
	tvec<fftw_complex> PCG_direction	  = fn_tvec_init<fftw_complex> ( rhsLen ); // x^{i} of PCG;
	tvec<fftw_complex> PCG_directionTrans = fn_tvec_init<fftw_complex> ( rhsLen );
	rhoJCplxNew.col			= nd;		rhoJCplxNew.row			= cplxReDofs;
	PCG_direction.col		= nd;		PCG_direction.row		= cplxReDofs;
	PCG_directionTrans.row	= nd;		PCG_directionTrans.col	= cplxReDofs;

	tvec<fftw_complex> rhoJCplxTrans	  = fn_tvec_init<fftw_complex> ( rhsLen );
	rhoJCplxTrans.row		= nd;		rhoJCplxTrans.col		= cplxReDofs;

	/* preparation for mass conservation; */
	tvec<fftw_complex> massVec = fn_tvec_init<fftw_complex> ( nd );
	for ( int j1 = 0; j1 < nd; j1++ )
	{
		int ind0 = j1*cplxReDofs;	// the first element along Fourier direction;
		massVec.val[j1][0] = ssysparam->scbndv->rhoJCplx.val[ind0][0];
		massVec.val[j1][1] = ssysparam->scbndv->rhoJCplx.val[ind0][1];
	}	

	/* initialization; */
	double step_size = step0;
	int subIter = 0;
	int	status = 1;
	double mu, delta;
	double hamilton, diffham;
	double hamiltonOld = energy.val[2];
	double res = fn_tvec_maxAbs_complex ( ssysparam->grad_err );
	char dataName[FILELEN];
	sprintf(dataName, "%s/sys_energy_error.dat", rsltDir);
	FILE *fdata = fopen(dataName, "a");

	int isDecay = 3;
	diffham		= 1.0;

	if ( isTest )
	{
		printf("mu_para = %.10e \t res = %.10e.\n", mu_para, res);
		char fwFile[FILELEN];
		sprintf(fwFile, "%s/newton_rhoJCplx%d.dat", rsltDir, 0);
		fn_tvec_save_complex ( ssysparam->scbndv->rhoJCplx, fwFile );
		sprintf(fwFile, "%s/newton_grad_err%d.dat", rsltDir, 0);
		fn_tvec_save_complex ( ssysparam->grad_err, fwFile );
	}

	/* obtain Hessian matrix; */
	fn_calc_system_hessian ( ssysparam, ssysparam->scbndv->rhoJCplx, iterator, isTest );

	while ( iterator < iter_max )
	{
		/* precondition; */
		double tmp = mu_para * res;	
		mu = tmp < 100 ? tmp : 100;
		mu = mu > 1.0e-20 ? mu : 1.0e-20;
		delta = fn_tvec_maxAbs_complex ( ssysparam->hessian );

		/* Pre-conditional Projected Conjugate Gradient method; */
		subIter = fn_calc_system_PCG ( ssysparam, PCG_direction, 
					interact_grad_trans, mu, delta, iterator, false );

		/* transpose; */
		fn_tvec_trans_complex ( PCG_directionTrans, PCG_direction, nd, cplxReDofs );

		/* rho = - \<x^{i}, b\> / \|x^{i}\|^{2};		b = grad_err; */
		fftw_complex rhoTmp0;
		fn_tvec_dotMultiplySum_complex ( ssysparam->grad_err, PCG_direction, rhoTmp0 );
		double rhoTmp1 = fn_tvec_norm_complex ( PCG_direction );
		double rho = -rhoTmp0[0] / pow(rhoTmp1,2);

		/* update 'mu_para'; */
		if ( fabs(rho) > 0.1 )
			mu_para = 0.1 * mu_para;
		else if ( fabs(rho) < 0.01 )
			mu_para = 10.0 * mu_para;
		else
			mu_para = 0.5 * mu_para;

		if ( isTest )
		{
			printf("mu = %.10e \t delta = %.10e \t rho = %.10e.\n", mu, delta, rho);
			printf("rhoTmp0 = %.5e, %.5e \t rhoTmp1 = %.5e.\n",
					rhoTmp0[0], rhoTmp0[1], rhoTmp1);
			printf("hamilton old = %.10e.\n\n", hamiltonOld);
			char fwFile[FILELEN];
			sprintf(fwFile, "%s/PCG_direction%d.dat", rsltDir, iterator);
			fn_tvec_save_complex ( PCG_direction, fwFile );
		}

		/* update 'rhoJCplx' and free energies; */
		step_size = step0;
		while ( step_size > 1.0e-6 )
		{
			/* update 'rhoJCplx' : rhoJCplx + step_size * PCG_directionTrans; */
			memcpy ( rhoJCplxNew.val, ssysparam->scbndv->rhoJCplx.val, 
						sizeof(fftw_complex) * rhsLen );
			fn_tvec_add_complex ( rhoJCplxNew, PCG_directionTrans, 1.0, step_size );
			for ( int j1 = 0; j1 < nd; j1++ )
			{
				int ind0 = j1*cplxReDofs;	// the first element along Fourier direction;
				rhoJCplxNew.val[ind0][0] = massVec.val[j1][0];
				rhoJCplxNew.val[ind0][1] = massVec.val[j1][1];
			}

			/* calculate free energy and entropy gradient; */
			fn_calc_system_energy ( ssysparam, rhoJCplxNew, energy, true, isTest );
			hamilton = energy.val[2];

			/* if decay; */
			diffham = hamilton - hamiltonOld;
			double diffhamCheck = 1.0e-6 * step_size * rhoTmp0[0];
			if ( diffham < diffhamCheck )
				break;
			else
				step_size = 0.5 * step_size;
//			printf("energy: %.10e \t %.10e \t %.10e.\n", energy.val[0], energy.val[1], energy.val[2]);
//			printf("diffhamCheck = %.5e \t step_size = %.5e \t diffham = %.5e.\n\n", diffhamCheck, step_size, diffham);
		}
		hamilton = energy.val[2];
		diffham = hamilton - hamiltonOld;
		hamiltonOld = hamilton;
		if ( diffham < 0 )
			isDecay = 3;
		else
			isDecay = 2;

		/* calculate gradient error; */
		/*
		memcpy ( ssysparam->grad_err.val, rhoJCplxNew.val, 
					sizeof(fftw_complex) * rhoJCplxNew.len );
		fn_tvec_add_complex ( ssysparam->grad_err, ssysparam->scbndv->rhoJCplx, 1.0, -1.0 );
		res = fn_tvec_maxAbs_complex ( ssysparam->grad_err );
		res /= step_size;
		*/

		/* calculate gradient error; */
		//fn_calc_system_energy ( ssysparam, rhoJCplxNew, energy, true, isTest );
		fn_obt_iter_rhs	( ssysparam, rhoJCplxNew, ssysparam->sfftv->gradient, iterator, isTest );
		fn_tvec_trans_complex ( rhoJCplxTrans, rhoJCplxNew, cplxReDofs, nd );
		fn_cvec_multiply_dCCSmat ( interact_grad_trans, rhoJCplxTrans, ssysparam->grad_err );
		fn_tvec_add_complex ( ssysparam->grad_err, ssysparam->iter_entropy_rhs, model_xi, 1.0 );
		for ( int i = 0; i < nd; i++ )
			fn_complex_setZero ( ssysparam->grad_err.val[i] );	// mass conservation;
		res = fn_tvec_maxAbs_complex ( ssysparam->grad_err );

		/* update 'rhoJCplx';	transpose for consistency; */
		memcpy ( ssysparam->scbndv->rhoJCplx.val, rhoJCplxNew.val,
					sizeof(fftw_complex) * rhoJCplxNew.len );

		/* update Hessian matrix;	interact_grad_trans; */
		fn_calc_system_hessian ( ssysparam, rhoJCplxNew, iterator, isTest );

		if ( isTest )
		{
			char fwFile[FILELEN];
			sprintf(fwFile, "%s/Newton_rhoJCplxNew%d.dat", rsltDir, iterator);
			fn_tvec_save_complex ( rhoJCplxNew, fwFile );
			sprintf(fwFile, "%s/Newton_iter_entropy_rhs%d.dat", rsltDir, iterator);
			fn_tvec_save_complex ( ssysparam->iter_entropy_rhs, fwFile );
			sprintf(fwFile, "%s/Newton_grad_err%d.dat", rsltDir, iterator);
			fn_tvec_save_complex ( ssysparam->grad_err, fwFile );
		}

		/* save data; */
		fprintf(fdata, "%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", 
						energy.val[0], energy.val[1], energy.val[2], res, isDecay);
		if ( iterator%print_step == 0 )
		{
			double mass = fn_obt_mass ( ssysparam, ssysparam->scbndv->rhoJCplx );
			printf("\nIterator %d: step_size = %.5e \t res = %.10e", iterator, step_size, res);
			printf(" \t diffham = %.10e \t isDecay = %d\n", diffham, isDecay);
			printf("\t Laplace term: % .10e \t Nonlinear term: % .10e,", 
					energy.val[0], energy.val[1]);
			printf("\t hamilton = % .10e\n", energy.val[2]);
			printf("\t Mass: %.10e\n", mass);
		}

		if ( (iterator > 0 ) && (iterator % save_step) == 0 )
			fn_save_system_phase (ssysparam, iterator);		// save data;

		iterator ++;
//		if ( (res < tol && iterator > 2) || fabs(diffham) < 1.0e-15 )
//		if ( res < tol && iterator > 2 )
//		if ( res < tol )
		if ( res < tol || fabs(diffham) < tolham )
			break;
	}

	/* output data about the stationary state; */
	fn_calc_system_energy ( ssysparam, ssysparam->scbndv->rhoJCplx, energy, false, isTest );
	hamilton = energy.val[2];
	fprintf(fdata, "%+.15E\t%+.15E\t%+.20E\t%+.15E\t%d\n", 
			energy.val[0], energy.val[1], energy.val[2], res, isDecay);
	printf("\n\t ***** the interface system\n");
	printf("\t ***** iterator = %d \t step size = %.3e\n", iterator, step_size);
	printf("\t ***** res = % .10e \t diffham = % .10e\n", res, diffham);
	printf("\t ***** Laplace term: % .15e,\t Nonlinear term: % .15e\n", 
			energy.val[0], energy.val[1]);
	printf("\t ***** hamilton = % .20e\n\n", energy.val[2]);
	fclose(fdata);

	/* Releases memory; */
	fn_tCCSmat_free<double>		( interact_grad_trans );
	fn_tvec_free<fftw_complex>	( rhoJCplxNew	);
	fn_tvec_free<fftw_complex>	( rhoJCplxTrans );
	fn_tvec_free<fftw_complex>	( PCG_direction );
	fn_tvec_free<fftw_complex>	( PCG_directionTrans );
	fn_tvec_free<fftw_complex>	( massVec );

	return iterator;
}
