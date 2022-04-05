/*! \file	MainCollect.cpp
 *
 *  \brief	Collection of main code to realize different functions;
 *
 */

#include "Head.h"
#include "Data.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"
#include "umfpack.h"

/*---------------------------------*/
/*--       Main Functions        --*/
/*---------------------------------*/

/**
 * \brief	Calculation of the stable interface structure;
 */
double fn_main_stable_interface	( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* true: save;	false: not save; */
	/* whether save densities; */
	sbulkparam1->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	sbulkparam1->srebndv->isSave = true; // 'fn_common_Fourier_GJP';
	sbulkparam2->srebndv->isSave = true; // 'fn_common_Fourier_GJP';
	/* variables for debug; */
	ssysparam->sGJPv->isTest	 = false; // 'fn_obt_system_gen_Jac_poly';
	sbulkparam1->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';
	sbulkparam1->srebndv->isTest = false; // 'fn_common_Fourier_GJP';
	sbulkparam2->srebndv->isTest = false; // 'fn_common_Fourier_GJP';
	ssysparam->scbndv->isTest	 = false; // 'fn_system_initial_value';
	bool iter_prepare_isTest	 = false; // 'fn_iter_prepare';
	bool iter_method_isTest		 = false; // 'fn_choose_method';


	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	/* obtain the stable bulk phases; */
	double hamilton1 = fn_obt_stable_bulk_phase ( sbulkparam1 );
	double hamilton2 = fn_obt_stable_bulk_phase ( sbulkparam2 );

	/* update 'x_range'; */
	if ( strcmp ( ssysparam->x_range_type, "direct" ) == 0 )
	{
		printf("x_range_type = %s.\n", ssysparam->x_range_type);
	}
	else if ( strcmp ( ssysparam->x_range_type, "cubePlane" ) == 0 )
	{
		double angle = fabs ( sbulkparam1->rotate_angle.val[0] ) / 180.0 * PI;
		double tmp1 = sin(angle) + cos(angle);
		tmp1 *= sbulkparam1->dirBox.val[0][0];
		ssysparam->x_range *= tmp1;
	}
	else
	{
		printf("Error x_range_type %s.\n", ssysparam->x_range_type);
	}
	ssysparam->sGJPv->x_range = ssysparam->x_range;

	/* prepare general Jacobi polynomials; */
	fn_obt_system_gen_Jac_poly ( ssysparam->sGJPv );

	/* project rhoCplx to the space of GJPs; */
	fn_project_Fourier_GJP ( sbulkparam1, ssysparam );
	fn_project_Fourier_GJP ( sbulkparam2, ssysparam );

	/* obtain the common 'rotateProjBoxMat'; */
	int status = fn_obt_commom_rotateProjBoxMat ( sbulkparam1, sbulkparam2, ssysparam );
	if ( status == 0 )
		return 0.0;

	/* save parameters (optimal computational box); */
	fn_param_save ( sbulkparam1, sbulkparam2, ssysparam, "opt" );

	/* to a common space; */
	fn_common_Fourier_GJP   ( sbulkparam1, sbulkparam2, ssysparam	);

	/* connect two bulk phases to construct initial value; */
	fn_system_initial_value ( sbulkparam1, sbulkparam2, ssysparam	);

	/* iteration preparation; */
	fn_iter_prepare	( ssysparam, iter_prepare_isTest );

	/* iteration by semi-implicit scheme (SIS); */
	double hamilton = fn_choose_method ( ssysparam, iter_method_isTest );

	/* memory free about the interface system; */
	fn_GJP_memory_free	  ( ssysparam->sGJPv );
	fn_bnd_memory_free	  ( ssysparam->scbndv, "part" );
	fn_system_memory_free ( ssysparam );

	return hamilton;
}


/**
 * \brief	Calculation of the stable bulk phase;
 */
double fn_main_stable_bulk ( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	/* obtain the stable bulk phases; */
	double hamilton = 0.0;
	if ( strcmp( main_type, "stable_bulk1" ) == 0 )
	{
		hamilton = fn_obt_stable_bulk_phase ( sbulkparam1 );
		fn_bulk_memory_free ( sbulkparam1, "project bulk" );
	}
	else if ( strcmp( main_type, "stable_bulk2" ) == 0 )
	{
		hamilton = fn_obt_stable_bulk_phase ( sbulkparam2 );
		fn_bulk_memory_free ( sbulkparam2, "project bulk" );
	}
	else
	{
		printf("Error use 'fn_main_stable_bulk'.\n");
		printf("'main_type' should be 'stable_bulk1' or 'stable_bulk2'.\n");
	}

	/* release memory; */
	fn_bulk_memory_free	( sbulkparam1, "bulk param" );
	fn_bulk_memory_free	( sbulkparam2, "bulk param" );

	return hamilton;
}


/**
 * \brief	Projection of the stable bulk phase;
 */
void fn_main_project_bulk ( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* true: save;	false: not save; */
	/* whether save densities; */
	sbulkparam1->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	/* variables for debug; */
	ssysparam->sGJPv->isTest	 = false; // 'fn_obt_system_gen_Jac_poly';
	sbulkparam1->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';


	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	double hamilton = 0.0;
	if ( strcmp( main_type, "proj_bulk1" ) == 0 )
	{
		/* obtain the stable bulk phases; */
		hamilton = fn_obt_stable_bulk_phase ( sbulkparam1 );

		/* prepare general Jacobi polynomials; */
		fn_obt_system_gen_Jac_poly ( ssysparam->sGJPv );

		/* project rhoCplx to the space of GJPs; */
		fn_project_Fourier_GJP ( sbulkparam1, ssysparam );

		/* release memory; */
		fn_bulk_memory_free ( sbulkparam1, "project bulk" );
		fn_bnd_memory_free  ( sbulkparam1->sbndv, "total" );
	}
	else if ( strcmp( main_type, "proj_bulk2" ) == 0 )
	{
		/* obtain the stable bulk phases; */
		hamilton = fn_obt_stable_bulk_phase ( sbulkparam2 );

		/* prepare general Jacobi polynomials; */
		fn_obt_system_gen_Jac_poly ( ssysparam->sGJPv );

		/* project rhoCplx to the space of GJPs; */
		fn_project_Fourier_GJP ( sbulkparam2, ssysparam );

		/* release memory; */
		fn_bulk_memory_free ( sbulkparam2, "project bulk" );
		fn_bnd_memory_free  ( sbulkparam2->sbndv, "total" );
	}
	else
	{
		printf("Error use 'fn_main_proj_bulk'.\n");
		printf("'main_type' should be 'proj_bulk1' or 'proj_bulk2'.\n");
	}

	/* release memory; */
	fn_bulk_memory_free	  ( sbulkparam1, "bulk param" );
	fn_bulk_memory_free	  ( sbulkparam2, "bulk param" );
	fn_GJP_memory_free	  ( ssysparam->sGJPv );
}


/**
 * \brief	Transform rhoJCplx (in the space of Fourier and GJPs) to 
 *			the common space (Fourier and GJPs);
 */
void fn_main_common_bulk ( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* true: save;	false: not save; */
	/* whether save densities; */
	sbulkparam1->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isSave	 = true; // 'fn_project_Fourier_GJP';
	sbulkparam1->srebndv->isSave = true; // 'fn_common_Fourier_GJP';
	sbulkparam2->srebndv->isSave = true; // 'fn_common_Fourier_GJP';
	/* variables for debug; */
	ssysparam->sGJPv->isTest	 = false; // 'fn_obt_system_gen_Jac_poly';
	sbulkparam1->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';
	sbulkparam2->sbndv->isTest	 = false; // 'fn_project_Fourier_GJP';
	sbulkparam1->srebndv->isTest = false; // 'fn_common_Fourier_GJP';
	sbulkparam2->srebndv->isTest = false; // 'fn_common_Fourier_GJP';

	
	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	/* obtain the stable bulk phases; */
	double hamilton1 = fn_obt_stable_bulk_phase ( sbulkparam1 );
	double hamilton2 = fn_obt_stable_bulk_phase ( sbulkparam2 );

	/* obtain the common 'rotateProjBoxMat'; */
	int status = fn_obt_commom_rotateProjBoxMat ( sbulkparam1, sbulkparam2, ssysparam );

	double hamilton = 0.0;
	if ( strcmp( main_type, "com_bulk1" ) == 0 )
	{
		/* prepare general Jacobi polynomials; */
		fn_obt_system_gen_Jac_poly ( ssysparam->sGJPv );

		/* project rhoCplx to the space of GJPs; */
		fn_project_Fourier_GJP ( sbulkparam1, ssysparam );

		printf(" <========== Transform to the common space ==========> \n\n");
		mytimer_t timer;
		timer.reset();
		timer.start();

		/* obtain the projection plane about the common 'rotateProjBoxMat'; */
		fn_get_system_projection_plane ( ssysparam );

		/**
		 * re-represent by the common 'rotateProjBoxMat';
		 *	including 'rhoJCplx', 'd0bnd', ...;
		 */
		fn_total_re_represent_common ( sbulkparam1, ssysparam );

		timer.pause();
		printf("\n\t ***** time cost of transformation to the common space: ");
		printf("%f seconds *****\n\n", timer.get_current_time());

		/* release memory; */
		fn_bulk_memory_free ( sbulkparam1, "project bulk" );
		fn_bnd_memory_free  ( sbulkparam1->sbndv, "total" );
		fn_bnd_memory_free  ( sbulkparam1->srebndv, "total" );
	}
	else if ( strcmp( main_type, "com_bulk2" ) == 0 )
	{
		/* prepare general Jacobi polynomials; */
		fn_obt_system_gen_Jac_poly ( ssysparam->sGJPv );

		/* project rhoCplx to the space of GJPs; */
		fn_project_Fourier_GJP ( sbulkparam2, ssysparam );

		printf(" <========== Transform to the common space ==========> \n\n");
		mytimer_t timer;
		timer.reset();
		timer.start();

		/* obtain the projection plane about the common 'rotateProjBoxMat'; */
		fn_get_system_projection_plane ( ssysparam );

		/**
		 * re-represent by the common 'rotateProjBoxMat';
		 *	including 'rhoJCplx', 'd0bnd', ...;
		 */
		fn_total_re_represent_common ( sbulkparam1, ssysparam );

		timer.pause();
		printf("\n\t ***** time cost of transformation to the common space: ");
		printf("%f seconds *****\n\n", timer.get_current_time());

		/* release memory; */
		fn_bulk_memory_free ( sbulkparam2, "project bulk" );
		fn_bnd_memory_free  ( sbulkparam2->sbndv, "total" );
		fn_bnd_memory_free  ( sbulkparam2->srebndv, "total" );
	}
	else
	{
		printf("Error use 'fn_main_proj_bulk'.\n");
		printf("'main_type' should be 'proj_bulk1' or 'proj_bulk2'.\n");
	}

	/* release memory; */
	fn_bulk_memory_free		( sbulkparam1, "bulk param" );
	fn_bulk_memory_free		( sbulkparam2, "bulk param" );
	fn_GJP_memory_free		( ssysparam->sGJPv );
	fn_tvec_free<int>		( ssysparam->NCpt );
	fn_tmat_free<double>	( ssysparam->rotateProjBoxMat );
	fn_tmat_free<int>		( ssysparam->sfftv->indKspace );
	fn_tmat_free<double>	( ssysparam->sfftv->projPlane );
}


/**
 * \brief	Calculation of the common projection matrix;
 */
void fn_main_common_rotateProjBoxMat ( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;
	
	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	/* obtain the common 'rotateProjBoxMat'; */
	int status = fn_obt_commom_rotateProjBoxMat ( sbulkparam1, sbulkparam2, ssysparam );

	/* release memory; */
	fn_bulk_memory_free		( sbulkparam1, "bulk param" );
	fn_bulk_memory_free		( sbulkparam2, "bulk param" );
	fn_tmat_free<double>	( ssysparam->rotateProjBoxMat );
}


/**
 * \brief	Display density according to the file about 'rhoJCplx';
 */
int fn_main_disp_system_density ( )
{
	/* introduce structure bodies; */
	stu_bulk_param sbulkparam1Type;
	stu_bulk_param *sbulkparam1 = &sbulkparam1Type;
	stu_bulk_param sbulkparam2Type;
	stu_bulk_param *sbulkparam2 = &sbulkparam2Type;
	stu_system_param ssysparamType;
	stu_system_param *ssysparam = &ssysparamType;

	/* input and save parameters; */
	fn_param_input ( sbulkparam1, sbulkparam2, ssysparam );

	FILE *fp;

	/* read 'd0JJ'; */
	char fname[FILELEN];
	sprintf(fname, "%s/d0GJP.dat", rsltDir); // file for saving parameters;
	fp = fopen(fname, "r");
	if ( fp == NULL )	// invalid file;
	{
		printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 0;
	}
	else
	{
		ssysparam->sGJPv->d0JJ = fn_tmat_read_double ( fname );
	}
	fclose(fp);

	/* read 'd0bnd'; */
	sprintf(fname, "%s/sys_d0bnd.dat", rsltDir); // file for saving parameters;
	fp = fopen(fname, "r");
	if ( fp == NULL )	// invalid file;
	{
		printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 0;
	}
	else	
	{
		ssysparam->scbndv->d0bnd = fn_tvec_read_complex ( fname );
	}
	fclose(fp);

	/* read 'x'; */
	sprintf(fname, "%s/x.dat", rsltDir); // file for saving parameters;
	fp = fopen(fname, "r");
	if ( fp == NULL )	// invalid file;
	{
		printf("Unable to read file '%s'. No such file or directory.\n", fname);
		return 0;
	}
	else	
	{
		ssysparam->sGJPv->x = fn_tvec_read_double ( fname );
	}
	fclose(fp);

	/* read the common 'rotateProjBoxMat'; */
	char paraFile[FILELEN];
	sprintf(paraFile, "%s/parameter_opt.dat", rsltDir);
	printf("parameter file: %s.\n\n", paraFile);
	fn_obt_com_projmat ( ssysparam, paraFile );	// read file;
	printf("\t ---> dimPhy = %d \t dimCpt = %d.\n", 
			ssysparam->dimRePhy, ssysparam->dimReCpt);
	printf("\t ---> The common projection matrix (read from %s) is \n", paraFile);
	fn_matrix_print ( ssysparam->rotateProjBoxMat );
	printf("\n");

	/* parameters; */
	ssysparam->scbndv->xlen		= ssysparam->sGJPv->d0JJ.row;
	ssysparam->scbndv->nd		= ssysparam->sGJPv->d0JJ.col;

	/* generate 'projPlane'; */
	fn_get_system_projection_plane ( ssysparam );
	ssysparam->scbndv->cplxDofs = ssysparam->cplxReDofs;

	/* read 'rhoJCplx'; */
	int	iter_max	= ssysparam->iter_max;
	int	save_step	= ssysparam->save_step;
	int	iterator	= 0;
	tvec<fftw_complex> rhoJCplx;
	iter_max = 0; // only "sys_rhoJCplx-1.dat";
	while ( iterator < iter_max )
	{
		if ( (iterator % save_step) == 0 )
		{
			sprintf(fname, "%s/sys_rhoJCplx%d.dat", rsltDir, iterator);
			fp = fopen(fname, "r");
			if ( fp == NULL )	// invalid file;
			{
				printf("Unable to read file '%s'. No such file or directory.\n", fname);
			}
			else
			{
				printf("Read file '%s'.\n", fname);
				rhoJCplx = fn_tvec_read_complex ( fname );
				fn_disp_system_density ( rhoJCplx, ssysparam, iterator );
				fn_tvec_free<fftw_complex>	( rhoJCplx );
			}
			fclose(fp);
		}
		iterator ++ ;
	}
	sprintf(fname, "%s/sys_rhoJCplx%d.dat", rsltDir, -1);
	printf("Read file '%s'.\n", fname);
	rhoJCplx = fn_tvec_read_complex ( fname );
	fn_disp_system_density ( rhoJCplx, ssysparam, -1 );

	/* memory free about the interface system; */
	fn_bulk_memory_free			( sbulkparam1, "bulk param" );
	fn_bulk_memory_free			( sbulkparam2, "bulk param" );
	fn_tmat_free<double>		( ssysparam->sGJPv->d0JJ );
	fn_tvec_free<fftw_complex>	( ssysparam->scbndv->d0bnd );
	fn_tvec_free<fftw_complex>	( rhoJCplx );
	fn_tvec_free<int>			( ssysparam->NCpt );
	fn_tvec_free<double>		( ssysparam->sGJPv->x );
	fn_tmat_free<double>		( ssysparam->rotateProjBoxMat );
	fn_tmat_free<double>		( ssysparam->sfftv->projPlane );
	fn_tmat_free<int>			( ssysparam->sfftv->indKspace );
	return 1;
}
