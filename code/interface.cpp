#include "Head.h"
#include "Data.h"
#include "DataOperators.h"
#include "Mytimer.h"
#include "functs.h"
#include "umfpack.h"

int main(int argc, char* argv[])
{
	printf("\n\n ====================================   START PROGRAM  ======================================\n");
	pid_t pid = getpid();	// getpid: self; getppid: parent; getpgid: group;
	printf("\t\t program PID : %d\n\n", pid);

	mytimer_t timer;
	timer.reset();
	timer.start();

	/* the flag of parameters; */
	if ( argc == 1 )
	{
		strcpy(paraDir, "bcc_test");	// default folder for parameters;
		para_flag = "0";				// default parameters;
	}
	else if ( argc == 2 )
	{
		strcpy(paraDir, "bcc_test");	// default folder for parameters;
		para_flag = argv[1];
	}
	else
	{
		strcpy(paraDir, argv[1]);
		para_flag = argv[2];
	}

	/* obtain 'main_type'; */
	fn_obt_main_type ( );

	/* different functions; */
	double hamilton;
	if		( strcmp( main_type, "stable_system" ) == 0 ) // stable interface system;
	{
		hamilton = fn_main_stable_interface ( );
	}
	else if ( strcmp( main_type, "stable_bulk1" ) == 0 ||
			  strcmp( main_type, "stable_bulk2" ) == 0 ) // stable bulk phases;
	{
		hamilton = fn_main_stable_bulk ( );
	}
	else if ( strcmp( main_type, "proj_bulk1" ) == 0 ||
			  strcmp( main_type, "proj_bulk2" ) == 0 )	// project bulk phases;
	{
		fn_main_project_bulk ( );
	}
	else if ( strcmp( main_type, "com_bulk1" ) == 0 ||
			  strcmp( main_type, "com_bulk2" ) == 0 )  // bulk phases in the common space;
	{
		fn_main_common_bulk	 ( );
	}
	else if ( strcmp( main_type, "com_projmat" ) == 0 ) // the common projection matrix;
	{
		fn_main_common_rotateProjBoxMat ( );
	}
	else if ( strcmp( main_type, "disp_density" ) == 0 ) // display densities according to files;
	{
		int status = fn_main_disp_system_density ( );
	}
	else
	{
		printf("Error 'main_type': %s.\n", main_type);
	}


	timer.pause();
	printf("\n\n\n\t\t time cost of program : %f seconds\n", timer.get_current_time());
	printf("\n\n ======================================   END PROGRAM  ======================================\n\n");
	return 1;
}
