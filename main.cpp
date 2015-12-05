#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "shock_data_types.hpp"	      /*tracer and shock data structers*/
#include "write_shock_catalogues.hpp" //correct catalogues
#include "blend_peaks.hpp"
#include "timer.h"

//#define SKIP_SEARCH
//#define SKIP_MERGE


//global variable definition

std::vector<tracer> tv;	/* tracers read in from file */
std::vector<tracer> td;	/* tracers sorted by density */
tracer tin;		          /* buffer for adding tracers to tracer vectors */

std::vector<shock>  s;   //shock catalogue
std::vector<shock>  bs;  //blended shock catalogue

std::vector<tracer> t;   //shock tracers
std::vector<tracer> bt;   //blended shock tracers

//main program

int main(int argc, char **argv)
{

  int n_peak_all;

  char	filebase[200];	//base filename (e.g., "turbulence")
  char	fdirbase[200];	//base directory name containing snaps (e.g., "find-dr")
  char	fdir_out[200];  //output directory
  char  fnamelist[200]; //file containing a list of peak catalogue directories
  char  fdircat[200];   //peak catalogue directory

  char  flist[200];     //peak catalogue list
  char  fdata[200];     //peak catalogue data
  int	  isnap;          //snapshot number
  int   ilist;
  int   nlist;

  FILE *fplist;         //file pointer for reading the list of peak catalogues

  long nt;              //total number of tracers in the peak catalogue


  //MPI wallclock timers
  double t_start_all;
  double t_end_all;
  double t_start;
  double t_end;


  //MPI rank, size, comm
  int rank;
  int np;

  //mpi world communicator
  MPI_Comm world = MPI_COMM_WORLD;




  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  ////	Begin program
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////



  //initialize MPI
  MPI_Init(&argc,&argv);  
  MPI_Comm_rank(world,&rank);
  MPI_Comm_size(world,&np);


  //start timer
  t_start_all = timer();


  //make directory
  sprintf(fdirbase,"%s",argv[1]);
  printf("Base directory: %s\n",fdirbase);

  //file containing the list of peak catalogue directories
  sprintf(fnamelist,"%s/%s",fdirbase,argv[2]);
  printf("List of peak catalogue directories: %s\n",fnamelist);

  //snapshot number
  isnap = atoi(argv[3]);
  printf("Snapshot number = %d.\n",isnap);

  //open the list of peak catalogue directories
  fplist = fopen(fnamelist,"r");
  fscanf(fplist,"%d\n",&nlist);

  //begin a loop over the peak catalogue directories
  //for(ilist=0;ilist<nlist;ilist++)
  for(ilist=0;ilist<15;ilist++)
  {
    printf("****************************\n");

    fscanf(fplist,"%s\n",fdircat);
    printf("ilist = %d, fdircat = %s\n",ilist,fdircat);

    //define the name of the peak list
    sprintf(flist,"%s/%s/peak.%04d.list",fdirbase,fdircat,isnap);
    printf("Reading catalogue list = %s.\n",flist);

    //define the name of the peak data
    sprintf(fdata,"%s/%s/peak.%04d.dat",fdirbase,fdircat,isnap);


    //read in the peak list
    read_shock_list(flist,&s);

    //print information about the shocks
    printf("Number of shocks = %ld\n",s.size());
    if(s.size()>0)
    {
      printf("Length of largest shock = %ld\n",s[0].l);

      //read in the peak data
      printf("Reading catalogue data = %s.\n",fdata);

      //resize the tracer array
      nt = 0;
      for(int si=0;si<s.size();si++)
        nt += s[si].l;
      t.resize(nt);

      //read the data into the tracer array
      read_shock_data(fdata,s,&t);
      printf("Total number of tracers = %ld (%ld).\n",t.size(),nt);


      //OK, now we blend the peaks
      blend_peaks(&bs, &bt, s, t);

    }

    //print information about the current peak list
    printf("Running number of blended peaks = %ld.\n",bs.size());

    //destroy the shock
    vector<shock>().swap(s);

    //destroy the tracers
    vector<tracer>().swap(t);
  }

  //close the list of peak catalogue directories
  fclose(fplist);

  //start timer
  t_end_all = timer();

  if(rank==0)
  {
    t_end = MPI_Wtime();
    printf("Total time = %fs.\n",t_end_all-t_start_all);
    printf("Done!\n");
    fflush(stdout);
  }

  MPI_Finalize();
  return 0;
 
}

