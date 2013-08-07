#include <octave/oct.h>
#include <iostream>
#include <stdio.h>
#include <stdbool.h>
#ifdef __cplusplus
extern "C"
{
#endif
#include <ninkasi_config.h>
#include <ninkasi.h>
#include <dirfile.h>
#include <readtod.h>
#ifdef __cplusplus
}  /* end extern "C" */
#endif

using namespace std;

/*--------------------------------------------------------------------------------*/
void *get_pointer(octave_value val)
{
  int64NDArray myptr=val.array_value();
  long myptr2=myptr(0,0);
  return (void *)myptr2;

}
/*--------------------------------------------------------------------------------*/

actData get_value(octave_value val)
{
  NDArray myptr=val.array_value();
  actData myval=(actData)myptr(0,0);
  return myval;

}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (set_skymap_cea_simple_c, args, nargout, "Turn an RA/DEC map into a CEA map.\n")
{
  MAP *map=(MAP *)get_pointer(args(0));

  if (args.length()>1)
    set_map_projection_cea_simple_keeppix(map);
  else
    set_map_projection_cea_simple(map);
  return octave_value_list();
  
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (set_skymap_cea_simple_predef_c,args, nargout, "Turn an RA/DEC map into a CEA map with pixelsize/pv2_1 defined.\n")
{
  MAP *map=(MAP *)get_pointer(args(0));
  actData pixsize=get_value(args(1));
  actData pv=get_value(args(2));
  set_map_projection_cea_simple_predef(map,pixsize,pv);
  return octave_value_list();
  

}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (allocate_ninkasi_skymap, args, nargout, "Make a ninkasi skymap, return pointer to it.\n")
{
  
  octave_value_list retval;


  NDArray pixsize_in=args(0).array_value();
  NDArray ramin_in=args(1).array_value();
  NDArray ramax_in=args(2).array_value();
  NDArray decmin_in=args(3).array_value();
  NDArray decmax_in=args(4).array_value();
  actData pixsize=pixsize_in(0,0);
  actData ramin=ramin_in(0,0);
  actData ramax=ramax_in(0,0);
  actData decmin=decmin_in(0,0);
  actData decmax=decmax_in(0,0);

  mprintf(stdout,"pixsize is %10.4f, lims are %12.5f %12.5f %12.5f %12.5f\n",pixsize,ramin,ramax,decmin,decmax);

  //  MAP *map=(MAP *)malloc(sizeof(MAP));
  MAP *map=(MAP *)calloc(sizeof(MAP),1);
  map->have_locks=false;
  map->pixsize=pixsize;
  map->ramin=ramin;
  map->ramax=ramax;
  map->decmin=decmin;
  map->decmax=decmax;
#ifdef ACTPOL
  map->pol_state[0]=1;
#endif

  map->nx=(map->ramax-map->ramin)/map->pixsize+1;
  map->ny=(map->decmax-map->decmin)/map->pixsize+1;
  
  map->npix=map->nx*map->ny;
  map->map=(actData *)malloc(sizeof(actData)*map->npix);
  map->projection=(nkProjection *)malloc(sizeof(nkProjection));
  map->projection->proj_type=NK_RECT;
  map->have_locks=false;
  
  
  long mapptr=(long)map;
  //const dim_vector dims(3,2);
  //int64NDArray myptr(dims);
  int64NDArray myptr(1);
  
  //myptr(1)=map->nx;
  //myptr(2)=map->ny;
  myptr(0)=mapptr;
  //myptr(1)=1;
  //myptr(2)=2;
  //return octave_value(myptr);
  //cout << myptr.length() << "\n";

  retval(0)=myptr;
  return retval;

  //return octave_value(mapptr);

#if 0
  maps.nmap=1;
  {
    MAP *mapvec=(MAP *)malloc(maps.nmap*sizeof(MAP));
    for (int i=0;i<maps.nmap;i++)
      mapvec[i].have_locks=0;
    maps.maps=&mapvec;
  }

  maps.maps[0]->pixsize=params.pixsize;  //30 arcsecond pixels
	      
  maps.maps[0]->ramin=tods.ramin;
  maps.maps[0]->ramax=tods.ramax;
  maps.maps[0]->decmin=tods.decmin;
  maps.maps[0]->decmax=tods.decmax;
#endif

}

/*--------------------------------------------------------------------------------*/


DEFUN_DLD (print_map_type, args, nargout, "Print the kind of map we have inside a pointer.\n")
{
  MAP *mymap=(MAP *)get_pointer(args(0));
  if (!mymap) {
    printf("Map missing.\n");
    return octave_value_list();
  }
  if (!mymap->projection) {
    printf("Projection missing from map.\n");
    return octave_value_list();
  }
  switch(mymap->projection->proj_type) {
  case NK_RECT: {
    printf("Have rectangular map in RA/DEC.\n");
    break;
  }
  case NK_CEA: {
    printf("Have rectangular map in RA/DEC.\n");
    break;
  }
  default: {
    printf("Unkown type.\n");
    break;
  }
    
  }
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/


DEFUN_DLD (make_weightmap_ninkasi, args, nargout, "Make a weight map.\n")
{

  printf("hello!\n");
  int64NDArray mytods=args(0).array_value();
  int64NDArray mymap=args(1).array_value();
  int ntod=mytods.length();
  printf("Have %d tods.\n",ntod);
  int nmap=mymap.length();
  printf("Have %d maps.\n",nmap);



  long mapptr=mymap(0);
  MAP *map=(MAP *)mapptr;
  MAPvec maps;
  maps.nmap=1;
  maps.maps=&map;

  //memset(map->map,0,sizeof(actData)*map->npix);
  clear_mapset(&maps);
  printf("Cleared map is %d by %d\n",map->nx,map->ny);

  int i;
  for (i=0;i<ntod;i++) {
    long todptr=mytods(i);
    mbTOD *mytod=(mbTOD *)todptr;
    printf("tod nrow and ncol are %d %d\n",mytod->nrow,mytod->ncol);
    allocate_tod_storage(mytod);
    assign_tod_value(mytod,1.0);
    tod2map(map,mytod,NULL);
    free_tod_storage(mytod);    
  }
   
  printf("Made weight map.\n");
  

  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (skymap2octave, args, nargout, "Turn a ninkasi skymap into an octave one.\n")
{
  MAP *mymap=(MAP *)get_pointer(args(0));

  //dim_vector dims(mymap->nx,mymap->ny,get_npol_in_map(mymap));
  int npol=get_npol_in_map(mymap);
  if (npol==1) {
    dim_vector dims(mymap->nx,mymap->ny);
    NDArray map(dims);
    double *mapvec=map.fortran_vec();
    memcpy(mapvec,mymap->map,mymap->npix*sizeof(double));
    return octave_value(map);
  }
  else {
    dim_vector dims(npol,mymap->nx,mymap->ny);    
    NDArray map(dims);
    double *mapvec=map.fortran_vec();
    memcpy(mapvec,mymap->map,mymap->npix*sizeof(double)*npol);
    return octave_value(map);
  }
  return octave_value_list();  //never get here.

}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (mytest, args, nargout, "Make a skymap in a form ninkasi understands.\n")
{
  int64NDArray myval=args(0).array_value();
  long myptr=myval(0,0);
  MAP *mymap=(MAP *)myptr;
  printf("pixsize is %14.4e\n",mymap->pixsize);
  printf("Hello from second function!\n");
  return octave_value(0);
}
