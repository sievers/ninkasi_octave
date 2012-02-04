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
#include <ninkasi_projection.h>
#ifdef USE_HEALPIX
#include <chealpix.h>
#define PI_OVER_TWO 1.5707963267948966

#endif

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


octave_value pointer_to_value(void *ptr)
{
  int64NDArray myptr(1);
  long myptr_asint=(long)ptr;
  myptr(0)=myptr_asint;

  return myptr;
}



/*--------------------------------------------------------------------------------*/

DEFUN_DLD (clear_map, args, nargout, "Clear a skymap.\n")
{
  MAP  *map=(MAP *)get_pointer(args(0));
  clear_map(map);
  return octave_value_list();
}


/*--------------------------------------------------------------------------------*/

DEFUN_DLD (make_map_copy, args, nargout, "Copy a sky map.\n")
{
  MAP  *map=(MAP *)get_pointer(args(0));
  MAP *copy=make_map_copy(map);
  
  return pointer_to_value(copy);


}


/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_map_type_c, args, nargout, "Find out what kind of map we have.\n")
{
  MAP  *map=(MAP *)get_pointer(args(0));
  if (!map->projection)
    return octave_value("radec");

  switch(map->projection->proj_type){
  case (NK_RECT):
    return octave_value("radec");
    break;
  case (NK_CEA):
    return octave_value("cea");
    break;
  case (NK_TAN):
    return octave_value("tan");
    break;
  case (NK_HEALPIX_RING):
    return octave_value("ring");
    break;
  case (NK_HEALPIX_NEST):
    return octave_value("nest");
    break;
  default:
    printf("unkown map type.\n");
    return octave_value_list();
    break;
  }
  
  //never get here.
  return octave_value_list();

}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (destroy_map, args, nargout, "Free a skymap.\n")


{
  MAP  *map=(MAP *)get_pointer(args(0));
  destroy_map(map);
  return octave_value_list();
}


/*--------------------------------------------------------------------------------*/
DEFUN_DLD (octave2skymap_c, args, nargout, "Copy an octave array into a skymap.\n")
{
  if (args.length()!=2) {
    printf("Error in octave2skymap.  requires exactly two args.\n");
    return octave_value_list();
  }
    

  NDArray map=args(0).array_value();
  MAP *mymap=(MAP *)get_pointer(args(1));
  //printf("mapsizes are %d %d\n",mymap->nx,mymap->ny);
  double *mapvec=map.fortran_vec();
  long i;
  for (i=0;i<mymap->npix;i++) {
    mymap->map[i]=mapvec[i];
  }
  return octave_value_list();
}


/*--------------------------------------------------------------------------------*/
DEFUN_DLD (map_axpy, args, nargout, "Add a map to another map.\n")
{
  MAP  *y=(MAP *)get_pointer(args(0));
  MAP  *x=(MAP *)get_pointer(args(1));
  actData val=get_value(args(2));
  map_axpy(y,x,val);

  return octave_value_list();
}


/*--------------------------------------------------------------------------------*/
DEFUN_DLD (map_times_map, args, nargout, "Get the dot product of two maps.\n")
{
  MAP  *y=(MAP *)get_pointer(args(0));
  MAP *x;
  if (args.length()==1)
    x=y;
  else
    x=(MAP *)get_pointer(args(1));

  actData val=map_times_map(x,y);
  Matrix vv(1,1);
  vv(0,0)=val;

  return octave_value(vv);
}

/*--------------------------------------------------------------------------------*/
char *get_char_from_arg(charMatrix ch)
{
  int nn=ch.length();
  char *s=ch.fortran_vec();
  char *ss=strndup(s,nn+1);
  ss[nn]='\0';
  return ss;
}


/*--------------------------------------------------------------------------------*/

DEFUN_DLD (write_simple_map_c, args, nargout, "Write a simple map to disk.\n")
{
  MAP  *y=(MAP *)get_pointer(args(0));
  char *fname=get_char_from_arg(args(1).char_matrix_value());
    
  readwrite_simple_map(y,fname,DOWRITE);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_pix_from_radec_c, args, nargout, "Say where ra/dec should be in pixels in a map.")
{
  MAP  *map=(MAP *)get_pointer(args(0));
  actData ra=get_value(args(1));
  actData dec=get_value(args(2));

  int rapix,decpix;
  double x,y;
  octave_value_list retval;

  switch(map->projection->proj_type){
  case (NK_CEA):
    radec2pix_cea(map,ra,dec,&rapix,&decpix);
    retval(0)=rapix;
    retval(1)=decpix;
    break;
  case (NK_TAN):
    radec2xy_tan(&x,&y,ra,dec,map->projection);
    retval(0)=x;
    retval(1)=y;
    break;
#ifdef USE_HEALPIX
  case (NK_HEALPIX_RING):
    long ipix_ring;
    ang2pix_ring(map->projection->nside,PI_OVER_TWO-dec,ra, &ipix_ring);
    retval(0)=ipix_ring;
    retval(1)=0; //in case it's expected to have two outputs
    break;
  case (NK_HEALPIX_NEST):
    long ipix_nest;
    ang2pix_nest(map->projection->nside,PI_OVER_TWO-dec,ra, &ipix_nest);
    retval(0)=ipix_nest;
    retval(1)=0; //in case it's expected to have two outputs
    break;
#endif
    
  }


  return retval;

}


/*--------------------------------------------------------------------------------*/
/*
DEFUN_DLD (get_pix_from_radec_c, args, nargout, "Say what x/y coordinates an ra/dec location is in a map.")
{

  MAP  *map=(MAP *)get_pointer(args(0));
  double ra=get_value(args(1));
  double dec=get_value(args(2));
  double x,y;
  printf("calling radec2xy.\n");
  radec2xy_tan(&x,&y,ra,dec,map->projection);

  octave_value_list retval;
  retval(0)=x;
  retval(1)=y;

  return retval;

  
}
*/
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_radec_from_pix_c, args, nargout, "Say what ra/dec a pixel is in a map.")
{
  MAP  *map=(MAP *)get_pointer(args(0));
  int rapix=get_value(args(1));
  int decpix=get_value(args(2));
  actData ra,dec;
  switch(map->projection->proj_type) {
  case (NK_CEA):
    pix2radec_cea(map,rapix,decpix,&ra,&dec);
    break;
  case (NK_TAN):
    printf("doing tangent\n");
    pix2radec_tan(map,rapix,decpix,&ra,&dec);
    break;
#ifdef USE_HEALPIX
  case (NK_HEALPIX_RING):
    pix2ang_ring(map->projection->nside,rapix,&dec,&ra);
    dec=PI_OVER_TWO-dec;
    break;
  case (NK_HEALPIX_NEST):
    pix2ang_nest(map->projection->nside,rapix,&dec,&ra);
    dec=PI_OVER_TWO-dec;
    break;
#endif
  }
  octave_value_list retval;
  retval(0)=ra;
  retval(1)=dec;

  return retval;

}




/*--------------------------------------------------------------------------------*/
#ifdef USE_HEALPIX
DEFUN_DLD (set_skymap_healpix_ring_c, args, nargout, "Set map to be a healpix ring.  Args are (map, nside).\n")
{
  if (args.length()!=2) {
    printf("Need exactly two arguments to set_skymap_healpix_ring_c.\n");
    return octave_value_list();
  }
  MAP *map=(MAP *)get_pointer(args(0));
  int nside=(int)get_value(args(1));
  set_map_projection_healpix_ring(map,nside);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (set_skymap_healpix_nest_c, args, nargout, "Set map to be a healpix nest.  Args are (map, nside).\n")
{
  if (args.length()!=2) {
    printf("Need exactly two arguments to set_skymap_healpix_nest_c.\n");
    return octave_value_list();
  }
  MAP *map=(MAP *)get_pointer(args(0));
  int nside=(int)get_value(args(1));
  set_map_projection_healpix_nest(map,nside);
  return octave_value_list();
}
#endif //use_healpix

/*--------------------------------------------------------------------------------*/


DEFUN_DLD (set_skymap_cea_predef_c, args, nargout, "Dial in parameters for a CEA map.  Args are (map,radelt,decdelt,rapix,decpix,pv,nra,ndec\n")
{
  if (args.length()<8) {
    printf("Need 8 argumetns in set_skymap_cea_predef_c.\n");
    return octave_value_list();
  }
  MAP *map=(MAP *)get_pointer(args(0));
  actData radelt=get_value(args(1));
  actData decdelt=get_value(args(2));

  actData rapix=get_value(args(3));
  actData decpix=get_value(args(4));

  actData pv=get_value(args(5));
  int nra=get_value(args(6));
  int ndec=get_value(args(7));

  set_map_projection_cea_predef(map,radelt,decdelt,rapix,decpix,pv,nra,ndec);
  return octave_value_list();
  
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (set_skymap_tan_explicit_c, args, nargout, "Dial FITS parameters for a TAN map.  Args are (map,rapix,decpix,radelt,decdelt,pv,ra_cent,dec_cent,nra,ndec)\n")
{
  if (args.length()<10) {
    printf("Need 10 arguments in set_skymap_tan_explicit_c.\n");
    return octave_value_list();
  }
  MAP *map=(MAP *)get_pointer(args(0));
  actData rapix=get_value(args(1));
  actData decpix=get_value(args(2));
  actData radelt=get_value(args(3));
  actData decdelt=get_value(args(4));
  actData pv=get_value(args(5));
  actData ra_cent=get_value(args(6));
  actData dec_cent=get_value(args(7));
  int nra=(int)get_value(args(8));
  int ndec=(int)get_value(args(9));
  printf("setting myself to tan now.\n");
  set_map_projection_tan_explicit(map,rapix,decpix,radelt,decdelt,pv,ra_cent,dec_cent,nra,ndec);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (set_skymap_tan_predef_c, args, nargout, "Dial in parameters for a TAN map.  Args are (map,pixsize,rapix,decpix,ra_cent,dec_cent,nra,ndec)\n")
{
  if (args.length()<8) {
    printf("Need 8 argumetns in set_skymap_tan_predef_c.\n");
    return octave_value_list();
  }
  MAP *map=(MAP *)get_pointer(args(0));
  actData pixsize=get_value(args(1));

  actData rapix=get_value(args(2));
  actData decpix=get_value(args(3));

  actData ra_cent=get_value(args(4));
  actData dec_cent=get_value(args(5));

  int nra=get_value(args(6));
  int ndec=get_value(args(7));

  set_map_projection_tan_predef(map,ra_cent,dec_cent,rapix,decpix,pixsize,nra,ndec);
  printf("ra crap here is %14.4e %14.4e\n",map->projection->rapix,map->projection->ra_cent);
  return octave_value_list();
  
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (deres_map,args,nargout,"Reduce the resolution of a map by a factor of 2.")
{
  MAP *map=(MAP *)get_pointer(args(0));


  if (!map->projection) {
    printf("Map doesn't have a projection.\n");
    return octave_value_list();
  }
  MAP *map2=deres_map(map);

  int64NDArray myptr(1);
  myptr(0)=(long)map2;


  return octave_value(myptr);
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (upres_map,args,nargout,"Increase the resolution of a map by a factor of 2.")
{
  MAP *map=(MAP *)get_pointer(args(0));


  if (!map->projection) {
    printf("Map doesn't have a projection.\n");
    return octave_value_list();
  }
  MAP *map2=upres_map(map);

  int64NDArray myptr(1);
  myptr(0)=(long)map2;


  return octave_value(myptr);
}
/*--------------------------------------------------------------------------------*/


DEFUN_DLD (get_skymap_cea_params_c, args, nargout, "Return current cea params from a map.")
{
  MAP *map=(MAP *)get_pointer(args(0));
  if (!map->projection) {
    printf("Map doesn't have a projection.\n");
    return octave_value_list();
  }
  if (map->projection->proj_type!=NK_CEA) {
    printf("Map is not in CEA format.\n");
    return octave_value_list();
  }
  
  octave_value_list retval;
  retval(0)=map->projection->rapix;
  retval(1)=map->projection->decpix;
  retval(2)=map->projection->radelt;
  retval(3)=map->projection->decdelt;
  retval(4)=map->projection->pv;
  return retval;
  
}


/*--------------------------------------------------------------------------------*/


DEFUN_DLD (get_skymap_rect_params_c, args, nargout, "Return current rectangular params from a map.")
{
  MAP *map=(MAP *)get_pointer(args(0));
  
  octave_value_list retval;
  retval(0)=map->ramin;
  retval(1)=map->ramax;
  retval(2)=map->decmin;
  retval(3)=map->decmax;
  retval(4)=map->pixsize;
  return retval;
  
}


/*--------------------------------------------------------------------------------*/


DEFUN_DLD (get_skymap_tan_params_c, args, nargout, "Return current tan params from a map.")
{
  MAP *map=(MAP *)get_pointer(args(0));
  if (!map->projection) {
    printf("Map doesn't have a projection.\n");
    return octave_value_list();
  }
  if (map->projection->proj_type!=NK_TAN) {
    printf("Map is not in TAN format.\n");
    return octave_value_list();
  }
  
  octave_value_list retval;
  retval(0)=map->projection->rapix;
  retval(1)=map->projection->decpix;
  retval(2)=map->projection->radelt;
  retval(3)=map->projection->decdelt;
  retval(4)=map->projection->pv;
  retval(5)=map->projection->ra_cent;
  retval(6)=map->projection->dec_cent;
  return retval;
  
}

