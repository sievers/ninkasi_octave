%mex test_mex_mat.c CFLAGS=' -std=c99  -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -fopenmp -pthread' -lgomp CC=gcc4


cflags=' -std=c99  -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -fopenmp -pthread ';

includes=' -I /cita/d/raid-sievers/sievers/util/trinity/fftw-3.2/include';


%eval(['mex tod2map.c CFLAGS=''' cflags ''' -lgomp CC=gcc4'])
fftw=' -L/cita/d/raid-sievers/sievers/util/trinity/fftw-3.2/lib -lfftw3 '
crud=['mex  fft_omp.c   CFLAGS=''' cflags includes ''' ' fftw '  -lgomp CC=gcc4'];eval(crud)


mex skymap2tod.c CFLAGS=' -std=c99  -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -fopenmp -pthread' -lgomp CC=gcc4
%mex tod2map.c CFLAGS=' -std=c99  -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -fopenmp -pthread' -lgomp CC=gcc4


%mex fft_omp.c CFLAGS=' -std=c99  -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -fopenmp -pthread' -lgomp CC=gcc4


