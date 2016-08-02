  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <string.h>

  #include <gsl/gsl_integration.h>
  #include <gsl/gsl_sf.h>
  #include <iniparser.h>

  #include "const.h"
  #include "varb.h"
  #include "mymath.h"
  #include "myerr.h"
  #include "matrix.h"
  #include "init.h"
  #include "power.h"
  #include "cospara.h"
  #include "myinterpolate.h"

  #include "glbvarb.h"





  void import_data(char *fn, void *d, size_t size, size_t count) {
    FILE *fp;

    if(!(fp=fopen(fn, "r"))) {
      printf("can't open file `%s`\n", fn); fflush(stdout);
      exit(0);
      }

    if(!(fread(d, size, count, fp)) ) {
      printf("File '%s' import Error.\n", fn); fflush(stdout);
      exit(0);
      }
    
    fclose(fp);
    return;
    }


  void write_data(char *fn, void *d, size_t size, size_t count){
    FILE *fp;

    if(!(fp=fopen(fn, "wb"))) {
      printf("can't open file `%s`\n", fn); fflush(stdout);
      exit(0);
      }

    fwrite(d, size, count, fp);

    fclose(fp);
    return;
    }
