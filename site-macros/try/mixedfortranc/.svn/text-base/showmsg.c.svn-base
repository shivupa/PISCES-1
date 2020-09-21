// C routine called from Fortran and prints a message
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void showmsg_(const char* c, int len) {
   char* msg = malloc(len+2);
   strncpy(msg,c,len);
   msg[len] = '\n';
   msg[len+1] = 0;
   printf(msg);
}
