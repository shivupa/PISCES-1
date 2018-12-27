#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include <mpi.h>
#include "tsin.h"



/*****************************************************
 *
 *  TSIN::ReadFromFile
 *
 *  read an input file 
 *
 *  lines beginning with # or empty lines are ignored
 *  all other lines are stored in the modul global char inpline
 *
 */
int TSIN::ReadFromFile(const char *fname, int debug)
{
  const int LineLength = 256;
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


  FILE   *flin;
  char   buffer[LineLength];
  int    i, l;
  int    nread = 0;

  /* check existence of input file */
  if ((flin = fopen(fname,"r")) == NULL) {
    if(rank==0)printf("Error in tsinreadfl opening %s\n", fname);
    exit(1);
  }
 
  /* read line by line allocating memory as required */
  nlines = 0;
  while (fgets(buffer, LineLength, flin)) {
    nread ++;
    l = strlen(buffer);
    if (buffer[0] == '#' || l == 1)
      continue;
    if (nlines > TSINMaxLines) {
     if(rank==0) printf("Error in rdinputfl: more than %d non-comment lines.\n", TSINMaxLines);
      exit(1);
    }
    line[nlines] = new char[l+1];
    strncpy(line[nlines], buffer, l+1);
    nlines ++;
  }
  fclose (flin);

  if (debug > 0)
    if(rank==0)printf("%d(%d) lines of %s have been read.\n", nread, nlines, fname);
  if (debug > 1) {
    for (i = 0; i < nlines; i++)
      if(rank==0)printf("%3d  %s", i, line[i]);
    if(rank==0)printf("EOF\n");
  }
  return(nlines);
}


/**********************************************
 *
 *  TSIN::GetGroup
 *
 *  search for group and return all the lines in inline
 *  memory for inlines is alloacted on the fly
 *
 *  return value is the number of lines copied to inlines
 */
int TSIN::GetGroup(const char *group, char ***inlines)
{
  char groupkey[TSINSigSigns+6];
  
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // look for Begingroup
  strcpy(groupkey,"Begin");
  strncat(groupkey,group,TSINSigSigns);
  int i1 = 0;
  while (i1 < TSIN::nlines) {
    if (strstr(TSIN::line[i1], groupkey) != NULL)
      break;
    i1 ++;
  }  
  if (i1 == TSIN::nlines) {
    if(rank==0)printf("Error in TSIN::GetGroup group:%s\n", group);
    if(rank==0) printf("Begin:group not found in input\n");
    exit(1);
  }

  // go to line after the group name and start seaching for Endgroup
  i1 ++;
  int i2 = i1;
  strcpy(groupkey,"End");
  strncat(groupkey,group,TSINSigSigns);
  while (i2 < TSIN::nlines) {
    if (strstr(TSIN::line[i2], groupkey) != NULL)
      break;
    i2 ++;
  }  
  if (i2 == TSIN::nlines) {
    if(rank==0)printf("Error in TSIN::GetGroup group:%s\n", group);
    if(rank==0)printf("End:group not found in input\n");
    exit(1);
  }

  int ngroup = i2 - i1;

  *inlines = new char* [ngroup];

  for (int k = 0; k < ngroup; ++k) {
    (*inlines)[k] = new char [strlen(line[i1+k])+1];
    strcpy((*inlines)[k], line[i1+k]);
  }

  return ngroup;

}



/**********************************************
 *
 *  TSIN::GetInt
 *
 *  search in group for key and return the integer if key exists
 *  if group or key do not exist retrun default
 *  if the key exists, scanf for key = 
 *
 */
int TSIN::GetInt(const char *group, const char *key, int defall)
{
  char *str = 0;
  char groupkey[TSINSigSigns+6];
  int i = 0;
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  
  // look for Begingroup
  strcpy(groupkey,"Begin");
  strncat(groupkey,group,TSINSigSigns);
  while (i < TSIN::nlines) {
    if ((str = strstr(TSIN::line[i], groupkey)) != NULL)
      break;
    i ++;
  }  
  if (i == TSIN::nlines)
    return(defall);  // EOF before group was found

  // each group name needs its line
  i ++;

  // look for key 
  strcpy(groupkey,"End");
  strncat(groupkey,group,TSINSigSigns);
  while (i < TSIN::nlines) {
    if ((str = strstr(TSIN::line[i], groupkey)) != NULL)  
      return(defall);  // end of group before key was found
    if ((str = strstr(TSIN::line[i], key)) != NULL)
      break;
    i ++;
  }
  if (i == TSIN::nlines)
    return(defall);  // EOF before key was found

  // now we have str pointing at the keyword
  int r;
  if (sscanf(str, "%*s = %i", &r) != 1) {
    if(rank==0)printf("Error in TSIN::GetInt  group:%s  key:%s\n", group, key);
    if(rank==0)printf("Could not get argument from line:%s\n", line[i]);
    exit(1);
  }
  return(r);
}


/**********************************************
 *
 *  TSIN::GetLong  same as GetINT
 *
 *  search in group for key and return the integer if key exists
 *  if group or key do not exist retrun default
 *  if the key exists, scanf for key = 
 *
 */
long int TSIN::GetLong(const char *group, const char *key, long int defall)
{
  char *str = 0;
  char groupkey[TSINSigSigns+6];
  int i = 0;
  
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // look for Begingroup
  strcpy(groupkey,"Begin");
  strncat(groupkey,group,TSINSigSigns);
  while (i < TSIN::nlines) {
    if ((str = strstr(TSIN::line[i], groupkey)) != NULL)
      break;
    i ++;
  }  
  if (i == TSIN::nlines)
    return(defall);  // EOF before group was found

  // each group name needs its line
  i ++;

  // look for key 
  strcpy(groupkey,"End");
  strncat(groupkey,group,TSINSigSigns);
  while (i < TSIN::nlines) {
    if ((str = strstr(TSIN::line[i], groupkey)) != NULL)  
      return(defall);  // end of group before key was found
    if ((str = strstr(TSIN::line[i], key)) != NULL)
      break;
    i ++;
  }
  if (i == TSIN::nlines)
    return(defall);  // EOF before key was found

  // now we have str pointing at the keyword
  long int r;
  if (sscanf(str, "%*s = %li", &r) != 1) {
    if(rank==0)printf("Error in TSIN::GetInt  group:%s  key:%s\n", group, key);
    if(rank==0)printf("Could not get argument from line:%s\n", line[i]);
    exit(1);
  }
  return(r);
}




/**********************************************
 *
 *  TSIN::GetDouble
 *
 *  search in group for key and return the integer if key exists
 *  if group or key do not exist retrun default
 *  if the key exists, scanf for key = 
 *
 */
double TSIN::GetDouble(const char *group, const char *key, double defall)
{
  char *str = 0;
  int i = 0;
  char groupkey[TSINSigSigns+6];
  
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // look for group
  strcpy(groupkey,"Begin");
  strncat(groupkey,group,TSINSigSigns);
  while (i < TSIN::nlines) {
    if ((str = strstr(TSIN::line[i], groupkey)) != NULL)
      break;
    i ++;
  }  
  if (i == TSIN::nlines)
    return(defall);  // EOF before group was found

  // each group name needs its line
  i ++;

  // look for key
  strcpy(groupkey,"End");
  strncat(groupkey,group,TSINSigSigns);
  while (i < TSIN::nlines) {
    if ((str = strstr(TSIN::line[i], groupkey)) != NULL)  
      return(defall);  // end of group before key was found
    if ((str = strstr(TSIN::line[i], key)) != NULL)
      break;
    i ++;
  }
  if (i == TSIN::nlines)
    return(defall);  // EOF before key was found

  // now we have str pointing at the keyword
  double r;
  if (sscanf(str, "%*s = %lf", &r) != 1) {
    if(rank==0)printf("Error in TSIN::GetDouble  group:%s  key:%s\n", group, key);
    if(rank==0)printf("Could not get argument from line:%s\n", line[i]);
    exit(1);
  }
  return(r);
}


/**********************************************
 *
 *  TSIN::GetString
 *
 *  search in group for key and return the string if key exists
 *  if group or key do not exist retrun default
 *  if the key exists, scanf for key =
 *
 */
void TSIN::GetString(const char *group, const char *key, const char *defall, char **rtn)
{
  char *str= 0;
  char groupkey[TSINSigSigns+6];
  int i = 0;

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // look for Begingroup
  strcpy(groupkey,"Begin");
  strncat(groupkey,group,TSINSigSigns);
  while (i < TSIN::nlines) {
    if ((str = strstr(TSIN::line[i], groupkey)) != NULL)
      break;
    i ++;
  }
  // EOF before group was found
  if (i == TSIN::nlines) {
    int n = strlen(defall) + 1;
    *rtn = new char[n];
    strcpy(*rtn, defall);
    return;
  }
  // each group name needs its line
  i ++;
  // look for key
  strcpy(groupkey,"End");
  strncat(groupkey,group,TSINSigSigns);
  while (i < TSIN::nlines) {
    if ((str = strstr(TSIN::line[i], groupkey)) != NULL) {
      int n = strlen(defall) + 1;
      *rtn = new char[n];
      strcpy(*rtn, defall);
      return;
    }
    if ((str = strstr(TSIN::line[i], key)) != NULL)
      break;
    i ++;
  }
  // EOF before key was found
  if (i == TSIN::nlines) {
    int n = strlen(defall) + 1;
    *rtn = new char[n];
    strcpy(*rtn, defall);
    return;
  }
  // now we have str pointing at the keyword
  char *s;
  int n = 0;
  s = new char[strlen(str)];
  while (str[0] != '\"')
    str ++;
  str ++;
  while (str[0] != '\"') {
    s[n] = str[0];
    n++;
    str ++;
  }
  s[n] = '\0';
  *rtn = new char[n];
  strcpy(*rtn, s);
  delete[] s;
  return;
}




/**********************************************
 *
 *  TSIN::GetIntArray
 *
 *  search in group for key 
 *  from this line try to read n integer into fill
 *
 *  if group or key do not exit with error
 *
 */
void TSIN::GetIntArray(const char *group, const char *key, int *fill, int n)
{
  char *str = 0;
  int i = 0;
  char groupkey[TSINSigSigns+6];
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


  // look for group
  strcpy(groupkey,"Begin");
  strncat(groupkey,group,TSINSigSigns);
  while (i < TSIN::nlines) {
    if ((str = strstr(TSIN::line[i], groupkey)) != NULL)
      break;
    i ++;
  }  
  if (i == TSIN::nlines) {
    if(rank==0)printf("Error in TSIN::GetIntArray group:%s, key:%s\n", group, key);
    if(rank==0)printf("group not found in input\n");
    exit(1);
  }


  // each group name needs its line
  i ++;

  // look for key 
  strcpy(groupkey,"End");
  strncat(groupkey,group,TSINSigSigns);
  while (i < TSIN::nlines) {
    if ((str = strstr(TSIN::line[i], groupkey)) != NULL) {
      if(rank==0)printf("Error in TSIN::GetIntArray group:%s, key:%s\n", group, key);
      if(rank==0)printf("End of group before key was found in input\n");
      exit(1);
    }
    if ((str = strstr(TSIN::line[i], key)) != NULL)
      break;
    i ++;
  }
  if (i == TSIN::nlines) {
    if(rank==0)printf("Error in TSIN::GetIntArray group:%s, key:%s\n", group, key);
    if(rank==0)printf("key not found in input\n");
    exit(1);
  }

  // now we have str pointing at the keyword
  // using scanf number by number
  // expected format is: key = n, l, k, ...
  // the last comma is vital

  // 1st argument
  if (sscanf(str, "%*s = %i,", &fill[0]) != 1) {
    if(rank==0)printf("Error in TSIN::GetIntArray  group:%s  key:%s\n", group, key);
    if(rank==0)printf("Could not get 1st argument from line:%s\n", line[i]);
    exit(1);
  }

  // (i+1)-th argument
  for (int k = 1; k < n; ++k) {
    // read till ","
    ++ str;
    while (*str != ',') {
      str ++;
      if (*str == '\0') {
	if(rank==0)printf("Error in TSIN::GetIntArray  group:%s  key:%s\n", group, key);
	if(rank==0)printf("Could not get argument %i from line:%s\n", k+1, line[i]);
	exit(1);
      }
    }
    if (sscanf(str, ", %i,", &fill[k]) != 1) {
    if(rank==0)printf("Error in TSIN::GetIntArray  group:%s  key:%s\n", group, key);
    if(rank==0)printf("Could not get argument %i from line:%s\n", k+1, line[i]);
    exit(1);
    }
  }

}

/**********************************************
 *
 *  TSIN::GetDoubleArray
 *
 *  search in group for key 
 *  from this line try to read n doubles into fill
 *
 *  if group or key do not exit with error
 *
 */
void TSIN::GetDoubleArray(const char *group, const char *key, double *fill, int n)
{
  char *str = 0;
  int i = 0;
  char groupkey[TSINSigSigns+6];
  
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // look for group
  strcpy(groupkey,"Begin");
  strncat(groupkey,group,TSINSigSigns);
  while (i < TSIN::nlines) {
    if ((str = strstr(TSIN::line[i], groupkey)) != NULL)
      break;
    i ++;
  }  
  if (i == TSIN::nlines) {
    if(rank==0)printf("Error in TSIN::GetDoubleArray group:%s, key:%s\n", group, key);
    if(rank==0)printf("group not found in input\n");
    exit(1);
  }

  // each group name needs its line
  i ++;

  // look for key 
  strcpy(groupkey,"End");
  strncat(groupkey,group,TSINSigSigns);
  while (i < TSIN::nlines) {
    if ((str = strstr(TSIN::line[i], groupkey)) != NULL) {
      if(rank==0)printf("Error in TSIN::GetDoubleArray group:%s, key:%s\n", group, key);
      if(rank==0)printf("End of group before key was found in input\n");
      exit(1);
    }
    if ((str = strstr(TSIN::line[i], key)) != NULL)
      break;
    i ++;
  }
  if (i == TSIN::nlines) {
    if(rank==0)printf("Error in TSIN::GetDoubleArray group:%s, key:%s\n", group, key);
    if(rank==0)printf("key not found in input\n");
    exit(1);
  }


  // now we have str pointing at the keyword
  // using scanf number by number
  // expected format is: key = n, l, k, ...
  // the last comma is vital

  // 1st argument
  if (sscanf(str, "%*s = %lf,", &fill[0]) != 1) {
    if(rank==0)printf("Error in TSIN::GetDoubleArray  group:%s  key:%s\n", group, key);
    if(rank==0)printf("Could not get 1st argument from line:%s\n", line[i]);
    exit(1);
  }

  // (i+1)-th argument
  for (int k = 1; k < n; ++k) {
    // read till ","
    ++ str;
    while (*str != ',') {
      str ++;
      if (*str == '\0') {
	if(rank==0)printf("Error in TSIN::GetDoubleArray  group:%s  key:%s\n", group, key);
	if(rank==0)printf("Could not get argument %i from line:%s\n", k+1, line[i]);
	exit(1);
      }
    }
    if (sscanf(str, ", %lf,", &fill[k]) != 1) {
    if(rank==0)printf("Error in TSIN::GetDoubleArray  group:%s  key:%s\n", group, key);
    if(rank==0)printf("Could not get argument %i from line:%s\n", k+1, line[i]);
    exit(1);
    }
  }

}

