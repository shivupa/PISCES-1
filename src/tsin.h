#ifndef PISCES_TSIN_H_
#define PISCES_TSIN_H_
//
//  This is a small line-oriented input parser 
// 
//  the idea is to put input into groups and to use keyword = expressions
//  groups begin with "Begingroupname" and end with "Endgroupname"
//  where groupname is the string that needs to be passed to the 
//  TSIN member functions
//
//  it works a bit like the Fortran namelists, but in addition it can
//  return all the lines between BeginXXX and EndXXX
//
//  it is far from foolproof; especially the default values of
//  GetInt and GetDouble can lead to irritation, when the bloody
//  program does not set this variable right
//  there is usually a cookup with group or keynames
//
//  public member functions are:
//
//  int ReadFromFile    reads an input file and stores it line by line
//  int GetGroup        returns all the lines of a group
//  int GetInt          gets an integer from a specific group/key
//  double GetDouble    gets a double 
//  void GetIntArray    gets a set of integers
//  void GetDoubleArray gets a set of doubles
//
//
//  an example for an input file is
//           
//  BeginPara
//     kflag = 1,
//     freqs = 2, 3, 4.5,
//  EndPara
//
//  BeginBasis
//    whatever the basis is
//  EndBasis
//
//  after ReadFromFile has been called, you can write k = GetInt("Para", "kflag", default_for_k)
//  or get the basis set input by calling GetGroup("Basis", basis_input)
//

const int TSINMaxLines = 1000;   // this can be large; only required memory is allocated
const int TSINSigSigns = 15;     // significant signs of group names

class TSIN {
  
public :
    
  int ReadFromFile(const char *fname, int debug); 
  int GetGroup(const char *group, char ***inlines);
  int GetInt(const char *group, const char *key, int defall);
  long int GetLong(const char *group, const char *key, long int defall);
  double GetDouble(const char *group, const char *key, double defall);
  void GetString(const char *group, const char *key, const char *defall, char **rtn);
  void GetIntArray(const char *group, const char *key, int *fill, int n);
  void GetDoubleArray(const char *group, const char *key, double *fill, int n);
  
private :

  int nlines;
  char *line[TSINMaxLines];

};
#endif // PISCES_TSIN_H_
