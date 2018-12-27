#include <iostream>

extern "C"
void FUN()
{
   std::cout << "called by fortran" << std::endl;
}

