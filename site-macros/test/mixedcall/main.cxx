#include <iostream>
#include "../../../../../safecast.hpp"
extern "C" void SUB();
extern "C" void SUB2FUN();


int main()
{
   SUB();
   SUB2FUN();
   return 0;
}
   
