#include "../../../src/C60.h"
#include <doctest/doctest.h>
//#include "../../../src/GetInput.h"
//#include "../../../src/Parameters.h"
#include "../../../src/tsin.h"
#include <fstream>
#include <iostream>
#include <streambuf>
#include <sstream>

TEST_CASE("C60 example regression test") {
  TSIN Input;
  Input.ReadFromFile("../../tests/regression_tests/c60/c60.inp", 5);
  /*
  Parameters InP;
  GetInputParameters(InP, Input);
*/
  std::ofstream file;
  file.open("c60_output.log");
  std::cout.rdbuf(file.rdbuf());
  C60_sp(Input);
  std::ifstream t("c60_output.log");
  std::string c60_output((std::istreambuf_iterator<char>(t)),
                         std::istreambuf_iterator<char>());
  size_t pos = c60_output.find("Eel"); 
  c60_output.erase(0, pos); 
  pos = c60_output.find("meV"); 
  c60_output.erase(c60_output.begin()+pos, c60_output.end()); 
  // String should now have in it
  // Eel = -0.00530703 Hartree = -144.412 
  pos = c60_output.find("="); 
  c60_output.erase(0, pos+1); 
  // -0.00530703 Hartree = -144.412 
  pos = c60_output.find("Hartree"); 
  size_t pos2 = c60_output.find("="); 
  c60_output.erase(c60_output.begin()+pos, c60_output.begin()+pos2+1); 
  // -0.00530703 -144.412 
  std::istringstream parser(c60_output);
  double har_E;
  double meV_E;
  parser >> har_E >> meV_E;
  CHECK(har_E == doctest::Approx(-0.00530703).epsilon(0.01));
  CHECK(meV_E == doctest::Approx(-144.412).epsilon(0.01));
}
