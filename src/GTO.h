//TS:
//
//  definition of the GTO class
//  this is just a plane old data structure
//
//  a GaussTypeOrbital is a Gaussian centered at the point x, y, z
//  it has an exponent a, an overall constant c, and the "angular momentum" lmn
//
//  f(X,Y,Z) = c * (X-x)^l * (Y-y)^m * (Z-z)^n * exp( -a * [(X-x)^2 + (Y-y)^2 + (Z-z)^2] )
//
//
//  as far as I can tell the identifier is set, but never used
//

class GTO
{
public:
   // parametrized constructor
   GTO(int id, double c, double a, int l, int m, int n, double x, double y, double z)
   : identifier(id), c(c), a(a), l(l), m(m), n(n), x(x), y(y), z(z)
   {}

public:
   int l, m, n;      // exponents of the (X-x)-like terms
   double c, a;      // overall coefficient and exponent
   double x, y, z;   // position of the GTO
   int identifier;   // connect GTO with the AtomSite
};

