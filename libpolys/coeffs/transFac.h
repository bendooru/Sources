#ifndef TRANSFAC_H
#define TRANSFAC_H

#include "factory/factory.h" // maybe less?
#include "polys/monomials/ring.h"
#include "polys/monomials/p_polys.h"

#include <list>

class transFac
{
  private:
    ring          f_r;         // polynomial ring where coefficient polynomials live
    CanonicalForm numerator;   // numerator of fraction
    CanonicalForm denominator; // factorized denominator of fraction

  public:
    transFac (const ring r);
    transFac (long n, const ring r);

    transFac (CanonicalForm n, const ring r);
    transFac (CanonicalForm n, CanonicalForm, const ring r);

    transFac (poly n, const ring r);
    transFac (poly n, poly d, const ring r);

    transFac (const transFac *l);

    ring getRing() const;
    CanonicalForm const& getNum()   const;
    CanonicalForm const& getDenom() const;

    poly getNumPoly() const;
    poly getDenomPoly() const;

    bool numeratorIsZero() const;
    bool numeratorIsOne() const;
    bool denominatorIsOne() const;

    void negateInplace();

    void normalize();
};

typedef transFac* pTransFac;

// needed at some points, as n_Init does not overload
number nfInit(poly, const coeffs);
int nftIsParam (number, const coeffs);
BOOLEAN n_transFacInitChar (coeffs, void*);

#endif
