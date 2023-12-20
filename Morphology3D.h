#ifndef MORPHO3D_H
#define MORPHO3D_H
//basic type
#include "DGtal/base/Common.h"
#include "DGtal/base/BasicTypes.h"
//typedef
#include "DGtal/helpers/StdDefs.h"
//image
#include "DGtal/images/ImageContainerBySTLVector.h"

using namespace DGtal;
using namespace Z3i;

typedef ImageContainerBySTLVector<Domain, bool> BinaryImage;

class Morphology3D{
  public:
      Morphology3D(){}
      static void erosion(BinaryImage& img,double r);
      static void dilatation(BinaryImage& img,double r);
      static void ouverture(BinaryImage& img,double r);
      static void fermeture(BinaryImage& img,double r);
  private:
      static BinaryImage augmentDomain(const BinaryImage& img, int extent);
      static BinaryImage reduceDomain(const BinaryImage& img, int extent);
};

#endif
