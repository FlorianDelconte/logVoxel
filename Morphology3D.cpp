#include "Morphology3D.h"
//surface
#include "DGtal/topology/helpers/Surfaces.h"
//frontiere
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/helpers/FrontierPredicate.h"
#include "DGtal/topology/LightExplicitDigitalSurface.h"
//image map
#include "DGtal/images/ImageContainerBySTLMap.h"
//fmm
#include "DGtal/geometry/volumes/distance/FMM.h"

BinaryImage
Morphology3D::augmentDomain(const BinaryImage& img, int extent){
  Z3i::Point lowerBound = img.domain().lowerBound() - extent;
  Z3i::Point upperBound = img.domain().upperBound() + extent;
  BinaryImage largedImg(Domain(lowerBound, upperBound));
  for (auto &p : img.domain()) {
      largedImg.setValue(p, img(p));
  }
  return largedImg;
}

BinaryImage
Morphology3D::reduceDomain(const BinaryImage& img, int extent) {
    Z3i::Point lowerBound = img.domain().lowerBound() + extent;
    Z3i::Point upperBound = img.domain().upperBound() - extent;
    BinaryImage reducedImg(Domain(lowerBound, upperBound));
    for (auto &p : reducedImg.domain()) {
        reducedImg.setValue(p, img(p));
    }
    return reducedImg;
}

void
Morphology3D::erosion(BinaryImage& img,double r){
  KSpace ks;
  ks.init( img.domain().lowerBound(), img.domain().upperBound(), true );
  KSpace::SCell bel;
  try {
    //getting a bel
    bel = Surfaces<KSpace>::findABel( ks, img, (unsigned int)img.domain().size() );
    trace.info() << "starting bel: " << bel <<" in domaine : "<<img.domain()<< std::endl;

  } catch (InputException& i) {
    trace.emphase() << "starting bel not found" << std::endl;
    exit(1);
  }
  trace.beginBlock("Implicit frontier...");
  //implicit frontier
  typedef functors::FrontierPredicate<KSpace, BinaryImage> SurfelPredicate;
  typedef LightExplicitDigitalSurface<KSpace, SurfelPredicate> Frontier;
  functors::SCellToIncidentPoints<KSpace> toIncidentPoints( ks );
  std::pair<Point,Point> bpair = toIncidentPoints( bel );
  SurfelPredicate surfelPredicate( ks, img, img( bpair.first ), img( bpair.second ) );
  Frontier frontier( ks, surfelPredicate, SurfelAdjacency<Z3i::KSpace::dimension>( true ), bel );
  trace.endBlock();
  //FMM
  DGtal::trace.beginBlock("FMM...");
  typedef ImageContainerBySTLMap<Domain,double> DistanceImage;
  typedef DigitalSetFromMap<DistanceImage> AcceptedPointSet;
  typedef FMM<DistanceImage, AcceptedPointSet, BinaryImage > FMM;
  DistanceImage imageDistance( img.domain(), 0.0 );
  AcceptedPointSet pointSet( imageDistance );
  functors::SCellToInnerPoint<KSpace> toInnerPoint( ks );
  for (Frontier::SurfelConstIterator it = frontier.begin(), itEnd = frontier.end();it != itEnd; ++it){
    pointSet.insert( toInnerPoint(*it) );
    imageDistance.setValue( toInnerPoint(*it), 0.5 );
  }
  FMM fmm( imageDistance, pointSet, img, img.domain().size(), r );
  fmm.compute();
  trace.info() << fmm << std::endl;
  trace.endBlock();
  //EROSION
  trace.beginBlock("Erosion");
  trace.info() << "Erosion of radius: " << r << std::endl;
  for (AcceptedPointSet::ConstIterator it = pointSet.begin(), itEnd = pointSet.end();it != itEnd; ++it){
    img.setValue(*it, 0);
  }
  trace.endBlock();
}

void
Morphology3D::dilatation(BinaryImage& img,double r){

  KSpace ks;
  ks.init( img.domain().lowerBound(), img.domain().upperBound(), true );
  KSpace::SCell bel;
  try {
    //getting a bel
    bel = Surfaces<KSpace>::findABel( ks, img, (unsigned int)img.domain().size() );
    trace.info() << "starting bel: " << bel <<" in domaine : "<<img.domain()<< std::endl;

  } catch (InputException& i) {
    trace.emphase() << "starting bel not found" << std::endl;
    exit(1);
  }

  trace.beginBlock("Implicit frontier...");
  typedef functors::FrontierPredicate<KSpace, BinaryImage> SurfelPredicate;
  typedef LightExplicitDigitalSurface<KSpace, SurfelPredicate> Frontier;
  functors::SCellToIncidentPoints<KSpace> toIncidentPoints( ks );
  std::pair<Point,Point> bpair = toIncidentPoints( bel );
  SurfelPredicate surfelPredicate( ks, img, img( bpair.first ), img( bpair.second ) );
  Frontier frontier( ks, surfelPredicate, SurfelAdjacency<KSpace::dimension>( true ), bel );
  trace.endBlock();

  trace.beginBlock("FMM...");
  typedef ImageContainerBySTLMap<Domain,double> DistanceImage;
  typedef DigitalSetFromMap<DistanceImage> AcceptedPointSet;
  typedef functors::NotPointPredicate<BinaryImage> BackgroundPredicate;
  typedef functors::BinaryPointPredicate<BackgroundPredicate,functors::DomainPredicate<Domain> > Predicate;
  BackgroundPredicate backgroundPredicate( img );
  functors::DomainPredicate<Domain> domainPredicate( img.domain());
  Predicate pred( backgroundPredicate, domainPredicate, functors::andBF2 );
  typedef FMM<DistanceImage, AcceptedPointSet, Predicate > FMM;
  DistanceImage imageDistance( img.domain(), 0.0 );
  AcceptedPointSet pointSet( imageDistance );
  functors::SCellToOuterPoint<KSpace> toOuterPoint( ks );
  for (Frontier::SurfelConstIterator it = frontier.begin(), itEnd = frontier.end();it != itEnd; ++it){
    pointSet.insert( toOuterPoint(*it) );
    imageDistance.setValue( toOuterPoint(*it), 0.5 );
  }
  FMM fmm( imageDistance, pointSet, pred, domainPredicate.domain().size(), r );
  fmm.compute();
  trace.info() << fmm << std::endl;
  trace.endBlock();
  //DILATATION
  trace.beginBlock("Dilatation");
  trace.info() << "Dilatation of radius: " << r << std::endl;
  for (AcceptedPointSet::ConstIterator it = pointSet.begin(), itEnd = pointSet.end();it != itEnd; ++it){
    img.setValue(*it, 1);
  }
  trace.endBlock();

}

void
Morphology3D::ouverture(BinaryImage& img,double r){
  Morphology3D::erosion(img,r);
  trace.info() << "_____________________________________________________________"<< std::endl;
  Morphology3D::dilatation(img,r);
}

void
Morphology3D::fermeture(BinaryImage& img,double r){
  //C'est ici que ca flanche, je maitrise pas tout ce code :/
  //J'me suis inspiré d'un tuto dgtal pour faire une erosion : https://dgtal.org/doc/stable/tutoFMMErosion.html
  //j'augmente la taille de l'image pour la dilatation et l'erosion.
  //en gros il faut grossir le domaine de l'image en fonction de r.
  img=augmentDomain( img, r+1);
  Morphology3D::dilatation(img,r);
  trace.info() << "_____________________________________________________________"<< std::endl;
  //pour l'erosion il faut ajouter 1 pcq on travaille sur le volume
  img=augmentDomain( img, 1);
  Morphology3D::erosion(img,r);
  //ici on retire les dimension
  img=reduceDomain( img, r+3);
}
