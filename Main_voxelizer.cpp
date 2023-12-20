//basic type
#include "DGtal/base/Common.h"
#include "DGtal/base/BasicTypes.h"
//IO
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/readers/GenericReader.h"
//surface
#include "DGtal/topology/helpers/Surfaces.h"
//CLI
#include "CLI11.hpp"
//local
#include "Morphology3D.h"
//images/sets
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/ConstImageAdapter.h"
//frontier
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/SurfelSetPredicate.h"
//namespace
using namespace DGtal;
using namespace Z3i;


//binary image
typedef ImageContainerBySTLVector<Domain, bool> BinaryImage;
//image output
typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3D;
//adaptateur binary-> output
typedef DGtal::functors::Rescaling<bool,unsigned char> RescalFCT;
typedef ConstImageAdapter<BinaryImage, Image3D::Domain, functors::Identity, unsigned char, RescalFCT> ImageAdapt;

std::set<Z3i::SCell>
convertToSet(std::vector<Z3i::SCell> v){
    // Declaring the set
    std::set<Z3i::SCell> s;
    // Traverse the Vector
    for (Z3i::SCell x : v) {
        // Insert each element
        // into the Set
        s.insert(x);
    }
    // Return the resultant Set
    return s;
}

/**
méthode calculer la boite englobante du nuage de point.

lpcl  : un nuage de point

return : la boite englobante du nuage de point sous forme d'une paire de point.
**/
std::pair<RealPoint, RealPoint>
computeBoundingBox(std::vector<RealPoint> lpcl) {
    // Initialize bounding box to the first point.
    RealPoint minPoint = lpcl[0];
    RealPoint maxPoint = lpcl[0];
    // Update bounding box to include all points.
    for (const auto& point : lpcl) {
        minPoint = RealPoint(std::min(minPoint[0], point[0]),
                                  std::min(minPoint[1], point[1]),
                                  std::min(minPoint[2], point[2]));
        maxPoint = RealPoint(std::max(maxPoint[0], point[0]),
                                  std::max(maxPoint[1], point[1]),
                                  std::max(maxPoint[2], point[2]));
    }
    return std::make_pair(minPoint, maxPoint);
}

/**
méthode pour voxeliser un nuage de point. (floor)

lpcl  : un nuage de point
cSize :  la taille des cellules cubique en metre

return : une image binaire 3D dont les cellule true sont celles qui ont des points dedans.
**/
BinaryImage
voxelize(std::vector<RealPoint> lpcl,double cSize)
{
  std::pair<RealPoint, RealPoint> bbxyz=computeBoundingBox(lpcl);
  int xSize=floor((bbxyz.second[0] - bbxyz.first[0])/cSize);
  int ySize=floor((bbxyz.second[1] - bbxyz.first[1])/cSize);
  int zSize=floor((bbxyz.second[2] - bbxyz.first[2])/cSize);
  trace.info()<<"taille du dommaine : X(0,"<<xSize<<"),Y(0,"<<ySize<<"),Z(0,"<<zSize<<")"<<std::endl;
  //creation du domaine 3D
  BinaryImage::Domain PCLDomain(Point(0,0,0),Point(xSize,ySize,zSize));
  //init l'image 3D
  BinaryImage img(PCLDomain);
  //boucle sur les point 3D
  for(unsigned int i = 0; i < lpcl.size(); i++){
          trace.progressBar(i, lpcl.size());
          RealPoint mpCurrent=lpcl[i];
          //calcul des positions du point dans la discretisation
          unsigned int indx=floor((mpCurrent[0]- bbxyz.first[0])/cSize);
          unsigned int indy=floor((mpCurrent[1]- bbxyz.first[1])/cSize);
          unsigned int indz=floor((mpCurrent[2]- bbxyz.first[2])/cSize);
          img.setValue(Point(indx,indy,indz),1);
  }trace.info()<<""<<std::endl;
  return img;
}

BinaryImage
augmentDomain(const BinaryImage& img, int extent){
  Z3i::Point lowerBound = img.domain().lowerBound() - extent;
  Z3i::Point upperBound = img.domain().upperBound() + extent;
  BinaryImage largedImg(Domain(lowerBound, upperBound));
  for (auto &p : img.domain()) {
      largedImg.setValue(p, img(p));
  }
  return largedImg;
}

BinaryImage
reduceDomain(const BinaryImage& img, int extent) {
    Z3i::Point lowerBound = img.domain().lowerBound() + extent;
    Z3i::Point upperBound = img.domain().upperBound() - extent;
    BinaryImage reducedImg(Domain(lowerBound, upperBound));
    for (auto &p : reducedImg.domain()) {
        reducedImg.setValue(p, img(p));
    }
    return reducedImg;
}

BinaryImage
fillInterior(const BinaryImage& logVox){
  BinaryImage logVox_filled = BinaryImage(logVox.domain());
  //init khalinky
  KSpace ks;
  ks.init( logVox.domain().lowerBound(), logVox.domain().upperBound(), true );
  //image to set
  trace.info() << "transform to set..."<<  std::endl;
  DigitalSet set3d (logVox.domain());
  SetFromImage<DigitalSet>::append<BinaryImage>(set3d, logVox, false, true);
  trace.info() << "taille set : "<<set3d.size()<<  std::endl;
  //extract boundary
  trace.info() << "compute all boundary..."<<  std::endl;
  SurfelAdjacency<3> sAdj( false );
  ////méthode1 : on prend tout les contour on extrait le plus gros
  /*
  std::vector<std::vector<SCell> > vectConnectedSCell;
  Surfaces<KSpace>::extractAllConnectedSCell(vectConnectedSCell,ks, sAdj, set3d, false);
  trace.info() << "number of boundary : "<<vectConnectedSCell.size()<<  std::endl;
  //select the max contour
  trace.info() << "find highest contour (=consider like log boundary)... "<<  std::endl;
  int idMax=-1;
  int maxSizeContour=-1;
  for (int i =0;i < vectConnectedSCell.size();i++){
    int currentSize=vectConnectedSCell[i].size();
    trace.info() << currentSize<<  std::endl;
    if(currentSize > maxSizeContour){
      maxSizeContour=currentSize;
      idMax=i;
    }
  }
  std::vector<SCell> trunkBoundary=vectConnectedSCell[idMax];
  std::set<SCell> boundary = convertToSet(trunkBoundary);
  */
  //méthode2 : on utilise une cellule dont on sait quelle appartient au coutour (une Bel Pour Boundary element)
  std::set<SCell> boundary;
  //getting a bel
  KSpace::SCell bel = Surfaces<KSpace>::findABel( ks, logVox, (unsigned int)logVox.domain().size());
  //track closed boundary
  Surfaces<KSpace>::trackClosedBoundary(boundary,ks, sAdj, set3d, bel);
  //fill interior
  trace.info() << "fill interior..."<<  std::endl;
  unsigned int nbInt =Surfaces<KSpace>::uFillInterior(ks, functors::SurfelSetPredicate<std::set<SCell>, SCell>(boundary),logVox_filled, true);
  trace.info() << "Interior size: " << nbInt<<  std::endl;

  return logVox_filled;
}


int
main(int argc,char **argv){
  //nom du fichier input (.xyz)
  std::string inputFileName="../dataPCL/part_1.xyz";
  //nom du fichier output(.vol)
  std::string outputFileName="test";
  //taille des cellules de la discretisation en metre
  double cellsSize=0.01;//en mètre
  //rayon de fermeture
  double rayon=2.;
  //bollean to char
  RescalFCT rescale (0, 1, 0,255);
  functors::Identity id;
  //gestion des paramètres
  CLI::App app;
  app.description("Allowed options are: ");
  app.add_option("-i,--input", inputFileName , "input pcl file name.");
  app.add_option("-c,--cSize", cellsSize , "taille (en m) des cellules de la discretisation");
  app.add_option("-r,--rayonMorpho", rayon , "distance euclidienne pour le calcul de la fermeture");
  app.add_option("-o,--output", outputFileName , "output vol file name.");
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  //READ
  trace.beginBlock("read pcl...");
  std::vector<Z3i::RealPoint> log = PointListReader<Z3i::RealPoint>::getPointsFromFile(inputFileName);
  trace.endBlock();
  //VOXELIZE
  trace.beginBlock("Voxelize...");
  BinaryImage imgLog = voxelize(log,cellsSize);
  trace.endBlock();
  GenericWriter<ImageAdapt>::exportFile(outputFileName+".vol", ImageAdapt (imgLog,
                                                                      imgLog.domain(),
                                                                      id, rescale ));
  //FERMETURE
  trace.beginBlock("Fast Marching Methods (FMM) : fermeture...");
  Morphology3D::fermeture(imgLog,rayon);
  trace.endBlock();
  GenericWriter<ImageAdapt>::exportFile(outputFileName+"_ferme.vol", ImageAdapt (imgLog,
                                                                      imgLog.domain(),
                                                                      id, rescale ));
  //REMPLISSAGE
  trace.beginBlock("Filled interrior...");
  imgLog=fillInterior(imgLog);
  trace.endBlock();
  GenericWriter<ImageAdapt>::exportFile(outputFileName+"_filled.vol", ImageAdapt (imgLog,
                                                                      imgLog.domain(),
                                                                      id, rescale ));

}
