# VOXELISE AND FILL A 3D LOG
## dependencies
DGtal version 1.3 or later, see [DGtal installation] (https://github.com/DGtal-team/DGtal).
## compilation instructions
```
mkdir build
cd build
cmake ..  -DDGtal_DIR=/path/to/DGtal
make
```
## example
```
./voxelizer -c 0.01
```
-c : voxels size in meter
