import os
import sys
from propagation import Propagator

WDIR = sys.argv[1]
FNIMG = sys.argv[2]
FNSEG = sys.argv[3]
NREF = int(sys.argv[4])
NSTART = int(sys.argv[5])
NSTOP = int(sys.argv[6])
TAG = sys.argv[7]
PATH_GREEDY = sys.argv[8]
PATH_VTKLEVELSET = sys.argv[9]

print("Propagating reference segmentation")

# Create a new Propagator
p = Propagator()

targetFrame = range(NSTART,NSTOP+1)

## Set parameters
p.SetTag(TAG)
p.SetInputImage(FNIMG)
p.SetReferenceSegmentation(FNSEG)
p.SetReferenceFrameNumber(NREF)
p.SetGreedyLocation(os.path.join(PATH_GREEDY,"greedy"))
p.SetVtkLevelSetLocation(os.path.join(PATH_VTKLEVELSET,"vtklevelset"))
p.SetTargetFrames(targetFrame)
p.SetOutputDir(os.path.join(WDIR, "output"))
p.SetSmoothingNumberOfIteration(5)
p.SetSmoothingPassband(0.05)

## Add additional mesh to warp
"""Reference mesh with empty string identifier is added by default and cannot be removed"""
#p.AddMeshToWarp('med', os.path.join(WDIR, 'segref.med.vtk'),False)
p.AddMeshToWarp('bnd', os.path.join(WDIR, 'segref.bnd.vtk'),False)

## check list of meshes to be warped
#p.GetWarpingList()
## remove meshe from the list
#p.RemoveMeshFromWarp('c')

### Optional Parameters for testing purpose
#p.SetFullResIterations('5x2x1')
#p.SetDilatedResIteration('5x2x1')
#p.SetGreedyThreads(6)

## Run propagation
"""
    - Set MeshWarpOnly to True to only warp meshes in the WarpingList based on existing 
      registration matrices. 
    - MeshWarpOnly mode is rely on previously generated registration matrices. Therefore, propagation
      with registration has to be completed before MeshWarpOnly run with same:
        - Reference frame
        - Output directory
        - Target frame
        - Tag
      Also do not move, delete, rename any files in the out/tmp folder.
      Or there will be missing file errors.
    - By default MeshWarpOnly is False, meaning propagator will run with full registration
      and mesh warping
"""
p.Run(MeshWarpOnly = False)
