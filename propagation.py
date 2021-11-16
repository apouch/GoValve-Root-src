import os
import vtk
import time
from Image4D import Image4D
from GreedyHelper import GreedyHelper
from shutil import copyfile

class Propagator:
    def __init__(self):
        # set default values
        self._fnimg = ""
        self._tag = "default"
        self._fref = -1
        self._outdir = "./out"
        self._greedyLocation = "greedy"
        self._vtklevelsetLocation = "vtklevelset"
        self._greedy = None
        self._fullResIteration = '100x30'
        self._dilatedResIteration = '100x100'
        self._metric_spec = 'SSD'
        self._threads = -1
        self._smoothingIter = 30
        self._smoothingPassband = 0.1
        self._meshWarpingList = {}
        self._isRegCompleted = False
        
    class MeshListItem:
        def __init__(self, filename, smooth):
            # file path of the mesh
            self._filename = filename
            # flat indicating if mesh to be smoothed after warping
            self._smooth = smooth

    def SetInputImage(self, _fnimg):
        self._fnimg = _fnimg

    def SetTag(self, _tag):
        self._tag = _tag

    def SetOutputDir(self, _outdir):
        self._outdir = _outdir

    def SetReferenceSegmentation(self, _fnsegref):
        self._fnsegref = _fnsegref
    
    def SetReferenceFrameNumber(self, _fref):
        self._fref = _fref
    
    def SetTargetFrames(self, _fnums):
        """Overrides existing target frames"""
        self._targetFrames = _fnums

    def SetTargetFrameRanges(self, _frange):
        """Overrides existing target frames"""
        self.__parseFrameRangeArray(_frange)

    def SetGreedyLocation(self, _greedyLoc):
        """
            Optional: Set the specific version of greedy for the propagation
            Default: run greedy from the path
        """
        self._greedyLocation = _greedyLoc

    def SetVtkLevelSetLocation(self, _vtklevelsetLoc):
        """
            Optional: Set the specific version of vtklevelset for the propagation
            Default: run vtklevelset from the path
        """
        self._vtklevelsetLocation = _vtklevelsetLoc

    def SetFullResIterations(self, _iter):
        """
            Optional: Set the Multi-Resolution Schedule (-n) parameter value for 
            full resolution greedy registration
            Default: 100x100
        """
        self._fullResIteration = _iter

    def SetDilatedResIteration(self, _iter):
        """
            Optional: Set the Multi-Resolution Schedule (-n) parameter value for
            dilated (forward & backward) greedy registration
            Default: 100x100
        """
        self._dilatedResIteration = _iter

    def SetMetricSpec(self, _metric_spec):
        """
            Optional: Set the Metric Specification (-m) parameter value for 
            greedy registration
            Default: SSD
        """
        self._metric_spec = _metric_spec

    def SetGreedyThreads(self, _threads):
        """ 
            Optional: Set the number of threads (-threads) greedy uses for registration
            Default: None (Unspecified)
        """
        self._threads = _threads

    def SetSmoothingNumberOfIteration(self, _iter):
        """
            Optional: Set number of iteration for taubin smoothing on mesh
            Default: 30
        """
        self._smoothingIter = _iter

    def SetSmoothingPassband(self, _passBand):
        """
            Optional: Set passband for taubin smoothing on mesh
            Type: double between 0 and 2
            Default: 0.1
        """
        self._smoothingPassband = _passBand

    def AddMeshToWarp(self, _id, _fnmesh, _smoothMesh):
        """
            Add additional mesh to be warped using the registration matrices
            _fnmesh: filename of the original mesh
            _id: Identifier to be added to the warped mesh
        """
        if _id == '':
            print('Empty string id is reserved for original reference mesh. Please use another id!')
            return
        
        if _id in self._meshWarpingList:
            print(f'Overriding existing mesh id = {_id} filename = {self._meshWarpingList[_id]._filename} with new mesh')

        self._meshWarpingList[_id] = self.MeshListItem(_fnmesh, _smoothMesh)


    def RemoveMeshFromWarp(self, _id):
        """
            Remove mesh with specified id from the warping list
        """
        if _id == '':
            print('Empty string id is reserved for original reference mesh and it cannot be removed!')
            return

        if _id not in self._meshWarpingList:
            print(f'Mesh with id {_id} does not exist in the warping list')
            return
        
        del self._meshWarpingList[_id]

    def GetWarpingList(self):
        """
            Return and print the warping list
        """
        print('----------------------------------')
        print('List of meshes to be warped')
        print('----------------------------------')
        for id in self._meshWarpingList:
            print(f'id: "{id}", filename: "{self._meshWarpingList[id]._filename}",\
                 smoothMesh: {self._meshWarpingList[id]._smooth}')

        return self._meshWarpingList


    def Run(self, MeshWarpOnly = False):
        self.__propagate(MeshWarpOnly)
        print('\r')
        if MeshWarpOnly:
            print('Mesh Warping Completed!')
        else:
            print('Propagation Completed!')


    def __parseFrameRangeArray(self, rangeArr):
        # todo: replace the placeholder
        self._targetFrames = []



    def __propagate(self, meshWarpingMode = False):
        #__propagate(fnimg, outdir = "", tag = "", seg_ref = "", self._targetFrames = [], fref = 0)

        """
        INPUTS:
        - fnimg: Image filename.
            - Supported formats:
                - (*.dcm) 4D Cartesian DICOM
                - (*.nii) 4D NIfTI
        - outdir: Name of existing directory for output files
        - tag: Study Identifier (e.g., 'bav08_root')
        - seg_ref: Filename of reference segmentation (nii)
        - self._targetFrames: List of frame numbers to propagate segementation to
        - fref: Reference segmentation frame

        Use full paths for all filenames.

        """

        """
        Process Inputs:
        - Load cartesian dicom image
        - Process reference segmentation
            - Dilate the reference segmentation
            - Create vtk mesh based on dialted reference segmentation
        - Extract and dialate 3D frame from the 4D image

        """
        # Validate files for the meshWarpingMode
        if meshWarpingMode:       
            print('==================================================')
            print('Propagation is running in Mesh Warping Mode!')
            print('==================================================')
            

        # Initialize Greedy Helper
        self._greedy = GreedyHelper(self._greedyLocation)

        # Validate Input Parameters
        if (self._fnimg == ""):
            raise RuntimeError("Input Image not set!")
        if (len(self._targetFrames) == 0):
            raise RuntimeError("Target Frames not set!")
        if (self._fref == -1):
            raise RuntimeError("Reference Frame Number not set!")
        if (self._fref not in self._targetFrames):
            # always including reference frame in the target frame
            self._targetFrames.append(self._fref)
            self._targetFrames.sort()
            print('Reference frame added to target frames')

        # create output directories
        meshdir = os.path.join(self._outdir, 'mesh')
        if not os.path.exists(self._outdir):
            # directory for output
            os.mkdir(self._outdir)
            # subdirectory for mesh outputs
            os.mkdir(meshdir)

        # recreate mesh dir in case it was removed manually
        if not os.path.exists(meshdir):
            os.mkdir(meshdir)

        tmpdir = os.path.join(self._outdir, "tmp")

        if not os.path.exists(tmpdir):
            os.mkdir(tmpdir)

        # performance logger
        perflog = {}
        timepoint = time.time()

        tag = self._tag
        fref = self._fref
        ref_ind = self._targetFrames.index(fref)


        if not meshWarpingMode:
            image = None
            
            # parse image type
            if self._fnimg.lower().endswith('.dcm'):
                print('Reading dicom image...')
                # Use the Dicom4D reader to create an Image4D object
                image = Image4D(self._fnimg, 'dicom')
                perflog['Dicom Loading'] = time.time() - timepoint
            elif self._fnimg.lower().endswith(('.nii.gz', '.nii')):
                print('Reading NIfTI image...')
                image = Image4D(self._fnimg, 'nifti')
                # Use the NIfTI reader to create an Image4D object
                perflog['NIfTI Loading'] = time.time() - timepoint
            else:
                print('Unknown image file type')
                return

            

            # Process reference segmentation
            # - Dilate the reference segmentation (mask)
            fn_mask_ref_srs = os.path.join(tmpdir, f'mask_{self._fref}_{self._tag}_srs.nii.gz')
            cmd = f'c3d -int 0 {self._fnsegref} -threshold 1 inf 1 0 \
                -dilate 1 10x10x10vox -resample 50% -o {fn_mask_ref_srs}'
            print("Dilating reference segmentation...")
            print(cmd)
            os.system(cmd)

            # - Create vtk mesh
            fn_mask_ref_vtk = os.path.join(tmpdir, f'mask_{self._fref}_{self._tag}.vtk')
            cmd = f'{self._vtklevelsetLocation} -pl {self._fnsegref} {fn_mask_ref_vtk} 1'
            print("Creating mask mesh...")
            print(cmd)
            os.system(cmd)

            # - add mesh to the meshWarpingList
            #   always smooth reference mesh
            self._meshWarpingList[''] = self.MeshListItem(fn_mask_ref_vtk, True)
        

            # - Create vtk mesh from the dilated mask
            fn_mask_ref_srs_vtk = os.path.join(tmpdir, f'mask_{self._fref}_{self._tag}_srs.vtk')
            cmd = f"{self._vtklevelsetLocation} -pl {fn_mask_ref_srs} {fn_mask_ref_srs_vtk} 1"
            print("Creating dilated mask mesh...")
            print(cmd)
            os.system(cmd)



            # Export 3D Frames
            # CartesianDicom.Export4D(os.path.join(outdir, 'img4D.nii.gz'))
            # Parallelizable

            
            timepoint = time.time()
            
            for i in self._targetFrames:
                fnImg = f'{tmpdir}/img_{i}_{tag}.nii.gz'
                fnImgRs = f'{tmpdir}/img_{i}_{tag}_srs.nii.gz'
                image.ExportFrame(i, fnImg)
                
                cmd = 'c3d ' + fnImg + ' -smooth 1mm -resample 50% \-o ' + fnImgRs
                print(cmd)
                os.system(cmd)
            
            perflog['Export 3D Frames'] = time.time() - timepoint
            
            
            # Initialize warp string array
            warp_str_array = [''] * len(self._targetFrames)
            
            # Propagate Forward
            print('---------------------------------')
            print('Propagating forward')
            print('---------------------------------')

            timepoint = time.time()

            for i in range(ref_ind, len(self._targetFrames) - 1):
                fCrnt = self._targetFrames[i]
                fNext = self._targetFrames[i + 1]
                fPrev = self._targetFrames[i - 1]

                print('Current Frame: ', self._targetFrames[i])

                if fCrnt == fref:
                    fn_mask_init = fn_mask_ref_srs
                else:
                    fn_mask_init = f'{tmpdir}/mask_{fPrev}_to_{fCrnt}_{tag}_srs_reslice_init.nii.gz'

                
                self.__propagation_helper(
                    work_dir = tmpdir,
                    crnt_ind = i,
                    is_forward = True,
                    mask_init = fn_mask_init,
                    warp_str_array = warp_str_array,
                    mask_ref_srs = fn_mask_ref_srs,
                    mask_ref_srs_vtk = fn_mask_ref_srs_vtk)
                
            perflog['Forward Propagation'] = time.time() - timepoint
            timepoint = time.time()

            # Propagate Backward
            print('---------------------------------')
            print('Propagating backward')
            print('---------------------------------')

            # - Clean up warp str array
            warp_str_array = [''] * len(self._targetFrames)

            for i in range(ref_ind, 0, -1):
                fCrnt = self._targetFrames[i]
                fNext = self._targetFrames[i - 1]
                fPrev = self._targetFrames[i + 1] if fCrnt != fref else -1

                print('Current Frame: ', self._targetFrames[i])

                if fCrnt == fref:
                    fn_mask_init = fn_mask_ref_srs
                else:
                    fn_mask_init = f'{tmpdir}/mask_{fPrev}_to_{fCrnt}_{tag}_srs_reslice_init.nii.gz'

                self.__propagation_helper(
                    work_dir = tmpdir,
                    crnt_ind = i,
                    is_forward = False,
                    mask_init = fn_mask_init,
                    warp_str_array = warp_str_array,
                    mask_ref_srs = fn_mask_ref_srs,
                    mask_ref_srs_vtk = fn_mask_ref_srs_vtk)
            
            perflog['Backward Propagation'] = time.time() - timepoint
            timepoint = time.time()

        # Propagate in Full Resolution
        if not meshWarpingMode:
            print('---------------------------------')
            print('Propagating in Full Resolution: ')
            print('---------------------------------\r')

        for i in range(0, len(self._targetFrames)):
            fCrnt = self._targetFrames[i]
            print('\r-----------------------------')
            print('Processing frame: ', fCrnt)
            print('-----------------------------\r')
            
            if fCrnt == fref:
                continue

            fn_img_fix = f'{tmpdir}/img_{fCrnt}_{tag}.nii.gz'
            fn_img_mov = f'{tmpdir}/img_{fref}_{tag}.nii.gz'

            # recursively build affine warp
            affine_warps = ''
            affine_warps_pts = ''

            if i > ref_ind:
                fPrev = self._targetFrames[i - 1]

                for j in range (ref_ind + 1, i + 1):
                    fn_regout_affine_init = f'{tmpdir}/affine_{self._targetFrames[j]}_to_{self._targetFrames[j - 1]}_srs_init.mat'
                    affine_warps = affine_warps + ' ' + fn_regout_affine_init + ',-1 '
                    affine_warps_pts = fn_regout_affine_init + ' ' + affine_warps_pts
            else:
                fPrev = self._targetFrames[i + 1]

                for j in range (ref_ind - 1, i - 1, -1):
                    fn_regout_affine_init = f'{tmpdir}/affine_{self._targetFrames[j]}_to_{self._targetFrames[j + 1]}_srs_init.mat'
                    affine_warps = affine_warps + ' ' + fn_regout_affine_init + ',-1 '
                    affine_warps_pts = fn_regout_affine_init + ' ' + affine_warps_pts

            #print("affine_warps: ", affine_warps)
            #print("affine_warps_pts: ", affine_warps_pts)

            # output file location
            fn_seg_reslice = f'{self._outdir}/seg_{fref}_to_{fCrnt}_{tag}_reslice.nii.gz'
            

            # transformation filenames
            fn_regout_deform = f'{tmpdir}/warp_{fref}_to_{fCrnt}.nii.gz'
            fn_regout_deform_inv = f'{tmpdir}/warp_{fref}_to_{fCrnt}_inv.nii.gz'

            if not meshWarpingMode:
                # full resolution mask for this frame
                mask_fix_srs = f'{tmpdir}/mask_{fPrev}_to_{fCrnt}_{tag}_srs_reslice_init.nii.gz'
                mask_fix = f'{tmpdir}/mask_{fPrev}_to_{fCrnt}_{tag}_reslice_init.nii.gz'
                print('Generating full res mask...')
                cmd = f'c3d -interpolation NearestNeighbor {fn_img_fix} {mask_fix_srs} -reslice-identity -o {mask_fix}'
                print(cmd)
                os.system(cmd)

                # trim to generate a reference frame
                fn_reference_frame = f'{tmpdir}/reference_{fPrev}_to_{fCrnt}_{tag}.nii.gz'
                cmd = f'c3d {mask_fix} -trim 0vox -o {fn_reference_frame}'
                print(cmd)
                os.system(cmd)

                # run registration and apply warp
                timepoint1 = time.time()

                print('Running Full Res Registration...')
                if self._fullResIteration != '100x100':
                    print(f'Using non-default parameters: -n {self._fullResIteration}')
                if self._metric_spec != 'SSD':
                    print(f'Using non-default parameter: -m {self._metric_spec}')

                self._greedy.run_reg(
                    img_fix = fn_img_fix,
                    img_mov = fn_img_mov,
                    affine_init = affine_warps,
                    regout_deform = fn_regout_deform,
                    regout_deform_inv = fn_regout_deform_inv,
                    reference_image = fn_reference_frame,
                    mask_fix = mask_fix,
                    multi_res_schedule = self._fullResIteration,
                    metric_spec = self._metric_spec,
                    threads = self._threads
                )

                perflog[f'Full Res Frame {fCrnt} - Registration'] = time.time() - timepoint1
                timepoint1 = time.time()

                print('Applying warp to segmentation...')
                self._greedy.apply_warp(
                    image_type = 'label',
                    img_fix = fn_img_fix,
                    img_mov = self._fnsegref,
                    img_reslice = fn_seg_reslice,
                    reg_affine = affine_warps,
                    reg_deform = fn_regout_deform
                )
                perflog[f'Full Res Frame {fCrnt} - Label Warp'] = time.time() - timepoint1

            timepoint1 = time.time()

            print('Applying warp to meshes...')
            mesh_warps = affine_warps_pts + fn_regout_deform_inv

            for id in self._meshWarpingList:
                meshItem = self._meshWarpingList[id]
                print(f'Mesh Warping {meshItem._filename}')

                fn_seg_reslice_vtk = os.path.join(meshdir, f'seg_{self._tag}_{id}_{fCrnt}.vtk')
                print(f'Output mesh filename: {fn_seg_reslice_vtk}')


                self._greedy.apply_warp(
                    image_type = 'mesh',
                    img_fix = fn_img_fix,
                    img_mov = meshItem._filename,
                    img_reslice = fn_seg_reslice_vtk,
                    reg_affine = mesh_warps
                )

                # Smooth warped mesh if smooth flag is True
                if (meshItem._smooth):
                    VTKHelper.SmoothMeshTaubin(fn_seg_reslice_vtk, fn_seg_reslice_vtk, self._smoothingIter, self._smoothingPassband)

                # Rename Point data (AMP)
                # VTKHelper.RenamePointData(fn_seg_reslice_vtk, 'Label')

            perflog[f'Full Res Frame {fCrnt} - Mesh Warp'] = time.time() - timepoint1

        # Processing reference mesh
        #   In the loop above only warped meshes (in target frames) are processed
        #   Following loop copy and smooth (if flagged) reference frame for each meshes in the list
        for id in self._meshWarpingList:
            # In mesh warping mode don't double smooth original reference mesh
            #   because original refernce mesh has already been smoothed during full mode run
            if meshWarpingMode and id == '':
                continue

            meshItem = self._meshWarpingList[id]

            # Target file name
            fn_seg_ref_vtk = os.path.join(meshdir, f'seg_{self._tag}_{id}_{self._fref}.vtk')

            # If mesh is flagged to be smoothed, smooth and export the mesh to the folder
            #   otherwise just copy the mesh without smoothing
            if (meshItem._smooth):
                VTKHelper.SmoothMeshTaubin(meshItem._filename, fn_seg_ref_vtk, self._smoothingIter, self._smoothingPassband)
            else:
                copyfile(meshItem._filename, fn_seg_ref_vtk)

        perflog['Full Res Propagation'] = time.time() - timepoint

        fn_perflog = os.path.join(self._outdir, 'perflog.txt')

        if os.path.exists(fn_perflog):
            fLog = open(fn_perflog, 'w')
        else:
            fLog = open(fn_perflog, 'x')

        for k in perflog:
            oneline = f'{k} : {perflog[k]}'
            print(oneline)
            fLog.write(oneline + '\n')

    def __propagation_helper(self, work_dir, crnt_ind, mask_init, \
        warp_str_array, mask_ref_srs, mask_ref_srs_vtk, is_forward = True):
        """
            Run propagation in one direction (forward or backward), 
            with downsampled images
        """

        fref = self._fref
        self._targetFrames = self._targetFrames
        tag = self._tag


        # forward propagation is in incremental order, backward is the reverse
        next_ind = crnt_ind + 1 if is_forward else crnt_ind - 1

        fCrnt = self._targetFrames[crnt_ind]
        fNext = self._targetFrames[next_ind]
        
        # current dilated image as fix
        fn_img_fix = f'{work_dir}/img_{fCrnt}_{tag}_srs.nii.gz'
        # next dilated image as moving image
        fn_img_mov = f'{work_dir}/img_{fNext}_{tag}_srs.nii.gz'
        
        # filenames of initial transformations
        fn_regout_affine = f'{work_dir}/affine_{fNext}_to_{fCrnt}_srs_init.mat'
        fn_regout_deform = f'{work_dir}/warp_{fNext}_to_{fCrnt}_srs_init.nii.gz'
        fn_regout_deform_inv = f'{work_dir}/warp_{fNext}_to_{fCrnt}_srs_init_inv.nii.gz'
        
        # call greedy to generate transformations
        self._greedy.run_reg(
            img_fix = fn_img_fix,
            img_mov = fn_img_mov,
            regout_affine = fn_regout_affine,
            regout_deform = fn_regout_deform,
            regout_deform_inv = fn_regout_deform_inv,
            mask_fix = mask_init,
            metric_spec=self._metric_spec,
            multi_res_schedule=self._dilatedResIteration,
            threads = self._threads
        )

        # Build warp string array recursively
        if fCrnt == fref:
            warp_str_array[crnt_ind] = f'{fn_regout_affine},-1 {fn_regout_deform_inv} '
        else:
            prev_ind = crnt_ind - 1 if is_forward else crnt_ind + 1
            warp_str_array[crnt_ind] = f'{warp_str_array[prev_ind]} {fn_regout_affine},-1 {fn_regout_deform_inv} '

        # Parallelizable
        fn_mask_init_reslice = f'{work_dir}/mask_{fCrnt}_to_{fNext}_{tag}_srs_reslice_init.nii.gz'
        fn_mask_init_reslice_vtk = f'{work_dir}/mask_{fCrnt}_to_{fNext}_{tag}_srs_reslice_init.vtk'

        # call greedy applying warp
        print('Applying warp to label...')
        self._greedy.apply_warp(
            image_type = 'label',
            img_fix = fn_img_mov,
            img_mov = mask_ref_srs,
            img_reslice = fn_mask_init_reslice,
            reg_affine = warp_str_array[crnt_ind]
        )

        print('Applying warp to mesh...')
        self._greedy.apply_warp(
            image_type = 'mesh',
            img_fix = fn_img_mov,
            img_mov = mask_ref_srs_vtk,
            img_reslice = fn_mask_init_reslice_vtk,
            reg_affine = warp_str_array[crnt_ind]
        )

class VTKHelper:
    @staticmethod
    def RenamePointData(fnMesh, newName):
        """A helper method that renames the point data array of a mesh"""
        _, x = os.path.splitext(fnMesh)
        assert(x == '.vtk')

        # Read vtk file
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(fnMesh)
        reader.Update()
        
        # Extract point data
        polyData = reader.GetOutput()
        pointData = polyData.GetPointData()

        # Extract the scalar array and set a new name
        scalarArr = pointData.GetScalars()
        scalarArr.SetName(newName)

        # Export and override existing file
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileName(fnMesh)
        writer.SetInputData(polyData)
        writer.Update()

    @staticmethod
    def SmoothMeshTaubin(fnMeshIn, fnMeshOut, numIter, passBand):
        """
            A helper method smoothes mesh using vtkWindowedSincPolyDataFilter
            From https://vtk.org/doc/nightly/html/classvtkWindowedSincPolyDataFilter.html#details
            ...
            " The total smoothing can be controlled by using two ivars. 

            The NumberOfIterations determines the maximum number of smoothing passes. The NumberOfIterations corresponds 
            to the degree of the polynomial that is used to approximate the windowed sinc function. Ten or twenty iterations 
            is all the is usually necessary. Contrast this with vtkSmoothPolyDataFilter which usually requires 100 to 200 
            smoothing iterations. vtkSmoothPolyDataFilter is also not an approximation to an ideal low-pass filter, which 
            can cause the geometry to shrink as the amount of smoothing increases.

            The second ivar is the specification of the PassBand for the windowed sinc filter. By design, the PassBand is 
            specified as a doubling point number between 0 and 2. Lower PassBand values produce more smoothing. A good 
            default value for the PassBand is 0.1 (for those interested, the PassBand (and frequencies) for PolyData are 
            based on the valence of the vertices, this limits all the frequency modes in a polyhedral mesh to between 0 and 2.)"
        """
        _, x = os.path.splitext(fnMeshIn)
        assert(x == '.vtk')
        _, x = os.path.splitext(fnMeshOut)
        assert(x == '.vtk')

        print(f'Mesh Smoothing: iter={numIter}, passBand={passBand}, Input:{fnMeshIn}, Output:{fnMeshOut}')

        # Read mesh input
        polyRdr = vtk.vtkPolyDataReader()
        polyRdr.SetFileName(fnMeshIn)
        polyRdr.Update()
        polyData = polyRdr.GetOutput()

        # Smooth
        sm = vtk.vtkWindowedSincPolyDataFilter()
        sm.SetInputData(polyData)
        sm.SetBoundarySmoothing(False)
        sm.SetNonManifoldSmoothing(False)
        sm.SetNormalizeCoordinates(True)
        sm.SetNumberOfIterations(numIter)
        sm.SetPassBand(passBand)
        sm.Update()

        # Write mesh output
        polyWtr = vtk.vtkPolyDataWriter()
        polyWtr.SetFileName(fnMeshOut)
        polyWtr.SetInputData(sm.GetOutput())
        polyWtr.Write()
