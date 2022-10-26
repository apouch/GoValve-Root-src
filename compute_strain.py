import numpy as np
import os
import sys
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import glob

def OrderList(flist,N,ref):
    '''
    Reorders the list so that the first file is the reference file, and returns new list of filenames and their IDs 
    Keyword arguments:
    flist -- list of filenames
    N -- number of files
    ref -- reference file name
    For example: for a list [file1.vtk, file2.vtk, file3.vtk, file4.vtk, file7.vtk, file8.vtk] and ref = file3.vtk
    it will return [file3.vtk file4.vtk file7.vtk file8.vtk file1.vtk file2.vtk], [3 4 7 8 1 2], and 3
    '''
    # Order filenames so that reference frame goes first
    Fno = np.zeros(N)
    FId = np.zeros(N)

    FListOrdered = [None]*N
    common = os.path.commonprefix(flist)
    for i, Fname in enumerate(flist):
        X = Fname.replace(common,'')
        X = X.replace('.vtk','')
        X = np.fromstring(X, dtype=int, sep=' ')
        #Get list of frame labels
        Fno[i] = X
        # Get label of reference frame
        if Fname==ref:
            refN = X
    # Sort fname labels
    Fno.sort()
    
    #Find Id of reference frame
    for i,X in enumerate(Fno):
        if X==refN:
            RefId = i

    # Sort reference area to be first in list
    FId[0:N-RefId] = Fno[RefId:N]
    FId[N-RefId:N]   = Fno[0:RefId]

    # Get list of file names in new order
    for i,F in enumerate(FId):
        for Fname in flist:
            X = Fname.replace(common,'')
            X = X.replace('.vtk','')
            X = np.fromstring(X, dtype=int, sep=' ')
            if X ==F:
                FListOrdered[i] = Fname

    return FListOrdered, FId, refN

def calDisp(polydata,Points,NP,RefPointsFixed,RefMids):
    '''
    Calculate displacement with respect to the reference points split into total, wall and root displacements
    Keyword arguments:
    polydata -- vtk object of the current timeframe
    Points -- numpy array of the current coordinates
    NP -- number of points
    RefPointsFixed -- numpy array of the reference coordinates, assumed to be fixed with geometric center at (0,0,0)
    Returns three numpy arrays:
    first -- total displacement
    second -- displacement of the wall relative to the fixed root
    third -- displacement of the whole root without the movement of the wall
    first == second + third
    '''
    #######################################
    # Calculate Displacements and Find Points Fixed Root
    TotDisp     = np.zeros((NP,3))
    WallDisp    = np.zeros((NP,3))
    RootDisp    = np.zeros((NP,3))
    PointsFixed = np.zeros((NP,3))

    # Define total displacement, relative to reference frame
    for i in range(NP):
        TotDisp[i,:] = polydata.GetPoint(i) - RefPointsFixed[i,:]

    # Find mid point of current frame
    Ranges = np.array(polydata.GetPoints().GetBounds())
    Mids = np.zeros((3))
    Mids[0] = np.mean(Ranges[0:2])
    Mids[1] = np.mean(Ranges[2:4])
    Mids[2] = np.mean(Ranges[4:6])


    for i in range(NP):
        # Define points with the centre fixed at the centre of the reference frame (RefMids)
        PointsFixed[i,:] = Points[i,:] - Mids + RefMids
        # Define wall displacement with a fixed root
        WallDisp[i,:]   = PointsFixed[i,:] - RefPointsFixed[i,:]
        # Define displacement of root, without wall movement
        RootDisp[i,:]   = TotDisp[i,:] - WallDisp[i,:]
    return TotDisp, WallDisp, RootDisp

def calAreaAndVol(polydata,ThickData,NC,NP):
    '''
    Calculate the wall area, wall volume, and the lumen volume
    Keyword arguments:
    polydata -- vtk object of the mesh
    ThickData -- numpy array of the thickness at the points
    NC -- number of cells
    NP -- number of points
    Returns three scalars: TotalWallArea, TotalWallVolume, TotalLumenVolume
    '''

    #######################################
    # Define Properties: Areas and Volumes

    TotalWallArea    = 0
    TotalWallVolume  = 0
    TotalLumenVolume = 0

    #Define Wall Area and Volume
    for i in range(NC):
        #Assign point id
        t = np.array([polydata.GetCell(i).GetPointId(j) for j in range(3)])
        #For each cell assign thickness for each point
        Thicks = [ThickData[t[j]] for j in range(3)]
        # Find and Sum Cell Areas and Volumes
        CA = polydata.GetCell(i).ComputeArea()
        CV = CA*(np.sum(Thicks))/3
        TotalWallArea   += CA
        TotalWallVolume += CV

    #Find Centre Points
    Ranges = np.array(polydata.GetPoints().GetBounds())
    Mids = np.zeros((3))
    Mids[0] = np.mean(Ranges[0:2])
    Mids[1] = np.mean(Ranges[2:4])
    Mids[2] = np.mean(Ranges[4:6])
    
    #Find edges of root mesh
    fedges = vtk.vtkFeatureEdges()
    fedges.BoundaryEdgesOn()
    fedges.FeatureEdgesOff()
    fedges.ManifoldEdgesOff()
    fedges.SetInputData(polydata)
    fedges.Update()
    ofedges = fedges.GetOutput()
    
    #Mark two edges separately
    Connect = vtk.vtkPolyDataConnectivityFilter()
    Connect.SetInputData(ofedges)
    Connect.SetExtractionModeToAllRegions()
    Connect.ColorRegionsOn()
    Connect.Update()
    Connect = Connect.GetOutput()

    #Choose specific edges
    Ring1 = vtk.vtkThreshold()
    Ring1.ThresholdByUpper(0.5)
    Ring1.SetInputData(Connect)
    Ring1.Update()
    Ring1 = Ring1.GetOutput()

    #Create cap
    Cap1 = vtk.vtkDelaunay2D()
    Cap1.SetProjectionPlaneMode(vtk.VTK_BEST_FITTING_PLANE)
    Cap1.SetInputData(Ring1)
    Cap1.Update()

    #Chose second ring
    Ring2 = vtk.vtkThreshold()
    Ring2.ThresholdByLower(0.5)
    Ring2.SetInputData(Connect)
    Ring2.Update()
    Ring2 = Ring2.GetOutput()

    #Chose second cap
    Cap2 = vtk.vtkDelaunay2D()
    Cap2.SetProjectionPlaneMode(vtk.VTK_BEST_FITTING_PLANE)
    Cap2.SetInputData(Ring2)
    Cap2.Update()

    #Find normals of cap 1 cells
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(Cap1.GetOutputPort())
    normals.ComputePointNormalsOff()
    normals.ComputeCellNormalsOn()
    normals.Update()
    Norm1 = normals.GetOutput()

    #Find centres of cap 1 cells
    vtkCenters=vtk.vtkCellCenters()
    vtkCenters.SetInputConnection(Cap1.GetOutputPort())
    vtkCenters.Update()

    #Get normals and points and their averages
    Norms = vtk_to_numpy(Norm1.GetCellData().GetNormals())
    Ps    = vtk_to_numpy(Norm1.GetPoints().GetData())
    NAvg1 = np.array([np.mean(Norms[:,i]) for i in range(3)])
    PAvg1 = np.array([np.mean(Ps[:,i]) for i in range(3)])

    #find direction of the root centre to the cap centre
    Cline1 = PAvg1-Mids
    u1 = NAvg1/np.linalg.norm(NAvg1)
    u2 = Cline1/np.linalg.norm(Cline1)
    angle = np.arccos(np.dot(u1,u2))

    #Check normals point inwards
    if angle<= (np.pi/2):
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputConnection(Cap1.GetOutputPort())
        normals.ComputePointNormalsOff()
        normals.ComputeCellNormalsOn()
        normals.FlipNormalsOn()
        normals.Update()
        Cap1 = normals.GetOutput()
    else:
        Cap1 = Cap1.GetOutput()

    #Repeat check for second cap
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(Cap2.GetOutputPort())
    normals.ComputePointNormalsOff()
    normals.ComputeCellNormalsOn()
    normals.Update()
    Norm2 = normals.GetOutput()

    vtkCenters=vtk.vtkCellCenters()
    vtkCenters.SetInputConnection(Cap2.GetOutputPort())
    vtkCenters.Update()

    Norms = vtk_to_numpy(Norm2.GetCellData().GetNormals())
    Ps    = vtk_to_numpy(Norm2.GetPoints().GetData())

    NAvg2 = np.array([np.mean(Norms[:,i]) for i in range(3)])
    PAvg2 = np.array([np.mean(Ps[:,i]) for i in range(3)])

    Cline2 = PAvg2 - Mids
    u1 = NAvg2/np.linalg.norm(NAvg2)
    u2 = Cline2/np.linalg.norm(Cline2)
    angle = np.arccos(np.dot(u1,u2))

    if angle<= (np.pi/2):
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputConnection(Cap2.GetOutputPort())
        normals.ComputePointNormalsOff()
        normals.ComputeCellNormalsOn()
        normals.FlipNormalsOn()
        normals.Update()
        Cap2 = normals.GetOutput()
    else:
        Cap2 = Cap2.GetOutput()

    polydata.ShallowCopy(polydata)
    Cap1.ShallowCopy(Cap1)
    Cap2.ShallowCopy(Cap2)

    #Attach Caps to root mesh
    appendFilter = vtk.vtkAppendPolyData()
    appendFilter.AddInputData(Cap1)
    appendFilter.AddInputData(polydata)
    appendFilter.AddInputData(Cap2)
    appendFilter.Update()

    cleanFilter = vtk.vtkCleanPolyData()
    cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
    cleanFilter.Update()
    
    #Get Lumen volume
    massUnion = vtk.vtkMassProperties()
    massUnion.SetInputConnection(cleanFilter.GetOutputPort())
    TotalLumenVolume = massUnion.GetVolume()

    return TotalWallArea, TotalWallVolume, TotalLumenVolume

def calStrains(polydata, RA, NC, l_Cell, c_Cell):
    '''
    Calculate strain invariants with respect to the reference basis
    Keyword arguments:
    polydata -- vtk object for the current time point
    RA -- numpy array of the reference basis vectors of size NC X 2 X 3
    NC -- number of cells
    l_Cell -- longitudinal unit vector 
    c_Cell -- circumferential unit vector 
    Returns two numpy arrays of invariant J and I1 at the cells, and the associated logngitudinal and circumferential strains
    Warning: l_Cell and c_Cell are mesh dependent, thus l_Strain and c_Strain are also mesh dependent
    '''
    I1 = np.zeros(NC)
    J  = np.zeros(NC)
    l_Strain = np.zeros(NC)
    c_Strain = np.zeros(NC)

    for i in range(NC):
        Cell = np.array(polydata.GetCell(i).GetPoints().GetData())

        # Define Reference Vectors
        A1 = RA[i,0,:]
        A2 = RA[i,1,:]

        #Define Deformed Vectors
        a1 = Cell[2,:]-Cell[0,:]
        a2 = Cell[2,:]-Cell[1,:]
        
        G      = np.array([[np.dot(A1,A1),np.dot(A1,A2)],[np.dot(A2,A1),np.dot(A2,A2)]])
        invG   = np.linalg.inv(G)
        
        A1dual = invG[0,0]*A1 + invG[0,1]*A2
        A2dual = invG[1,0]*A1 + invG[1,1]*A2
        F   = np.outer(a1,A1dual) + np.outer(a2,A2dual)
        C   = np.dot(F.T,F)
        C2  = np.dot(C,C)

        # Define principle strains and J
        trC   = C.trace()
        trC2  = C2.trace()
        I1[i] = trC
        J[i]  = np.sqrt((trC**2-trC2)/2)
        
        # Calculate Strains
        l_Strain[i] = np.dot(l_Cell[i,:],np.dot(C,l_Cell[i,:]))
        c_Strain[i] = np.dot(c_Cell[i,:],np.dot(C,c_Cell[i,:]))
        
    return J, I1, l_Strain, c_Strain


def calVectors(polydata,Pts,Cells,NP,NC,CellIds,STJMid,VAJMid):
    '''
    Finds the directional vectors for each point. 
    Keyword arguments:
    polydata -- vtk object for the current time point
    Pts -- vtk object for the current time point
    NP -- number of points
    NC -- number of cells
    CellIds -- list of point ids for each cell
    Returns the normal vectors, n, the longitudinal vectors, l, and the circumferential vectors, c, 
    for both the points and the cells (with suffices *_Pt and *_Cell respectively) 
    Warning: the longitudinal (and by asscociation circumferential) vectors are dependent on the mesh structure.
    '''
    
    #Get Mid-Line Vector
    MidLine = np.subtract(STJMid,VAJMid)
    MidLine /= np.linalg.norm(MidLine)

    #Use vtk to get point normals
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(polydata)
    normals.SetFlipNormals(True)
    normals.Update()
    
    n_Pt = vtk_to_numpy(normals.GetOutput().GetPointData().GetNormals())
    
    l_Pt = np.zeros((NP,3))
    c_Pt = np.zeros((NP,3))
    for i in range(NP):
        # Project Midline onto plane of normal
        l_Pt[i] = MidLine - np.multiply(np.dot(n_Pt[i],MidLine),n_Pt[i])/(np.dot(n_Pt[i],n_Pt[i]))
        # Find circumferential vector from normal and longitudinal
        c_Pt[i] = np.cross(l_Pt[i],n_Pt[i])
        
        # Normalise Vectors
        n_Pt[i] /= np.linalg.norm(n_Pt[i])
        l_Pt[i] /= np.linalg.norm(l_Pt[i])
        c_Pt[i] /= np.linalg.norm(c_Pt[i])
        
    n_Cell = np.zeros((NC,3)) 
    l_Cell = np.zeros((NC,3))
    c_Cell = np.zeros((NC,3))
    
    for i in range(NC):
        #Get normal to cell
        n_Cell[i] = np.cross(np.subtract(Cells[i,0],Cells[i,2]),np.subtract(Cells[i,0],Cells[i,1]))
        # Project midline onto plane of normal
        l_Cell[i] = MidLine - np.multiply(np.dot(n_Cell[i],MidLine),n_Cell[i])/(np.dot(n_Cell[i],n_Cell[i]))
        # Find circumferential vector from normal and longitudinal
        c_Cell[i] = np.cross(l_Cell[i],n_Cell[i])
        
        # Normalise Vectors
        n_Cell[i] /= np.linalg.norm(n_Cell[i])
        l_Cell[i] /= np.linalg.norm(l_Cell[i])
        c_Cell[i] /= np.linalg.norm(c_Cell[i])
    
    return n_Pt, l_Pt, c_Pt, n_Cell, l_Cell, c_Cell

def ProcessData(flist,ref,FT,OF=None,CF=None,prefix='Strains/',FixAndRotate=True,opformat='vtp'):
    '''
    Process all the files in a given list, including calculation of strains with respect to a reference
    Keyword arguments:
    flist -- list of filenames
    ref -- name of the reference file
    FT -- time between frames (or framerate)
    OF -- file number of the frame when AV opens
    CF -- file number of the frame when AV closes
    prefix -- prefix to the directory where the resulting files are stored (default 'Strains/')
    FixAndRotate -- whether the points should be fixed and rotated or not (default True)
    opformat -- format of the output files, either vtp (default) or vtk
    The function creates new .vtp files with the following data added:
    calculated strains, curvatures, cell areas and volumes, and field data of total wall area, total wall volume, total lumen volume, AV open/closed status
    
    Then the function returns the following numpy arrays (which have been reordered using OrderList):
    WallArea, WallVol, LumenVol, Time, Pts, WallAreaRatio, WallVolRatio, LumenVolRatio
    and N (the number of frames)
    '''
    print('##############################')
    print('Starting ProcessData Function')
    print('##############################')

    reader = vtk.vtkPolyDataReader()
    #################################
    # Define Reference Vectors
    #################################
    print('Reading Reference Frame:', ref)

    # Read the source file.
    reader.SetFileName(ref)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()

    polydata = reader.GetOutput()
    RefPoints = vtk_to_numpy(polydata.GetPoints().GetData())
    # Get Number of Points and Cells
    NP = polydata.GetNumberOfPoints()
    NC = polydata.GetNumberOfCells()
    if not polydata.GetPointData().GetArray(0):
        print('Warning Wall Thickness data is missing')
        ThickData = np.zeros((NP,3))
    else:
        ThickData = vtk_to_numpy(polydata.GetPointData().GetArray('Thickness'))
    #Default is set to fix the mesh to the reference position
    RefPointsFixed = RefPoints

    # Find mid point of current frame
    Ranges = np.array(polydata.GetPoints().GetBounds())
    RefMids = np.zeros((3))
    RefMids[0] = np.mean(Ranges[0:2])
    RefMids[1] = np.mean(Ranges[2:4])
    RefMids[2] = np.mean(Ranges[4:6])

    #Define initial frame of wall displacement for later use in motion calculation
    _, WallDispPF, _ = calDisp(polydata,RefPoints,NP,RefPointsFixed,RefMids)

    # Define empty array for reference basis
    RA = np.zeros((NC,2,3))
    Cells = np.zeros((NC,3,3))
    for i in range(NC):
        Cells[i,:,:] = vtk_to_numpy(polydata.GetCell(i).GetPoints().GetData())
        # Define refernce cell vectors
        RA[i,0,:] = Cells[i,2,:]-Cells[i,0,:]
        RA[i,1,:] = Cells[i,2,:]-Cells[i,1,:]


    ###########################################
    # Define Data At Every Time Frame
    ###########################################
    #Numer of frames
    N = len(flist)

    # Empty arrays for deformed frames

    Time          = np.zeros(N)
    ValvePosition = np.zeros(N)
    WallArea      = np.zeros(N)
    WallVol       = np.zeros(N)
    LumenVol      = np.zeros(N)
    WallAreaRatio = np.zeros(N)
    WallVolRatio  = np.zeros(N)
    LumenVolRatio = np.zeros(N)
    AvgJ          = np.zeros(N)
    AvgI1         = np.zeros(N)
    AvgJRatio     = np.zeros(N)
    AvgI1Ratio    = np.zeros(N)
    Pts           = np.zeros((N,NP,3))
    TotalMotion   = np.zeros((N,NP))
    
    # Re-order filename list to start with reference frame
    FListOrdered, FId, refN = OrderList(flist,N,ref)

    # Analyse each frame
    for X,Fname in enumerate(FListOrdered):
        if OF is not None and CF is not None:
            OpenX  = float(OF)-float(refN)
            CloseX = float(CF)-float(refN)
            if float(X)>=OpenX and float(X)<CloseX:
                ValvePosition[X] = 1

        print('Running Frame:', X+1)

        # Read the source file.
        reader.SetFileName(Fname)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()

        ######################################
        # Define File Data
        polydata  = reader.GetOutput()
        
        # Get time of frames
        if FId[0]<=FId[X]:
            Time[X] = float(FId[X])*float(FT) 
        else :
            if FId[X]-FId[X-1]<0:
                Time[X] = Time[X-1]+1*float(FT)
            elif FId[X]-FId[X-1]>0:
                Time[X] = Time[X-1]+(FId[X]-FId[X-1])*float(FT)

        if not polydata.GetPointData().GetArray(0): 
            print('Warning Wall Thickness data is missing')
            ThickData = np.zeros((NP,3))
        else:
            ThickData = vtk_to_numpy(polydata.GetPointData().GetArray('Thickness'))
        
        #Check that the number of points and cells are consistent
        if NP != polydata.GetNumberOfPoints():
            raise ValueError('Error: the number of points between the reference and deformed frames differ')
        if NC != polydata.GetNumberOfCells():
            raise ValueError('Error: the number of cells between the reference and deformed frames differ')

        # Save Point coordinates to array to be saved
        Pts[X,:,:] = vtk_to_numpy(polydata.GetPoints().GetData())

        Cells = np.zeros((NC,3,3))
        CellIds = np.zeros((NC,3))
        for i in range(NC):
            Cells[i,:,:] = vtk_to_numpy(polydata.GetCell(i).GetPoints().GetData())
            for j in range(3):
                CellIds[i,j] = (polydata.GetCell(i).GetPointIds().GetId(j))
        
        #########################################
        # Calculate Displacements
        TotDisp, WallDisp, RootDisp = calDisp(polydata,Pts[X,:,:],NP,RefPointsFixed,RefMids)
        
        #########################################
        # Calculate Total Wall Area and Volume, and Lumen Volume
        TotalWallArea, TotalWallVolume, TotalLumenVolume = calAreaAndVol(polydata,ThickData,NC,NP)

        # Save areas and volumes for each frame
        WallArea[X]       = TotalWallArea
        WallVol[X]        = TotalWallVolume
        LumenVol[X]       = TotalLumenVolume

        #########################################
        # Calculate Motion Vector (difference between current and previous frame of wall displacement)
        MotionVector = WallDisp - WallDispPF
        Motion = np.zeros(NP)

        #Define Motion magnitude
        for i in range(NP):
            Motion[i] = np.sqrt(sum(MotionVector[i,j]**2 for j in range (3)))
        
        #Define cumulative motion
        if X==0:
            TotalMotion[X,:] += Motion
        else :
            TotalMotion[X,:] = TotalMotion[X-1,:] + Motion

        #Update previous frame of wall disp for use in next frame
        WallDispPF = WallDisp
        
        #########################################
        #Mark Location of interatrial septum
        InterAtrialSeptum = np.zeros(NP)
        InterAtrialSeptum[0:25] = 1
        
        #Get STJ and VAJ Points
        STJ       = vtk_to_numpy(polydata.GetPointData().GetArray('STJ'))
        VAJ       = vtk_to_numpy(polydata.GetPointData().GetArray('VAJ'))
        
        STJPts = np.zeros((int(sum(STJ)),3))
        VAJPts = np.zeros((int(sum(VAJ)),3))
        k=0
        l=0
        for i in range(NP):
            if STJ[i] == 1:
                STJPts[k]=Pts[X,i,:]
                k+=1
            if VAJ[i] == 1:
                VAJPts[l]=Pts[X,i,:]
                l+=1
                
        STJMid = [sum(STJPts[:,0])/len(STJPts[:,0]),sum(STJPts[:,1])/len(STJPts[:,1]),sum(STJPts[:,2])/len(STJPts[:,2])]
        VAJMid = [sum(VAJPts[:,0])/len(VAJPts[:,0]),sum(VAJPts[:,1])/len(VAJPts[:,1]),sum(VAJPts[:,2])/len(VAJPts[:,2])]
        
        #########################################
        #Define Normal, Longitudinal, and Circumferential vectors
        n_Pt, l_Pt, c_Pt, n_Cell, l_Cell, c_Cell = calVectors(polydata,Pts[X,:,:],Cells,NP,NC,CellIds,STJMid,VAJMid)

        #########################################
        #Define Properties: J and I1
        J, I1, Long_Strain, Circ_Strain = calStrains(polydata,RA,NC, l_Cell, c_Cell)

        # Save areas and volumes for each frame
        AvgJ[X]       = np.mean(J)
        AvgI1[X]      = np.mean(I1)

        #########################################
        # Define Curvature

        # Define Gaussian Curvature
        curvgauss = vtk.vtkCurvatures()
        curvgauss.SetInputConnection(reader.GetOutputPort())
        curvgauss.SetCurvatureTypeToGaussian()
        curvgauss.Update()
        CurvG = curvgauss.GetOutput()
        CurvG.GetPointData().GetScalars()
        NumOfArrs = CurvG.GetPointData().GetNumberOfArrays()
        # Get curvature data that has been saved to new array
        CDG = np.asarray(CurvG.GetPointData().GetArray(NumOfArrs-1))  

        # Define Mean Curvature
        curvmean = vtk.vtkCurvatures()
        curvmean.SetInputConnection(reader.GetOutputPort())
        curvmean.SetCurvatureTypeToMean()
        curvmean.Update()
        CurvM = curvmean.GetOutput()
        CurvM.GetPointData().GetScalars()
        NumOfArrs = CurvM.GetPointData().GetNumberOfArrays()
        # Get curvature data that has been saved to new array
        CDM = np.asarray(CurvM.GetPointData().GetArray(NumOfArrs-1))

        #####################################
        # Add Data to Files and Write Files

        # Add Cell Data
        CellData  = [I1, J, Long_Strain, Circ_Strain]
        CellNames = ['I1_Cell','J_Cell','Long_Strain_Cell', 'Circ_Strain_Cell'] 
        for i in range(len(CellNames)) :
            arrayCell = vtk.util.numpy_support.numpy_to_vtk(CellData[i], deep=True)
            arrayCell.SetName(CellNames[i])
            dataCells = polydata.GetCellData()
            dataCells.AddArray(arrayCell)

        # Convert Cell Data to Point Data
        c2p = vtk.vtkCellDataToPointData()
        c2p.AddInputData(polydata)
        c2p.Update()
        c2p.GetOutput()
        NumOfArr = c2p.GetPolyDataOutput().GetPointData().GetNumberOfArrays()

        for i in range(NumOfArr):
            if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'I1_Cell':
                I1_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
            if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'J_Cell':
                J_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
            if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Long_Strain_Cell':
                l_Strain_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))
            if c2p.GetPolyDataOutput().GetPointData().GetArrayName(i) == 'Circ_Strain_Cell':
                c_Strain_Pt = vtk_to_numpy(c2p.GetPolyDataOutput().GetPointData().GetArray(i))

        # Add Point Data
        PointData = [CDG,CDM,I1_Pt,J_Pt,Motion,InterAtrialSeptum,TotalMotion[X],l_Strain_Pt,c_Strain_Pt]
        PointNames = ['Curv_Gaussian','Curv_Mean','I1_Pt','J_Pt','Motion','IAS','Total_Motion','Long_Strain_Pt','Circ_Strain_Pt']
        for i in range(len(PointNames)) :
            arrayPoint = vtk.util.numpy_support.numpy_to_vtk(PointData[i], deep=True)
            arrayPoint.SetName(PointNames[i])
            dataPoints = polydata.GetPointData()
            dataPoints.AddArray(arrayPoint)
            dataPoints.Modified()

        # Add Vector Data on Points
        VectorData = [TotDisp,WallDisp,RootDisp,n_Pt,l_Pt,c_Pt]
        VectorNames = ['Displacement_Total','Displacement_Wall','Displacement_Root','Normal_Pt','Longtudinal_Pt','Circumferential_Pt']
        for i in range(len(VectorNames)) :
            arrayVector = vtk.util.numpy_support.numpy_to_vtk(VectorData[i], deep=True)
            arrayVector.SetName(VectorNames[i])
            dataVectors = polydata.GetPointData()
            dataVectors.AddArray(arrayVector)
            dataVectors.Modified()

        # Add Vector Data on Cells
        VectorData = [n_Cell,l_Cell,c_Cell]
        VectorNames = ['Normal_Cell','Longtudinal_Cell','Circumferential_Cell']
        for i in range(len(VectorNames)) :
            arrayVector = vtk.util.numpy_support.numpy_to_vtk(VectorData[i], deep=True)
            arrayVector.SetName(VectorNames[i])
            dataVectors = polydata.GetCellData()
            dataVectors.AddArray(arrayVector)
            dataVectors.Modified()

        # Add Field Data
        if OF is not None and CF is not None:
            FieldData = [WallArea[X], WallVol[X],LumenVol[X],ValvePosition[X]]
            FieldNames = ['Wall_Area','Wall_Volume','Lumen_Volume','Valve_Position'] 
        else:
            FieldData = [WallArea[X], WallVol[X],LumenVol[X]]
            FieldNames = ['Wall_Area','Wall_Volume','Lumen_Volume'] 

        for i in range(len(FieldNames)) :
            arrayField = vtk.util.numpy_support.numpy_to_vtk(FieldData[i], deep=True)
            arrayField.SetName(FieldNames[i])
            dataFields = polydata.GetFieldData()
            dataFields.AddArray(arrayField)
            dataFields.Modified() 

        # Save points to be the reference frame with centre fixed at (0,0,0)
        if FixAndRotate == True:
            for i in range(NP):
                ptNew = RefPointsFixed[i,:]
                polydata.GetPoints().SetPoint(i, ptNew)

        #################################
        # Write data to vtp files
        if opformat == 'vtp':
            fname = os.path.join(prefix,os.path.split(os.path.splitext(Fname)[0])[1] + '.vtp')
        elif opformat == 'vtk':
            fname = os.path.join(prefix,os.path.split(os.path.splitext(Fname)[0])[1] + '.vtk')
        else:
            raise ValueError("Only vtp and vtk output formats are allowed")
        directory = os.path.dirname(fname)
        if not os.path.exists(directory):
            os.makedirs(directory)

        if opformat == 'vtp':
            writer = vtk.vtkXMLDataSetWriter()
        elif opformat == 'vtk':
            writer = vtk.vtkDataSetWriter()
        else:
            raise ValueError("Only vtp and vtk output formats are allowed")
        
        writer.SetFileName(fname)
        writer.SetInputData(polydata)
        writer.Write()

    WallAreaRatio[:]  = WallArea[:]/WallArea[0]
    WallVolRatio[:]   = WallVol[:]/WallVol[0]
    LumenVolRatio[:]  = LumenVol[:]/LumenVol[0]
    AvgJRatio[:]      = AvgJ[:]/AvgJ[0]
    AvgI1Ratio[:]     = AvgI1[:]/AvgI1[0]

    return WallArea, WallVol, LumenVol, Time, Pts, WallAreaRatio, WallVolRatio, LumenVolRatio, AvgJ, AvgI1, AvgJRatio, AvgI1Ratio, TotalMotion, N, FId

if __name__=='__main__':
    
    FixAndRotate = True
    
    WDIR = sys.argv[1]
    FT = float(sys.argv[2])
    OF = int(sys.argv[3])
    CF = int(sys.argv[4])
    refN = int(sys.argv[5])
        
    print("Computing root strain")
    print("Frame time:",FT)
    
    fnames = sorted(glob.glob(os.path.join(WDIR,'*med*.vtk')))
    fdir = os.path.dirname(fnames[0])

    common = os.path.commonprefix(fnames)
    for Fname in list(fnames):
        X = Fname.replace(common,'')
        X = X.replace('.vtk','')
        X = np.fromstring(X, dtype=int, sep=' ')
        X=X[0]
        if X==refN:
            ref=Fname
    NX = len(fnames)

    if not fnames:
        print(WDIR," is empty")
    else:
        fdir = os.path.dirname(fnames[0])
        # Check Directory
        if not os.path.exists(fdir):
            print('Error: Path does not exist:', fdir)
            sys.exit()
        WallArea, WallVol, LumenVol, Time, Pts, WallAreaRatio, WallVolRatio, LumenVolRatio, AvgJ, AvgI1, AvgJRatio, AvgI1Ratio, TotalMotion, N, FId = ProcessData(flist=fnames,ref=ref,FT=FT,OF=OF,CF=CF,prefix=os.path.join(WDIR,'Strains'),opformat='vtp')
   
    ###################################
    # Save data
    DataLocation = os.path.join(WDIR,'Data.npz')
    np.savez(DataLocation,Time=Time,Pts=Pts,WallArea=WallArea,WallVol=WallVol, LumenVol=LumenVol, WallAreaRatio=WallAreaRatio, WallVolRatio=WallVolRatio, LumenVolRatio=LumenVolRatio, AvgJ=AvgJ,AvgI1=AvgI1, AvgJRatio=AvgJRatio, AvgI1Ratio=AvgI1Ratio, N=N, OF=OF, CF=CF,refN = refN)

    CSVDataOriginal = np.zeros((NX,15))
    CSVDataStandardA = np.zeros((33,15))
    CSVDataStandardB = np.zeros((68,15))
    CSVDataStandard = np.zeros((101,15))

    CSVDataOriginal[:,0]  = Time
    CSVDataOriginal[:,1]  = WallArea
    CSVDataOriginal[:,2]  = WallVol
    CSVDataOriginal[:,3]  = LumenVol
    CSVDataOriginal[:,4]  = WallAreaRatio
    CSVDataOriginal[:,5]  = WallVolRatio
    CSVDataOriginal[:,6]  = LumenVolRatio
    CSVDataOriginal[:,7]  = AvgJ
    CSVDataOriginal[:,8]  = AvgI1
    CSVDataOriginal[:,9]  = AvgJRatio
    CSVDataOriginal[:,10] = AvgI1Ratio
    CSVDataOriginal[:,11] = np.full(NX,N)
    CSVDataOriginal[:,12] = np.full(NX,OF)
    CSVDataOriginal[:,13] = np.full(NX,CF)
    CSVDataOriginal[:,14] = np.full(NX,refN)

    TimeOpen  = np.linspace(0,320,33)
    TimeClose = np.linspace(330,1000,68)
    TimeStandard  = np.concatenate((TimeOpen,TimeClose))

    Oid = int(OF) - int(refN)
    Cid = int(CF) - int(refN)

    # CSVDataStandard[:,0] = np.concatenate((TimeStandard,TimeStandard[0]),axis=None)
    CSVDataStandard[:,0] = TimeStandard
    for i in range(1,15):
        CSVDataStandardA[:,i] = np.interp(np.linspace(Time[0],Time[Cid],33),Time[0:Cid+1],CSVDataOriginal[0:Cid+1,i])
        CSVDataStandardB[:,i] = np.interp(np.linspace(Time[Cid],Time[NX-1]+Time[0],68),np.concatenate((Time[Cid:NX],Time[NX-1]+Time[0]),axis=None),np.concatenate((CSVDataOriginal[Cid:NX,i],CSVDataOriginal[0,i]),axis=None))
        CSVDataStandard[:,i]  = np.concatenate((CSVDataStandardA[:,i],CSVDataStandardB[:,i]),axis=None)
    
    # With the Standardised data saved, reorder the original data to be in the file order
    # Find Id of reference frame
    for i,X in enumerate(FId):
        if X==np.min(FId):
            RefId = i

    # Sort reference area to be first in list
    FIdOriginal = np.zeros(N)
    FIdOriginal[0:N-RefId] = [int(i) for i in range(RefId,N)]
    FIdOriginal[N-RefId:N] = [int(i) for i in range(0,RefId)]
    
    for i,X in enumerate(FId[FIdOriginal.astype(int)]):
        if X==refN:
            refId=i
        if X==OF:
            OFId = i
        if X==CF:
            CFId = i
    
    CSVDataOriginal = np.zeros((NX+1,15))
    CSVDataOriginal[:,0]    = np.concatenate((FId[FIdOriginal.astype(int)]*FT,FId[FIdOriginal.astype(int)][NX-1]*FT+FId[FIdOriginal.astype(int)][0]*FT),axis=None)
    CSVDataOriginal[:NX,1]  = np.array(WallArea)[FIdOriginal.astype(int)]
    CSVDataOriginal[:NX,2]  = np.array(WallVol)[FIdOriginal.astype(int)]
    CSVDataOriginal[:NX,3]  = np.array(LumenVol)[FIdOriginal.astype(int)]
    CSVDataOriginal[:NX,4]  = np.array(WallAreaRatio)[FIdOriginal.astype(int)]
    CSVDataOriginal[:NX,5]  = np.array(WallVolRatio)[FIdOriginal.astype(int)]
    CSVDataOriginal[:NX,6]  = np.array(LumenVolRatio)[FIdOriginal.astype(int)]
    CSVDataOriginal[:NX,7]  = np.array(AvgJ)[FIdOriginal.astype(int)]
    CSVDataOriginal[:NX,8]  = np.array(AvgI1)[FIdOriginal.astype(int)]
    CSVDataOriginal[:NX,9]  = np.array(AvgJRatio)[FIdOriginal.astype(int)]
    CSVDataOriginal[:NX,10] = np.array(AvgI1Ratio)[FIdOriginal.astype(int)]
    CSVDataOriginal[:NX,11] = np.full(NX,N)
    CSVDataOriginal[:NX,12] = np.full(NX,OF)
    CSVDataOriginal[:NX,13] = np.full(NX,CF)
    CSVDataOriginal[:NX,14] = np.full(NX,refN)
    
    for i in range(1,15):
        CSVDataOriginal[:,i] = np.concatenate((CSVDataOriginal[:NX,i],CSVDataOriginal[0,i]),axis=None)
        
    CSVDataOriginalFile = os.path.join(WDIR,'CSVDataOriginal.csv')
    CSVDataStandardFile = os.path.join(WDIR,'CSVDataStandard.csv')

    np.savetxt(CSVDataOriginalFile, CSVDataOriginal, delimiter=",",header="Time, Wall Area,Wall Volume, Lumen Volume, Wall Area Ratio, Wall Volume Ratio, Lumen Volume Ratio, Avg J,Avg I1, Avg J Ratio, Avg I1 Ratio, No. of Frames, Open Frame, Closing Frame, Ref Frame")
    np.savetxt(CSVDataStandardFile, CSVDataStandard, delimiter=",",header="Time, Wall Area, Wall Volume, Lumen Volume, Wall Area Ratio, Wall Volume Ratio, Lumen Volume Ratio, Average J, Average I1, Average J Ratio, Average I1 Ratio, No. of Frames, Open Frame, Closing Frame, Ref Frame")

    Fig1name = 'NormalisedDataOriginalTime.png'
    Fig2name = 'NormalisedDataStandardisedTime.png'
    Fig3name = 'RawDataOriginalTime.png'
    Fig4name = 'RawDataStandardisedTime.png'

    # Figure 1   
    Plot_Ymax = np.amax(np.concatenate((CSVDataOriginal[:,4],CSVDataOriginal[:,5],CSVDataOriginal[:,6]),axis=None))
    Plot_Ymin = np.amin(np.concatenate((CSVDataOriginal[:,4],CSVDataOriginal[:,5],CSVDataOriginal[:,6]),axis=None))
    Plot_Yrange = Plot_Ymax - Plot_Ymin
    ymin = Plot_Ymin - Plot_Yrange*0.2
    ymax = Plot_Ymax + Plot_Yrange*0.1
    
    OpenFrameTime = np.ones(100)*CSVDataOriginal[OFId,0]
    CloseFrameTime = np.ones(100)*CSVDataOriginal[CFId,0]
    OpenFrameRange = np.linspace(ymin,np.max(np.concatenate((CSVDataOriginal[OFId,4],CSVDataOriginal[OFId,5],CSVDataOriginal[OFId,6]),axis=None)),100)
    CloseFrameRange = np.linspace(ymin,np.max(np.concatenate((CSVDataOriginal[CFId,4],CSVDataOriginal[CFId,5],CSVDataOriginal[CFId,6]),axis=None)),100)
    
    plt.figure(num=1,figsize=(7,4))
    plt.plot(OpenFrameTime,OpenFrameRange, '--',c='0.65',label='Open')
    plt.plot(CloseFrameTime,CloseFrameRange,  '--',c='0.35',label='Close')
    plt.plot(CSVDataOriginal[:,0], CSVDataOriginal[:,4],'r',label='Wall Area Ratio')
    plt.plot(CSVDataOriginal[:,0], CSVDataOriginal[:,5],'b',label='Wall Volume Ratio')
    plt.plot(CSVDataOriginal[:,0],CSVDataOriginal[:,6],'k',label='Lumen Volume Ratio')
    plt.set_title=('Raw Time Data')
    plt.xlabel('Time (ms)',size=20,fontname='Times')
    plt.ylabel('Ratio ',size=20,fontname='Times')
    plt.xticks(size=16,fontname='Times')
    plt.yticks(size=16,fontname='Times')
    lgd = plt.legend(bbox_to_anchor=(1.04,1), loc="upper left",prop={'family':'Times','size':'14'})
    plt.ylim(ymin, ymax)
    plt.savefig(Fig1name, format='png', bbox_inches='tight')


    # Figure 2
    # For standardised time, get index for time moment that is equivalent to the valve opening time
    if CSVDataOriginal[OFId,0]>CSVDataOriginal[refId,0]:
        OF_Stan = int(100*(CSVDataOriginal[OFId,0]-CSVDataOriginal[refId,0])/(np.max(CSVDataOriginal[:,0])-np.min(CSVDataOriginal[:,0])))
    elif CSVDataOriginal[OFId,0]<CSVDataOriginal[refId,0]:
        OF_Stan = int(100*(CSVDataOriginal[refId,0]-CSVDataOriginal[OFId,0])/(np.max(CSVDataOriginal[:,0])-np.min(CSVDataOriginal[:,0])))
    elif OFId == refId:
        OF_Stan = 0
        
    Plot_Ymax = np.amax(np.concatenate((CSVDataStandard[:,4],CSVDataStandard[:,5],CSVDataStandard[:,6]),axis=None))
    Plot_Ymin = np.amin(np.concatenate((CSVDataStandard[:,4],CSVDataStandard[:,5],CSVDataStandard[:,6]),axis=None))
    Plot_Yrange = Plot_Ymax - Plot_Ymin
    ymin = Plot_Ymin - Plot_Yrange*0.2
    ymax = Plot_Ymax + Plot_Yrange*0.1
    
    OpenFrameTime = np.ones(100)*CSVDataStandard[OF_Stan,0]
    CloseFrameTime = np.ones(100)*CSVDataStandard[33,0]
    OpenFrameRange = np.linspace(ymin,np.max(np.concatenate((CSVDataStandard[OF_Stan,4],CSVDataStandard[OF_Stan,5],CSVDataStandard[OF_Stan,6]),axis=None)),100)
    CloseFrameRange = np.linspace(ymin,np.max(np.concatenate((CSVDataStandard[33,4],CSVDataStandard[33,5],CSVDataStandard[33,6]),axis=None)),100)

    plt.figure(num=2,figsize=(7,4))
    plt.plot(OpenFrameTime,OpenFrameRange, '--',c='0.65',label='Open')
    plt.plot(CloseFrameTime,CloseFrameRange,  '--',c='0.35',label='Close')
    plt.plot(CSVDataStandard[:,0], CSVDataStandard[:,4],'r',label='Wall Area Ratio')
    plt.plot(CSVDataStandard[:,0],CSVDataStandard[:,5],'b',label='Wall Volume Ratio')
    plt.plot(CSVDataStandard[:,0],CSVDataStandard[:,6],'k',label='Lumen Volume Ratio')
    plt.set_title=('Standardised Time Data')
    plt.xlabel('Time (ms)',size=20,fontname='Times')
    plt.ylabel('Ratio ',size=20,fontname='Times')
    plt.xticks(size=16,fontname='Times')
    plt.yticks(size=16,fontname='Times')
    lgd = plt.legend(bbox_to_anchor=(1.04,1), loc="upper left",prop={'family':'Times','size':'14'})
    plt.ylim(ymin, ymax)
    plt.savefig(Fig2name, format='png', bbox_inches='tight')


    # Figure 3
    Plot_YmaxA = np.amax(CSVDataOriginal[:,1]/1000)
    Plot_YminA = np.amin(CSVDataOriginal[:,1]/1000)
    Plot_YrangeA = Plot_YmaxA - Plot_YminA
    yminA = Plot_YminA - Plot_YrangeA*0.2
    ymaxA = Plot_YmaxA + Plot_YrangeA*0.1
    
    Plot_YmaxB = np.amax(CSVDataOriginal[:,2]/1000)
    Plot_YminB = np.amin(CSVDataOriginal[:,2]/1000)
    Plot_YrangeB = Plot_YmaxB - Plot_YminB
    yminB = Plot_YminB - Plot_YrangeB*0.2
    ymaxB = Plot_YmaxB + Plot_YrangeB*0.1
    
    Plot_YmaxC = np.amax(CSVDataOriginal[:,3]/1000)
    Plot_YminC = np.amin(CSVDataOriginal[:,3]/1000)
    Plot_YrangeC = Plot_YmaxC - Plot_YminC
    yminC = Plot_YminC - Plot_YrangeC*0.2
    ymaxC = Plot_YmaxC + Plot_YrangeC*0.1
    
    OpenFrameTime    = np.ones(100)*CSVDataOriginal[OFId,0]
    CloseFrameTime   = np.ones(100)*CSVDataOriginal[CFId,0]
    OpenFrameRangeA  = np.linspace(yminA,np.max(CSVDataOriginal[OFId,1]/1000),100)
    CloseFrameRangeA = np.linspace(yminA,np.max(CSVDataOriginal[CFId,1]/1000),100)
    OpenFrameRangeB  = np.linspace(yminB,np.max(CSVDataOriginal[OFId,2]/1000),100)
    CloseFrameRangeB = np.linspace(yminB,np.max(CSVDataOriginal[CFId,2]/1000),100)
    OpenFrameRangeC  = np.linspace(yminC,np.max(CSVDataOriginal[OFId,3]/1000),100)
    CloseFrameRangeC = np.linspace(yminC,np.max(CSVDataOriginal[CFId,3]/1000),100)


    OpenFrameTime = np.ones(100)*CSVDataOriginal[OFId,0]
    CloseFrameTime = np.ones(100)*CSVDataOriginal[CFId,0]
    
    plt.figure(num=3,figsize=(25, 6))
    plt.subplot(1,3,1)
    plt.plot(OpenFrameTime,OpenFrameRangeA, '--',c='0.65')
    plt.plot(CloseFrameTime,CloseFrameRangeA,  '--',c='0.35')
    plt.plot(CSVDataOriginal[:,0], CSVDataOriginal[:,1]/1000,'r')
    plt.xlabel('Time (ms)',size=20,fontname='Times')
    plt.ylabel('Wall Area (x10$^3$ mm$^2$)',size=20,fontname='Times')
    plt.xticks(size=16,fontname='Times')
    plt.yticks(size=16,fontname='Times')
    plt.ylim(yminA, ymaxA)
    plt.subplot(1,3,2)
    plt.plot(OpenFrameTime,OpenFrameRangeB, '--',c='0.65')
    plt.plot(CloseFrameTime,CloseFrameRangeB,  '--',c='0.35')
    plt.plot( CSVDataOriginal[:,0],CSVDataOriginal[:,2]/1000,'b')
    plt.xlabel('Time (ms)',size=20,fontname='Times')
    plt.ylabel('Wall Volume (x10$^3$ mm$^3$)',size=20,fontname='Times')
    plt.xticks(size=16,fontname='Times')
    plt.yticks(size=16,fontname='Times')
    plt.ylim(yminB, ymaxB)
    plt.subplot(1,3,3)
    plt.plot(OpenFrameTime,OpenFrameRangeC, '--',c='0.65')
    plt.plot(CloseFrameTime,CloseFrameRangeC,  '--',c='0.35')
    plt.plot(CSVDataOriginal[:,0],CSVDataOriginal[:,3]/1000,'k')
    plt.set_title=('Raw Time Data')
    plt.xlabel('Time (ms)',size=20,fontname='Times')
    plt.ylabel('Lumen Volume (x10$^3$ mm$^3$)',size=20,fontname='Times')
    plt.xticks(size=16,fontname='Times')
    plt.yticks(size=16,fontname='Times')
    plt.ylim(yminC, ymaxC)
    plt.savefig(Fig3name, format='png', bbox_inches='tight')
    
    
    # Figure 4
    Plot_YmaxA = np.amax(CSVDataStandard[:,1]/1000)
    Plot_YminA = np.amin(CSVDataStandard[:,1]/1000)
    Plot_YrangeA = Plot_YmaxA - Plot_YminA
    yminA = Plot_YminA - Plot_YrangeA*0.2
    ymaxA = Plot_YmaxA + Plot_YrangeA*0.1
    
    Plot_YmaxB = np.amax(CSVDataStandard[:,2]/1000)
    Plot_YminB = np.amin(CSVDataStandard[:,2]/1000)
    Plot_YrangeB = Plot_YmaxB - Plot_YminB
    yminB = Plot_YminB - Plot_YrangeB*0.2
    ymaxB = Plot_YmaxB + Plot_YrangeB*0.1
    
    Plot_YmaxC = np.amax(CSVDataStandard[:,3]/1000)
    Plot_YminC = np.amin(CSVDataStandard[:,3]/1000)
    Plot_YrangeC = Plot_YmaxC - Plot_YminC
    yminC = Plot_YminC - Plot_YrangeC*0.2
    ymaxC = Plot_YmaxC + Plot_YrangeC*0.1
    
    OpenFrameTime = np.ones(100)*CSVDataStandard[OF_Stan,0]
    CloseFrameTime = np.ones(100)*CSVDataStandard[33,0]
    OpenFrameRangeA  = np.linspace(yminA,np.max(CSVDataStandard[OF_Stan,1]/1000),100)
    CloseFrameRangeA = np.linspace(yminA,np.max(CSVDataStandard[33,1]/1000),100)
    OpenFrameRangeB  = np.linspace(yminB,np.max(CSVDataStandard[OF_Stan,2]/1000),100)
    CloseFrameRangeB = np.linspace(yminB,np.max(CSVDataStandard[33,2]/1000),100)
    OpenFrameRangeC  = np.linspace(yminC,np.max(CSVDataStandard[OF_Stan,3]/1000),100)
    CloseFrameRangeC = np.linspace(yminC,np.max(CSVDataStandard[33,3]/1000),100)
    
    plt.figure(num=4,figsize=(25, 6))
    plt.subplot(1,3,1)
    plt.plot(OpenFrameTime,OpenFrameRangeA, '--',c='0.65')
    plt.plot(CloseFrameTime,CloseFrameRangeA,  '--',c='0.35')
    plt.plot(CSVDataStandard[:,0], CSVDataStandard[:,1]/1000,'r')
    plt.xlabel('Time (ms)',size=20,fontname='Times')
    plt.ylabel('Wall Area (x10$^3$ mm$^2$)',size=20,fontname='Times')
    plt.xticks(size=16,fontname='Times')
    plt.yticks(size=16,fontname='Times')
    plt.ylim(yminA, ymaxA)
    plt.subplot(1,3,2)
    plt.plot(OpenFrameTime,OpenFrameRangeB, '--',c='0.65')
    plt.plot(CloseFrameTime,CloseFrameRangeB,  '--',c='0.35')
    plt.plot( CSVDataStandard[:,0],CSVDataStandard[:,2]/1000,'b')
    plt.xlabel('Time (ms)',size=20,fontname='Times')
    plt.ylabel('Wall Volume (x10$^3$ mm$^3$)',size=20,fontname='Times')
    plt.xticks(size=16,fontname='Times')
    plt.yticks(size=16,fontname='Times')
    plt.ylim(yminB, ymaxB)
    plt.subplot(1,3,3)
    plt.plot(OpenFrameTime,OpenFrameRangeC, '--',c='0.65')
    plt.plot(CloseFrameTime,CloseFrameRangeC,  '--',c='0.35')
    plt.plot(CSVDataStandard[:,0],CSVDataStandard[:,3]/1000,'k')
    plt.set_title=('Raw Time Data')
    plt.xlabel('Time (ms)',size=20,fontname='Times')
    plt.ylabel('Lumen Volume (x10$^3$ mm$^3$)',size=20,fontname='Times')
    plt.xticks(size=16,fontname='Times')
    plt.yticks(size=16,fontname='Times')
    plt.ylim(yminC, ymaxC)
    plt.savefig(Fig4name, format='png', bbox_inches='tight')

    os.system("pdflatex Summary.tex")
    pdfname = os.path.join(WDIR,'Summary.pdf')
    os.rename("Summary.pdf", pdfname )
    
    Fig1name = os.path.join(WDIR,'NormalisedDataOriginalTime.png')
    Fig2name = os.path.join(WDIR,'NormalisedDataStandardisedTime.png')
    Fig3name = os.path.join(WDIR,'RawDataOriginalTime.png')
    Fig4name = os.path.join(WDIR,'RawDataStandardisedTime.png')
    
    os.rename("NormalisedDataOriginalTime.png", Fig1name)
    os.rename("NormalisedDataStandardisedTime.png", Fig2name)
    os.rename("RawDataOriginalTime.png", Fig3name)
    os.rename("RawDataStandardisedTime.png", Fig4name)
