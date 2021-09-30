import struct
import numpy as np
import nibabel as nib
import gzip

class Image4D:
    """
    """
    def __init__ (self, fnimg, _type = 'dicom'):
        self.Filename = fnimg
        self.ImageType = _type

        # Read file into buffer
        self.Data = ImageData()
        if (_type == 'dicom'):
            self.__loadDicom()
        elif (_type == 'nifti'):
            self.__loadNIfTI()
        else:
            print(f'Image4D: Unknown Image Type "{_type}"!')
            
        
    # Load Dicom Data
    def __loadDicom (self):
        data = self.Data

        # Get Header Info
        f = open(self.Filename, 'rb')
        
        # Skip first 128 Bytes
        f.seek(128)
        
        # DICOM header
        DICM = f.read(4).decode('utf-8')
        print("Header: ", DICM)

        # The data tag (7fe0, 0010)
        dataTag = (0x7fe0, 0x0010)
        
        LoopCnt = 0

        # Initialize current tag
        tag = (0x0000, 0x0000)

        while (tag != dataTag and LoopCnt < 200):
            LoopCnt += 1

            tag = (int.from_bytes(f.read(2), byteorder = 'little'), int.from_bytes(f.read(2), byteorder = 'little'))

            code = f.read(2).decode('utf-8')
            n = int.from_bytes(f.read(2), byteorder = 'little')
            # print('(', '{:04x}'.format(tag[0]), '{:04x}'.format(tag[1]), ')', 'code = ', code, '\t n = ', n)
            
            if (tag[0] == 0x0018):
                if (tag[1] == 0x602c):
                    data.deltaX = struct.unpack('<d', f.read(n))[0] * 10
                    # print('deltaX = ', data.deltaX)
                elif (tag[1] == 0x602e):
                    data.deltaY = struct.unpack('<d', f.read(n))[0] * 10
                    # print('deltaY = ', data.deltaY)
                elif (tag[1] == 0x1063):
                    # print('deltaT size = ', n)
                    data.deltaT = float(f.read(n).decode('utf-8'))
                    # print('deltaT = ', data.deltaT)
                else:
                    f.read(n)
            elif (tag[0] == 0x0028):
                if (tag[1] == 0x0008):
                    data.dimT = int(f.read(n).decode('utf-8'))
                    # print('dimT = ', data.dimT)
                elif (tag[1] == 0x0010):
                    data.dimY = int.from_bytes(f.read(n), byteorder = 'little')
                    # print('dimY = ', data.dimY)
                elif (tag[1] == 0x0011):
                    data.dimX = int.from_bytes(f.read(n), byteorder = 'little')
                    # print('dimX = ', data.dimX)
                else:
                    f.read(n)
            elif (tag[0] == 0x3001):
                if (tag[1] == 0x1001):
                    data.dimZ = int.from_bytes(f.read(n), byteorder = 'little')
                    # print('dimZ = ', data.dimZ)
                elif (tag[1] == 0x1003):
                    data.deltaZ = struct.unpack('<d', f.read(n))[0] * 10
                    # print('deltaZ = ', data.deltaZ)
                else:
                    f.read(n)
            else:
                f.read(n)

            if (code == 'OB'):
                f.read(6)

        # Load Data
        n = int.from_bytes(f.read(4), byteorder = 'little')
        expectedSize = data.dimX * data.dimY * data.dimZ * data.dimT
        
        # Data length should match total # of voxels
        assert(n == expectedSize)

        raw = f.read(n)
        assert(len(raw) == expectedSize)
        
        # Import buffer into a numpy array
        # We need to reverse the dimension order (t*z*y*x) to load the 1D buffer,
        # and transpose the array (x*y*z*t) to adjust to the order that numpy 
        # organizes array dimensions
        bufferArr = np.frombuffer(raw, np.uint8) \
            .reshape((data.dimT, data.dimZ, data.dimY, data.dimX))
        bufferArr = np.transpose(bufferArr)

        data.voxels = bufferArr
        data.printInfo()

        f.close()

    # Load NIfTI Data
    def __loadNIfTI (self):
        data = self.Data

        # Load 4D nifti data from file
        buffer = nib.load(self.Filename)
        #print(buffer)
        
        hdr = buffer.header
        data.affine = buffer.affine
        dim = hdr.get_data_shape()
        assert(len(dim) == 4)

        data.dimX = dim[0]
        data.dimY = dim[1]
        data.dimZ = dim[2]
        data.dimT = dim[3]

        pixdim = hdr['pixdim']
        data.deltaX = pixdim[1]
        data.deltaY = pixdim[2]
        data.deltaZ = pixdim[3]
        data.deltaT = pixdim[4]

        print("Voxel datatype: ", buffer.get_data_dtype())
        data.voxels = np.array(buffer.get_fdata()).astype(buffer.get_data_dtype())

        data.printInfo()


    # Export all frames from buffer
    def Export4D (self, outfn):
        print("Exporting 4D Image to: ", outfn)
        data = self.Data
        affine = data.GetAffine()
        img = nib.Nifti1Image(data.voxels, affine)
        # print(img)

        # Save image to the file
        nib.save(img, outfn)


    # Export specific frame
    def ExportFrame (self, frameNum, outfn):
        print("Exporting Frame ", frameNum, " to: ", outfn)
        data = self.Data
        affine = data.GetAffine()
        #print(affine)
        # Array index (0 based) of the framenum (1 based) is always 1 smaller
        voxelArr = np.transpose(np.transpose(data.voxels)[frameNum - 1])
        img = nib.Nifti1Image(voxelArr, affine)

        # Save image to the file
        nib.save(img, outfn)


    def printInfo (self):
        print("Class: Dicom4DReader")
        print("Filename: ", self.Filename)


class ImageData:
    """
    """
    def __init__ (self):
        # Voxel Data
        # - Type: np array with shape [dimX, dimY, dimZ, dimT]
        self.voxels = None
        # Dimension
        self.dimX = None
        self.dimY = None
        self.dimZ = None
        self.dimT = None
        # Spacing
        self.deltaX = None
        self.deltaY = None
        self.deltaZ = None
        self.deltaT = None # frame time

        # Affine
        self.affine = None
    

    def GetAffine (self):
        if self.affine is not None:
            return self.affine
        else:
            # Export in LPS by default
            return np.diag([self.deltaX, self.deltaY, -self.deltaZ, 1])

    def printInfo (self):
        print("Class: ImageData")
        print("Dimension: ", self.dimX, "x", self.dimY, "x", self.dimZ, "x", self.dimT)
        print("Spacing: [", self.deltaX, ", ", self.deltaY, ", ", self.deltaZ, ", ", self.deltaT, "]")

    