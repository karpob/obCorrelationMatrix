from scipy.io import FortranFile
import numpy as np

class gsiCovariance:
    def __init__(self, filename ):
        self.fname = filename
    def setName(self, filename):
        """
        set the name of the file to read/write
        """
        self.fname = filename
    def read(self):
        """
        read data fortran binary. 
        """
        f = FortranFile(self.fname)
        nchan, type = f.read_ints(np.int32)
        self.channelIdx = f.read_ints(np.int32)
        if(type == 8):
            self.R = f.read_reals(np.float64).reshape((nchan,nchan), order="F")
        else:
            self.R = f.read_reals(np.float32).reshape((nchan,nchan), order="F")
        f.close()
    def get(self):
        """
        Return the actual values from the file.
        """
        return self.channelIdx, self.R
    def set(self, idx, R):
        """
        set the channel idx and R for the file.
        """
        self.channelIdx = idx
        self.R = R   
    def write(self):
        """
        write the binary file. 
        """
        f = FortranFile( self.fname, 'w' )
        nchan = len(self.channelIdx)
        f.write_record( np.array([nchan, 8], dtype=np.int32) )
        f.write_record(self.channelIdx)
        f.write_record(self.R.T)
        f.close()

iasi = gsiCovariance('iasi_metop_sea_rcov.bin')
iasi.read()
idx,R = iasi.get()
print(idx,R)
print(R.max(), R.min())
iasi.setName('my_iasi.bin')
iasi.write()

airs = gsiCovariance('airs_aqua_sea_rcov.bin')
airs.read()
idx,R = airs.get()
print(idx,R)
print(R.max(),R.min())
airs.setName('my_airs.bin')
airs.write()


