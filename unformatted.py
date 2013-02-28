"""
    This module defines a class FortranBinary for interaction with binary files 
    generated by FORTRAN unformatted I/O
"""
import struct, sys

#
#class file:
#    eorskip = 8
#    def __init__(self, name, status=None):
#        self.name = name
#        if status == 'new':
#            self.file = open(name,'w')
#            self.data = ""
#            self.size = 0
#            self.loc = 0
#        else:
#            self.file = open(name,'rw')
#            self.data = self.file.read()
#            self.size = len(self.data)
#            self.loc = 0
#
#    def rewind(self):
#        self.loc = 0
#
#    def close(self):
#        self.file.close()
#
#    def writebuf(self,c,data):
#        n=len(data)
#        for i in range(n):
#            self.data += struct.pack(c,data[i])
#        self.size = len(self.data)
#        self.loc = self.size
#        self.file.write(self.data)
#
#    def readbuf(self, n, c, eor=False):
#        start, stop = self.loc,self.loc+struct.calcsize(c*n)
#        vec = struct.unpack(c*n,self.data[start:stop])
#        self.loc = stop
#        if eor:
#            self.loc += file.eorskip
#        return vec
#
#    def find(self,label):
#        import string
#        bytes = struct.unpack('c'*self.size,self.data)
#        filestring = "".join(bytes)
#        self.loc = string.find(filestring,label)
#        if self.loc < 0:
#            raise NameError("No " + label + " on " + self.name)
#
class EOF(Exception):
    """Defines a simple end of file exception"""
    def __init__(self, filename):
        Exception.__init__(self)
        self.filename = filename
    def __str__(self):
        return "Trying to read past end of file: %s" % self.filename

class FortranBinary():
    """Class for binary files compatible with Fortran Unformatted I/O"""
    pad = 4
    def __init__(self, name, status=None):
        self.name = name
        if status == 'new':
            self.file = open(name,'wb')
        else:
            try:
                self.file = open(name, 'rwb', 10)
            except(IOError):
                print "%s: file %s not found" % (__name__ + '.' + self.__class__.__name__, name)
                raise IOError
                sys.exit(1)
        self.data = None
        self.loc = 0
        self.reclen = 0

    def __iter__(self):
        return self

    def next(self):
        return self.readrec()

    def readrec(self):
        """Read a Fortran record"""
        head = self.file.read(self.pad)
        if len(head) != self.pad:
            raise EOF(self.name)
        size = struct.unpack('i', head)[0]
        self.data = self.file.read(size)
        self.reclen = size
        tail = self.file.read(self.pad)
        assert head == tail
        self.loc = 0
        return self.data


    def readbuf(self, n, c):
        """Read data from current record"""
        start, stop = self.loc, self.loc+struct.calcsize(c*n)
        vec = struct.unpack(c*n, self.data[start:stop])
        self.loc = stop
        return vec

    def find(self, label):
        """Find string label in file"""
        try:
            for rec in self:
                if label in rec: break
        except(EOF):
            print "%s not found on file %s" % (label, self.name)
            sys.exit(1)
        return rec

    def close(self):
        """Close file"""
        self.file.close()

if __name__ == "__main__":
    pass
