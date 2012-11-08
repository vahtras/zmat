import struct,sys
class file:
   eorskip=8
   def __init__(self,name,status=None):
      self.name=name
      if status == 'new':
         self.file=open(name,'w')
         self.data=""
         self.size=0
         self.loc=0
      else:
         self.file=open(name,'rw')
         self.data=self.file.read()
         self.size=len(self.data)
         self.loc=0
   def rewind(self):
      self.loc=0
   def close(self):
      self.file.close()
   def writebuf(self,c,data):
      n=len(data)
      for i in range(n):
         self.data+=struct.pack(c,data[i])
      self.size=len(self.data)
      self.loc=self.size
      self.file.write(self.data)
   def readbuf(self,n,c,eor=False):
      #print "readbuf n c",n,c
      start,stop=self.loc,self.loc+struct.calcsize(c*n)
      #print "start,stop",start,stop
      vec=struct.unpack(c*n,self.data[start:stop])
      self.loc=stop
      if eor:
	 self.loc+=file.eorskip
      return vec
   def find(self,label):
      import string
      bytes=struct.unpack('c'*self.size,self.data)
      filestring="".join(bytes)
      self.loc=string.find(filestring,label)
      if self.loc < 0:
         raise NameError("No " + label + " on " + self.name)
class EOF(Exception):
   def __init__(self,filename):
      self.filename=filename
   def __str__(self):
      return "Trying to read past end of file: %s"%self.filename
class fortranbinary:
   pad=4
   def __init__(self,name,status=None):
      self.name=name
      if status == 'new':
         self.file=open(name,'wb')
      else:
	 try:
            self.file=open(name,'rwb',10)
	 except(IOError):
	    print "%s: file %s not found"
	    sys.exit(1)
      self.data=None
      self.loc=0
   def readrec(self):
      head=self.file.read(self.pad)
      if len(head) != self.pad:
	 raise EOF(self.name)
      size=struct.unpack('i',head)[0]
      self.data=self.file.read(size)
      self.reclen=size
      tail=self.file.read(self.pad)
      assert head == tail
      self.loc=0
      return self.data
   def readbuf(self,n,c,eor=False):
      #print "readbuf n c",n,c
      start,stop=self.loc,self.loc+struct.calcsize(c*n)
      #print "start,stop",start,stop
      vec=struct.unpack(c*n,self.data[start:stop])
      self.loc=stop
      return vec
   def find(self,label):
      import string
      while True:
	 try:
	    rec=self.readrec()
	    if string.find(rec,label) > -1:
	       break
	 except EOF:
	    print "%s not found on file %s"%(label,self.name)
	    sys.exit(1)
      return rec

   def close(self):
      self.file.close()

if __name__ == "__main__":
   if 0:
      sirifc=file("SIRIFC")
      sirifc.find("SIR IPH ")
      dbl=sirifc.readbuf(4,'d')
      print dbl
   if 0:
      sirifc=largefile("SIRIFC")
      sirifc.find("SIR IPH ")
      dbl=sirifc.readbuf(4,'d')
      print dbl
   if 0:
      # write
      a=file("yo",'new')
      x=[1.0, 2.0, 3.0]
      print "Writing x",x
      a.writebuf('d',x)
      a.close()
      # read
      b=file("yo")
      y=b.readbuf(3,'d')
      print "Reading y",y
      b.close()
   if 0:
      # write
      a=file("yo",'new')
      s="abcdefghijklmnopqrstuvxyz"
      a.writebuf('c',s)
      a.close()
      # read
      b=largefile("yo")
      #b.rewind()
      y=b.readbuf(25,'c')
      print y
      b.close()
   if 0:
      #1ew
      a=fortranbinary("EOR")
      r1=a.readrec()
      print r1
      r2=a.readrec()
      print r2
      try:
	 r3=a.readrec()
      except EOF,err:
	 print err
      a.close()
   if 1:
      a=fortranbinary("EOR")
      rec=a.find("there")
      print rec
      yo=a.find("yo")
      print yo


