import full
class matrix:
   def __init__(self,nrow,ncol):
      self.rowblocks=len(nrow)
      self.colblocks=len(ncol)
      self.nrow=nrow
      self.ncol=ncol
      self.subblock=[]
      self.irow=[]
      self.icol=[]
      for i in range(self.rowblocks):
	 self.subblock.append([])
	 self.irow.append(sum(self.nrow[:i]))
	 self.icol.append(sum(self.ncol[:i]))
	 for j in range(self.colblocks):
	    self.subblock[i].append([])
	    self.subblock[i][j]=full.matrix((nrow[i],ncol[j]))
   def __str__(self):
      retstr=""
      for i in range(self.rowblocks):
	 for j in range(self.colblocks):
	    retstr+="\nBlock (%d,%d)\n"%(i+1,j+1) + str(self.subblock[i][j])
      return retstr
   def T(self):
      new=matrix(self.ncol,self.nrow)
      for i in range(self.rowblocks):
	 for j in range(self.colblocks):
	    new.subblock[i][j]=self.subblock[j][i].transpose()
      return new
   def __mul__(self,other):
      new=matrix(self.nrow,other.ncol)
      for i in range(self.rowblocks):
	 for j in range(other.colblocks):
	    if self.nrow[i]*other.ncol[j]:
	       for k in range(self.colblocks):
		  new.subblock[i][j]=self.subblock[i][k]*other.subblock[k][j]
      return new
   def random(self):
      for i in range(self.rowblocks):
	 for j in range(self.colblocks):
	    self.subblock[i][j].random()
   def unblock(self):
      #print "unsubblock r c",self.nrow,self.ncol
      nrows=sum(self.nrow)
      ncols=sum(self.ncol)
      new=full.matrix((nrows,ncols))
      for i in range(self.rowblocks):
	 for j in range(self.colblocks):
	    try:
	       new[self.irow[i]:self.irow[i]+self.nrow[i],self.icol[j]:self.icol[j]+self.ncol[j]]=self.subblock[i][j]
	    except ValueError:
	       print "ValueError:"
	       print "i,j",i,j
	       print "irow,nrow",self.irow[i],self.nrow[i]
	       print "icol,ncol",self.icol[j],self.ncol[j]
	       print "lhs", new[self.irow[i]:self.irow[i]+self.nrow[i],self.icol[j]]
	       print "rhs",self.subblock[i][j]
	       import sys
	       sys.exit(1)
      return new

if __name__ == "__main__":
   nbas=(4,2,1)
   norb=(2,1)
   from full import header
   #
   # init
   #
   if 1:
      header('init')
      f=matrix(nbas,norb)
   if 1:
      header('random')
      f.random()
   if 1:
      header('__str__')
      print f
   if 1:
      header('unblock')
      print f.unblock()

