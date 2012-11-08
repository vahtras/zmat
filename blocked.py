import full
class matrix:
    """ Blocked matrix class based on lists of full matrices"""

    def __init__(self,nrow,ncol):
        """ Constructur of the class. 
        Example
        >>> TR=matrix([2,1],[2,1])

        # Dimensions of subblocks
        >>> print TR.nrow, TR.ncol
        [2, 1] [2, 1]

        # Offset dimensions of subblocks
        >>> print TR.irow, TR.icol
        [0, 2] [0, 2]

        # Row and column rank equal
        >>> matrix([3,2],[1])
        Traceback (most recent call last):
        AssertionError

        """

        assert ( len(nrow) == len(ncol) )
        self.nblocks=len(nrow)
        self.nrow=nrow
        self.ncol=ncol
        self.subblock=[]
        self.irow=[]
        self.icol=[]
        for i in range(self.nblocks):
            self.subblock.append(full.matrix((nrow[i],ncol[i])))
            self.irow.append(sum(self.nrow[:i]))
            self.icol.append(sum(self.ncol[:i]))

    def __str__(self):
        """ Formatted output based on full matrix class
        >>> TR=matrix([1],[1])
        >>> print TR
        <BLANKLINE>
        Block 1
        <BLANKLINE>
         (1, 1) 
                      Column   1
        <BLANKLINE>
        """

        retstr=""
        for i in range(self.nblocks):
            if (self.nrow[i]*self.ncol[i]):
                retstr+="\nBlock %d\n"%(i+1) + str(self.subblock[i])
        return retstr

    def __getitem__(self,n):
        """ Index argument returns subblock
        Example
        >>> M = matrix([2], [2])
        >>> print M[0]
        <BLANKLINE>
         (2, 2) 
                      Column   1    Column   2
        <BLANKLINE>
        """
        return self.subblock[n]

    def random(self):
        """ Fill matrix subblocks with random numbers
        Example
        >>> M = matrix([2], [2]).random()
        >>> assert M.subblock[0][0,0] < 1 and M.subblock[0][0,0] > 0
        >>> assert M.subblock[0][1,0] < 1 and M.subblock[0][1,0] > 0
        >>> assert M.subblock[0][0,1] < 1 and M.subblock[0][0,1] > 0
        >>> assert M.subblock[0][1,1] < 1 and M.subblock[0][1,1] > 0
        """
        for i in range(self.nblocks):
            self.subblock[i].random()
        return self

    def __mul__(self,other):
        """Multiplication blockwise
        Example:
        >>> M = matrix([6, 5], [4, 2]) * matrix([4, 2], [3, 1])
        >>> print M.nrow, M.ncol
        [6, 5] [3, 1]
        """
        
        new=matrix(self.nrow,other.ncol)
        for i in range(self.nblocks):
           if self.nrow[i]:
              new.subblock[i]=self.subblock[i]*other.subblock[i]
        return new

    def __rmul__(self,other):
        """Scalar multiplication
        Example:
        >>> A = matrix([1],[1]); A.subblock[0][0,0] = 1; print 2*A
        <BLANKLINE>
        Block 1
        <BLANKLINE>
         (1, 1) 
                      Column   1
               1      2.00000000
        <BLANKLINE>
        """
        new=matrix(self.nrow,self.ncol)
        for i in range(self.nblocks):
           if self.nrow[i]:
               new.subblock[i]=other*self.subblock[i]
        return new

    def __add__(self,other):
        """Addition
        Example:
        >>> A = matrix([1], [1]); A.subblock[0][0, 0] = 2; print A + A
        <BLANKLINE>
        Block 1
        <BLANKLINE>
         (1, 1) 
                      Column   1
               1      4.00000000
        <BLANKLINE>
        """
        new=matrix(self.nrow,self.ncol)
        for i in range(self.nblocks):
           if self.nrow[i]:
               new.subblock[i]=self.subblock[i]+other.subblock[i]
        return new

    def __sub__(self,other):
        """Subtraction
        Example:
        >>> A = matrix([1], [1]); A.subblock[0][0, 0] = 2; print A - A
        <BLANKLINE>
        Block 1
        <BLANKLINE>
         (1, 1) 
                      Column   1
        <BLANKLINE>
        """
        new=matrix(self.nrow,self.ncol)
        for i in range(self.nblocks):
            if self.nrow[i]:
                new.subblock[i]=self.subblock[i]-other.subblock[i]
        return new

    def __neg__(self):
        """Negation
        Example:
        >>> A = matrix([1], [1]); A.subblock[0][0, 0] = 3; print -A
        <BLANKLINE>
        Block 1
        <BLANKLINE>
         (1, 1) 
                      Column   1
               1     -3.00000000
        <BLANKLINE>
        """
        new=matrix(self.nrow,self.ncol)
        for i in range(self.nblocks):
            if self.nrow[i]:
                new.subblock[i]=-self.subblock[i]
        return new

    def __div__(self,other):
        "Solve linear equations"""
        new=matrix(self.nrow,self.ncol)
        for i in range(self.nblocks):
            if self.nrow[i]:
                if isinstance(other,self.__class__):
                    new.subblock[i]=self.subblock[i]/other.subblock[i]
                else:
                    new.subblock[i]=self.subblock[i]/other
        return new

    def __rdiv__(self,other):
        """Inversion
        Example:
        >>> M = matrix([2,1], [2,1])
        >>> M.subblock[0][0, 0] = 2.0
        >>> M.subblock[0][1, 1] = 1.0
        >>> M.subblock[1][0, 0] = 4.0
        >>> print 1/M
        <BLANKLINE>
        Block 1
        <BLANKLINE>
         (2, 2) 
                      Column   1    Column   2
               1      0.50000000    0.00000000
               2      0.00000000    1.00000000
        <BLANKLINE>
        Block 2
        <BLANKLINE>
         (1, 1) 
                      Column   1
               1      0.25000000
        <BLANKLINE>
        """
        return unit(self.nrow,other)/self

    def pack(self):
        for i in range(self.nblocks):
            assert ( self.nrow[i] == self.ncol[i] )
        new=triangular(self.nrow)
        for i in range(self.nblocks):
            new.subblock[i]=self.subblock[i].pack()
        return new

    def unblock(self):
        #print "unblock r c",self.nrow,self.ncol
        nrows=sum(self.nrow)
        ncols=sum(self.ncol)
        new=full.matrix((nrows,ncols))
        for i in range(self.nblocks):
            try:
                new[self.irow[i]:self.irow[i]+self.nrow[i],self.icol[i]:self.icol[i]+self.ncol[i]]=self.subblock[i]
            except ValueError:
                print "blocked.matrix.unblock:ValueError:"
                print "i row col",i,self.nrow[i],self.ncol[i]
                print "lhs",new
                print "rhs",self.subblock[i]
                import sys
                sys.exit(1)
        return new

    def T(self):
        """Transpose
        Example:
        >>> M=matrix([2],[2]); M.subblock[0][0,1]=1; M.subblock[0][1,0]=2
        >>> print M, M.T()
        <BLANKLINE>
        Block 1
        <BLANKLINE>
         (2, 2) 
                      Column   1    Column   2
               1      0.00000000    1.00000000
               2      2.00000000    0.00000000
        <BLANKLINE>
        Block 1
        <BLANKLINE>
         (2, 2) 
                      Column   1    Column   2
               1      0.00000000    2.00000000
               2      1.00000000    0.00000000
        <BLANKLINE>
        """
        new=matrix(self.ncol,self.nrow)
        for i in range(self.nblocks):
            new.subblock[i]=self.subblock[i].T
        return new

    def func(self,f):
         """ Blockwise function of matrix"""
         new=matrix(self.ncol,self.nrow)
         for i in range(self.nblocks):
             new.subblock[i]=self.subblock[i].func(f)
         return new

    def tr(self):
        """Sum blockwise traces
        Example:
        >>> M = matrix([2, 1], [2, 1])
        >>> M.subblock[0][0, 0] = 3
        >>> M.subblock[0][1, 1] = 2
        >>> M.subblock[1][0, 0] = 1
        >>> print M.tr()
        6.0
        """

        sum=0
        for i in range(self.nblocks):
            if self.nrow[i]:
                sum+=self.subblock[i].tr()
        return sum

    def eigvec(self):
        u=matrix(self.nrow,self.nblocks*[1])
        v=matrix(self.nrow,self.ncol)
        for i in range(self.nblocks):
           u.subblock[i],v.subblock[i]=self.subblock[i].eigvec()
        return u,v

    def qr(self):
        q=matrix(self.nrow,self.nrow)
        r=matrix(self.nrow,self.ncol)
        for i in range(self.nblocks):
          q.subblock[i],r.subblock[i]=self.subblock[i].qr()
        return q,r

    def gs(self,S):
        new=matrix(self.nrow,self.ncol)
        Sbl=S.block(self.nrow,self.nrow)
        for i in range(self.nblocks):
           new.subblock[i]=self.subblock[i].gs(Sbl.subblock[i])
        return new

         
      
def unit(nbl,factor=1):
   new=matrix(nbl,nbl)
   for i in range(len(nbl)):
      if nbl[i]:
         new.subblock[i]=full.unit(nbl[i],factor)
   return new

class triangular:
   def __init__(self,dim):
      self.nblocks=len(dim)
      self.dim=dim
      self.subblock=[]
      for i in range(self.nblocks):
         self.subblock.append(full.triangular((dim[i],dim[i])))
   def __str__(self):
      retstr=""
      for i in range(self.nblocks):
         if (self.dim[i]):
            retstr+="\nBlock %d\n"%(i+1) + str(self.subblock[i])
      return retstr
   def random(self):
      for i in range(self.nblocks):
         self.subblock[i].random()
   def __add__(self,other):
      new=triangular(self.dim)
      for i in range(self.nblocks):
         new.subblock[i]=self.subblock[i]+other.subblock[i]
      return new
   def __sub__(self,other):
      new=triangular(self.dim)
      for i in range(self.nblocks):
         new.subblock[i]=self.subblock[i]-other.subblock[i]
      return new
   def unpack(self):
      new=matrix(self.dim,self.dim)
      for i in range(self.nblocks):
         new.subblock[i]=self.subblock[i].unpack()
      return new
   def unblock(self):
      return self.unpack().unblock().pack()
            
      
      
if __name__ == "__main__":
   nbas=(2,1,0)
   norb=(2,1,0)
   from full import header
   #
   # init
   #
   if 1:
      header('init')
      f=matrix(nbas,norb)
   if 0:
      header('random')
      f.random()
   #
   # print
   #
   if 0:
      header('str')
      print f
   if 0:
      header('tr')
      print f.T()

   if 0:
      header('mul')
      print f.T()*f
   if 0:
      header('rmul')
      print f,2*f
   if 0:
      header('add')
      print f+f
   if 0:
      header('sub')
      print f-f
   if 0:
      header('neg')
      print -f
   if 0:
      header('unit')
      print unit((3,2),2)
   if 0:
      header('div')
      print f/f
      print f/2
   if 0:
      header('rdiv')
      print 2./f
   #
   if 0:
      header('unblock')
      print f.unblock()
   if 0:
      header('func')
      import math
      f=lambda x: 1.0/math.sqrt(x)
      a=matrix((3,8),(3,8)).random()
      print a
      S=a.T()*a
      T=S.func(f)
      print "S T TST",S,T,T*S*T
      S=S.unblock()
      T=T.unblock()
      print "S T TST",S,T,T*S*T
   if 0:
      header ('func2')
      import math
      f=lambda x: 1.0/math.sqrt(x)
      a=full.matrix((11,11)).random()
      S=a.T*a
      Sbl=S.block((3,8),(3,8))
      Tbl=Sbl.func(f)
      print "S T TST",Sbl,Tbl,Tbl*Sbl*Tbl
#
# Blocked triangular tests
#
   #
   # init
   #
   if 0:
      header('init')
      tr=triangular((3,2,1))
   #
   # print
   #
   if 0:
      header('print')
      print tr
   #
   #
   #
   if 0:
      header('random')
      tr.random()
      print tr

   #print tr+tr
   #g=matrix(norb,norb)
   #print g.pack()
   #print tr.triangular2square()
   #dim=[2,1]
   #f=triangular(dim)
   ##print f
   #print f.unblock()
