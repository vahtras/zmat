import math,types,sys
import numpy

class matrix(numpy.ndarray):
   """ A subclass of numpy.ndarray for matrix syntax and better printing """
   fmt="%14.8f"
   order='F'

   def __new__(subtype, shape,order=order,fmt=None):
      """Constructor ..."""            
      #print "matrix.__new__",subtype,shape,order,fmt
      obj=numpy.zeros(shape,order=matrix.order).view(subtype)
      if fmt is None: obj.fmt=matrix.fmt
      return obj

   def __array_finalize__(self,obj):
      #print "matrix.__array_finalize__",type(self),type(obj)
      if obj is None: return
      self.fmt=getattr(obj,'fmt',matrix.fmt)
      self._I=None

   def debug(self):
      """ Called for various error conditinos"""
      print "type ",type(self)
      print "shape ",self.shape
      print "dtype", self.dtype
      print "strides", self.strides
      print "order", self.order
      print "fmt",self.fmt

   def __str__(self):
      """Output formatting of matrix object, inspired by Dalton OUTPUT subroutine
      Example:
      >>> M=matrix((2,2)); M[0,0]=M[1,1]=1; print M
      <BLANKLINE>
       (2, 2) 
                    Column   1    Column   2
             1      1.00000000    0.00000000
             2      0.00000000    1.00000000
      <BLANKLINE>
      """

      #self.debug()
      retstr='\n %s \n'%str(self.shape)
      if len(self.shape) == 1:
         r=self.shape[0]
         #retstr+="(%d)\n"%(r)
         if 0:
            for i in range(r):
               retstr=retstr + self.fmt % self[i]
               retstr=retstr+'\n'
         else:
            columnsperblock=1
            fullblocks=1
            trailblock=0
            for b in range(fullblocks):
               crange=range(b*columnsperblock,(b+1)*columnsperblock)
               retstr+=" "*10
               for j in crange:retstr+="    Column%4d"%(j+1)
               retstr+='\n'
               for i in range(r):
                  rownorm=math.fabs(self[i])
                  if rownorm > 1e-8:
                     retstr+="%8d  "%(i+1)
                     retstr=retstr + self.fmt % self[i]
                     retstr+='\n'
               retstr+='\n'
      elif len(self.shape) == 2:
         r,c=self.shape
         #retstr+="(%d,%d)\n"%(r,c)
         if 0:
            for i in range(r):
               for j in range(c):
                  retstr=retstr + self.fmt % self[i,j]
               retstr=retstr+'\n'
         else:
            columnsperblock=5
            fullblocks=c/columnsperblock
            trailblock=c%5
            for b in range(fullblocks):
               crange=range(b*columnsperblock,(b+1)*columnsperblock)
               retstr+=" "*10
               for j in crange:retstr+="    Column%4d"%(j+1)
               retstr+='\n'
               for i in range(r):
                  rownorm=self[i,crange].norm2()
                  if rownorm > 1e-8:
                     retstr+="%8d  "%(i+1)
                     for j in crange:
                        retstr=retstr + self.fmt % self[i,j]
                     retstr+='\n'
               retstr+='\n'
            crange=range(fullblocks*columnsperblock,c)
            if trailblock:
               retstr+=" "*10
               for j in crange:retstr+="    Column%4d"%(j+1)
               retstr+='\n'
               for i in range(r):
                  rownorm=self[i,crange].norm2()
                  if rownorm > 1e-8:
                     retstr+="%8d  "%(i+1)
                     for j in crange:
                        retstr=retstr + self.fmt % self[i,j]
                     retstr=retstr+'\n'
      elif len(self.shape) > 2:
         r,c=self.shape[:2]
         hishape=self.shape[2:]
         losize=r*c
         if losize == 0:
            return "\nZero dimension\n"
         hisize=self.size/losize
         altshape=(r,c,hisize)
         #print altshape
         alt=self.reshape(altshape,order=matrix.order)
         #
         # Given linear index n find tuple idx=(i0,i1,i2...) for actual shape hishape=(n0,n1,n2..)
         #
         for n in range(hisize):
            k=n
            idx=[]
            for i in range(len(hishape),1,-1):
               #
               # integer divide by n0*...n(i-2), save remainder
               #
               dp=numpy.asarray(hishape[:i-1]).prod()
               idx.append(k/dp)
               #print "idx", idx
               k=k%dp
            idx.append(k)
            idx.reverse()

            retstr+=str(idx) + str(alt[:,:,n])
      elif len(self.shape) == 0: #returned by numpy.max
         return self.fmt%self.sum()
         

      return retstr

   def __mul__(self,other):
      """Matrix multiplication
      Example:
      >>> M = matrix((2,2)); M[0,0]=M[1,1]=2.0; M[0,1]=M[1,0]=1.0
      >>> print M*M
      <BLANKLINE>
       (2, 2) 
                    Column   1    Column   2
             1      5.00000000    4.00000000
             2      4.00000000    5.00000000
      <BLANKLINE>

      >>> print 2*M
      <BLANKLINE>
       (2, 2) 
                    Column   1    Column   2
             1      4.00000000    2.00000000
             2      2.00000000    4.00000000
      <BLANKLINE>

      >>> print M*2
      <BLANKLINE>
       (2, 2) 
                    Column   1    Column   2
             1      4.00000000    2.00000000
             2      2.00000000    4.00000000
      <BLANKLINE>
      """

      if isinstance(other,self.__class__):
         try:
            return numpy.dot(self,other)
         except ValueError:
            print "full.matrix.__mul__:ValueError",self.shape,other.shape
            raise ValueError
      else:
         return other*self

   def x(self,other):
      """Outer product
      Example: (1,0).T*(0,1)
      >>> v1=matrix((2,))
      >>> v2=matrix((2,))
      >>> v1[0]=v2[1]=1.0
      >>> print v1
      <BLANKLINE>
       (2,) 
                    Column   1
             1      1.00000000
      <BLANKLINE>
      <BLANKLINE>

      >>> print v2
      <BLANKLINE>
       (2,) 
                    Column   1
             2      1.00000000
      <BLANKLINE>
      <BLANKLINE>

      >>> print v1.x(v2)
      <BLANKLINE>
       (2, 2) 
                    Column   1    Column   2
             1      0.00000000    1.00000000
      <BLANKLINE>
      """

      a=self.flatten(matrix.order)
      b=other.flatten(matrix.order)
      c=numpy.outer(self,other).reshape(self.shape+other.shape).view(matrix)
      return c

   def __div__(self,other):
      """Solution of linear equation/inversion
      Example:
      >>> A=matrix((2,2)).random(); x=matrix((2,1)).random()
      >>> b=A*x
      >>> print x-b/A
      <BLANKLINE>
       (2, 1) 
                    Column   1
      <BLANKLINE>
      """

      if isinstance(other,self.__class__):
         new=numpy.linalg.solve(other,self)
      else:
         new=(1.0/other)*self
      return new

   def scatteradd(self,other,rows=None,columns=None):
      self.scatter(other,rows,columns,add=1)

   def scatter(self,other,rows=None,columns=None,add=0):
      #
      # Usage scatter copies element into other matrix
      # with indices other[rows[k],columns[l]]=self[k,l]
      #
      #print "scatter:self,ind,n",self,ind,n
      #print self.cdim
      #print len(ind)
      if not add:
         other.clear()
      r,c=self.shape
      if rows and columns:
         assert(r == len(rows))
         assert(c == len(columns))
         for i in range(r):
            for j in range(c):
               other[rows[i],columns[j]]+=self[i,j]
      else:
         if rows:
            assert(r == len(rows))
            for i in range(r):
               other[rows[i],:]+=self[i,:]
         if columns:
            assert(c == len(columns))
            for j in range(c):
               other[:,columns[j]]+=self[:,j]
      return

   def pack(self,other,rows=None,columns=None,add=0):
      #
      # Usage pack copies element into other matrix
      # with indices other[k,columns[l]=self[rows[k],columns[l]
      #
      if not add:
         other.clear()
      if rows and columns:
         assert(other.rdim == len(rows))
         assert(other.cdim == len(columns))
         for i in range(other.rdim):
            for j in range(other.cdim):
                  other[i,j]+=self[rows[i],columns[j]]
      else:
         if rows:
            assert(self.rdim == len(rows))
            for i in range(self.rdim):
               other[i,:]+=self[rows[i],:]
         if columns:
            assert(self.cdim == len(columns))
            for j in range(self.cdim):
               other[:,j]+=self[:,columns[j]]
      return

   def __sob__(self,other):
      return self + (-other)
#  def __pos__(self):
#     new=self.__class__(self.rdim,self.cdim)
#     new=+self
      return new

   def __nog__(self):
      new=self.__class__(self.rdim,self.cdim)
      new=-self
      return new

   def __and__(self,other):
      if len(self.shape) == 1:
         return (self.T*other)
      else:
         return self.ravel(self.order)*other.ravel(other.order)

   def __xor__(self,other):
      return self*other-other*self

   def inv(self):
      if self._I is None:
         r,c=self.shape
         assert r == c
         self._I=unit(r)/self
      return self._I
   I = property(fget=inv)
   def __rdiv__(self,other):
      r,c=self.shape
      return unit(r,other)/self
   def tr(self):
      return self.trace()

   def det(self):
      return numpy.linalg.det(self)
   def minor(self,i,j):
      r,c=self.shape
      rows=range(r)
      cols=range(c)
      rows.remove(i)
      cols.remove(j)
      # bug
      #print self[rows,cols]
      # do rows and cols separately
      redr=self[rows,:]
      #print 'redr',redr
      redc=redr[:,cols]
      #print 'redc',redc
      return redc

   def cofactor(self):
      r,c=self.shape
      assert(r == c)
      new=matrix((r,c))
      for i in range(r):
         for j in range(c):
            new[i,j]=(-1)**(i+j)*self.minor(i,j).det()
      return new
   def eig(self):
      """
      Return sorted eigenvalues as a column matrix
      """
      r,c=self.shape; assert r == c;
      eigvals=numpy.linalg.eigvals(self)
      p=eigvals.argsort()
      return eigvals[p].view(matrix)
   def eigvec(self):
      """
      Return eigenvalue/eigenvector pair, sorted
      by eigenvalue
      """
      r,c=self.shape; assert r == c;
      U,V=numpy.linalg.eig(self)
      p=U.argsort()
      #
      # Note that eig returns U as ndarray but V as its subclass
      #
      return U[p].view(matrix),V[:,p]
   def qr(self):
      """
      Return eigenvalue/eigenvector pair, sorted
      by eigenvalue
      """
      r,c=self.shape; assert r == c;
      Q,R=numpy.linalg.qr(self)
      #
      # Note that eig returns U as ndarray but V as its subclass
      #
      return Q,R
   def normalize(self,S=None):
      if S is None:
         norm=1/math.sqrt(self&self)
      else:
         norm=1/math.sqrt(self&(S*self))
      self*=norm
   def GS(self,S,T=False):
      r,c=self.shape
      new=matrix((r,c))
      new[:,0]=self[:,0]
      new[:,0].normalize(S)
      for i in range(1,c):
         P=new[:,:i]*new[:,:i].T*S
         new[:,i]=self[:,i]-P*self[:,i]
         new[:,i].normalize(S)
      #
      # relation new=self*T
      #          self.T*new=self.T*self*T
      #          (self.T*self).inv()*self.T*new = T (QR?)
      #
      if T:
         return (new.T*S*self).inv()
      else:
         return new
   def GST(self,S):
      return self.GS(S,T=True)
         

   def lowdin(self,S):
      invsq=lambda x: 1.0/math.sqrt(x)
      Si12=scipy.linalg.sqrtm(S.I)
      return

   def sqrt(self):
      from scipy.linalg import sqrtm
      return sqrtm(self).real.view(matrix)

   def invsqrt(self):
      from scipy.linalg import sqrtm
      return sqrtm(self.I).real.view(matrix)

   def vec2diag(self):
      n=self.shape[0]
      new=matrix((n*n))
      new[:n*n:n+1]=self[:]
      return new.reshape((n,n))

   def func(self,f):
      if False:
         import scipy.linalg
         return scipy.linalg.funm(self,f)
      else:
         val,vec=numpy.linalg.eig(self)
         #print self,val,vec
         #print "full.func:normaltest",numpy.allclose(self.T*self,self*self.T)
         #print "full.func:unitarytest",vec.T*vec,vec*vec.T
         fval=val.view(matrix)
         #print fval
         n=len(val)
         new=matrix((n,n))
         for i in range(n):
            fval[i]=f(val[i])
            new[i,i]=fval[i]
         #print numpy.mat(numpy.diag(f(val)))
         #print fval
         #print new
         #return vec*new*vec.T
         #not all eigenvectors are returns as unitary
         #even for symmetric matrices: switch to inverse
         return vec*new*vec.inv()
   def exp(self):
      r,c=self.shape
      new=unit(r)
      termnorm=new&new
      term=new*1
      i=0
      while (termnorm > 1e-8):
         i=i+1
         term*=self/i
         new+=term
         termnorm=math.sqrt(term&term)
         #print term,termnorm
      return new
         
   def random(self):
      if 0:
         import random
         for i in range(self.size):
            self.flat[i] = random.random()
      else:
         self.flat=numpy.random.random(self.size)
      return self
   def sym(self):
      """return symmetrized matrix"""
      return .5*(self + self.T)
   def antisym(self):
      return .5*(self - self.T)
   def pack(self,anti=False):
      #raise Exception("test")
      r,c=self.shape
      assert r == c
      _t=triangular(self.shape)
      fac=1
      if anti is True: fac=-1
      for i in range(r):
         for j in range(i+1):
            _t[i,j]=.5*(self[i,j]+fac*self[j,i])
      return _t
   def lower(self):
      assert self.cdim == self.rdim
      _t=triangular(self.shape)
      r,c=self.shape
      for i in range(r):
         for j in range(i+1):
            _t[i,j]=self[i,j]
      return _t
   def fold(self):
      assert self.cdim == self.rdim
      _t=triangular(self.cdim)
      for i in range(self.cdim):
         for j in range(i):
            _t[i,j]=math.sqrt(2)*self[i,j]
         _t[i,i]=self[i,i]
      return _t
   def norm2(self):
      return numpy.linalg.norm(self)
   def block(self,rdim,cdim):
      assert len(rdim) == len(cdim)
      import blocked
      new=blocked.matrix(rdim,cdim)
      rstart=0
      cstart=0
      for i in range(len(rdim)):
         new.subblock[i]=self[rstart:rstart+rdim[i],cstart:cstart+cdim[i]]
         rstart+=rdim[i]
         cstart+=cdim[i]
      return new
   def subblocked(self,rdim,cdim):
      import subblocked
      new=subblocked.matrix(rdim,cdim)
      rstart=0
      for i in range(new.rowblocks):
         cstart=0
         for j in range(new.colblocks):
            #print "rows",rstart,rstart+rdim[i],"cols",cstart,cstart+cdim[j]
            new.subblock[i][j]=self[rstart:rstart+rdim[i],cstart:cstart+cdim[j]]
            cstart+=cdim[j]
         rstart+=rdim[i]
      return new
   def clear(self):
      self*=0.0

   def cross(self, *args):
     """ With out argument return 3x3 matrix Ax, with a vector return ordinary
         cross produce AxB"""
     assert self.shape == (3,)
     if not args:
         new=matrix((3,3))
         new[0,1]=-self[2]
         new[1,0]=self[2]
         new[0,2]=self[1]
         new[2,0]=-self[1]
         new[1,2]=-self[0]
         new[2,1]=self[0]
         return new
     else:
         c = matrix(3)
         a = self
         b, = args
         c[0] = a[1]*b[2] - a[2]*b[1]
         c[1] = a[2]*b[0] - a[0]*b[2]
         c[2] = a[0]*b[1] - a[1]*b[0]
         return c

   def dist(self,other):
      return (self-other).norm2()
   def angle3(self,B,C):
      return (self-B).angle(C-B)
   def angle3d(self,B,C):
      return (self-B).angle(C-B)*180/math.pi
   def angle(self,other):
      dot=self&other
      cos2a=dot*dot/((self&self)*(other&other))
      if (cos2a > 1):
         if cos2a-1>1e-14:
            print "angle:self",self
            print "angle:other",other
            print "angle:dot",dot
            print "angle:cos2a=1+%20.14e"%(cos2a-1)
            raise ValueError
         else:
            #print "full.matrix.angle:cosa reset to 1 due to numerical roundoff error"
            cos2a = 1
      if dot > 0:
         cosa=math.sqrt(cos2a)
      else:
         cosa=-math.sqrt(cos2a)
      return math.acos(cosa)
   def angled(self,other):
      return self.angle(other)*180/math.pi

   def rot(self, angle, vec, origin=None):
      p = vec[:]/vec.norm2()
      if origin is None:
         so = self
      else:
         so = self - origin
       
      sp = p*(p&so)
      sq = so - sp
      sr = p.cross(sq)
      #print "pqr",p, q, r
      self[:] = sp + sq*math.cos(angle) + sr*math.sin(angle)
      if origin is not None:
         self[:] += origin
      return self

   def dihedral(self,r3,r2,r1):
      b3=self-r3
      b2=r3-r2
      b1=r2-r1
      b1xb2=b1.cross()*b2
      b2xb3=b2.cross()*b3
      n2=b2.norm2()
      return math.atan2(
         n2*b1&b2xb3,
         b1xb2&b2xb3
         )
   def dihedrald(self,r3,r2,r1):
      return self.dihedral(r3,r2,r1)*180/math.pi
   def __or__(self,other):
      assert self.rdim == other.rdim
      assert self.cdim == other.cdim
      c=matrix((self.rdim,self.cdim))
      for i in range(self.rdim):
         for j  in range(self.cdim):
            c[i,j]=self[i,j]*other[i,j]
      return c

   def svd(self):
      """Compact singular value decomposition
      input n,p
      output u(n,p)
             s(p,p)
             v(p,p)
      """
      u,s,vt=numpy.linalg.svd(self,full_matrices=0)
      s=numpy.diag(s).view(matrix)
      return  u,s,vt.T

   def sum(self, **kwargs):
      """Was previously handled by numpy.sum but 
         as of ubuntu 12.04 numpy.sum returns matrix([sum]) rather than sum
      """
      return numpy.sum(self.view(numpy.ndarray), **kwargs)

def unit(n,factor=1):
   vec=matrix((n*n,))
   vec[:n*n:n+1]=factor
   return vec.reshape((n,n))
def permute(select,n):
   complement=range(n)
   for i in range(n):
      if i in select:
         complement.remove(i)
   permlist=select+complement
   new=matrix((n,n))
   for i in range(n):
      new[permlist[i],i]=1
   return new
def init(nestlist):
   return numpy.array(nestlist).view(matrix).T
#
# init should be generalized with matrix.order, now 'F' assumed -> traspose
#

class triangular(numpy.ndarray):
   def __new__(subtype, shape,anti=False,fmt=None):
      #print "triangular.__new__:",subtype,shape,anti,fmt
      tshape=((shape[0]*(shape[0]+1))/2,)
      obj=matrix(tshape).view(subtype)
      obj.sshape=shape
      obj.anti=anti
      if fmt is None: obj.fmt=matrix.fmt
      return obj
   def __array_finalize__(self,obj):
      #print "triangular.__array_finalize__:",type(obj)
      if obj is None: return
      self.dim=int(math.sqrt(0.25+2*obj.size))
      #print "dim",self.dim
      #print "shape",self.shape
      self.sshape=getattr(obj,'sshape',(self.dim,self.dim))
      self.anti=getattr(obj,'anti',False)
      self.fmt=getattr(obj,'fmt',matrix.fmt)
      
   @staticmethod
   def init(arr):
      n = int(round(-0.5 + math.sqrt(0.25+2*len(arr))))
      # should test for valid n
      new = triangular((n, n))
      ij = 0
      for i in range(n):
         for j in range(i+1):
             new[i, j] = arr[ij]
             ij += 1
      return new



   def __str__(self):
      retstr="\n"
      r,c=self.sshape
      for i in range(r):
         for j in range(i+1):
            retstr+=self.fmt%(self[i,j])
         retstr+="\n"
      return retstr
   def __getitem__(self,args):
      #print "__getitem__",args
      vec=self.view(matrix)
      i,j=args
      if self.anti and i < j:
         ij=j*(j+1)/2+i
         return -vec[ij]
      else:
         ij=i*(i+1)/2+j
         return vec[ij]
   def __setitem__(self,args,value):
      vec=self.view(matrix)
      i,j=args
      if i < j and self.anti:
         ij=j*(j+1)/2+i
         vec[ij]=-value
      else:
         ij=i*(i+1)/2+j
         vec[ij]=value
   def unpack(self):
      n=self.sshape[0]
      new=matrix((n,n))
      try:
         import pdpack
         if self.anti:
            new=pdpack.daptge(self,new)
         else:
            new=pdpack.dsptsi(self,new)
      except ImportError:
         for i in range(n):
            new[i,i]=self[i,i]
            for j in range(i):
               new[i,j]=self[i,j]
               if self.anti:
                  new[j,i]=-self[i,j]
               else:
                  new[j,i]=self[i,j]
      return new
   def __mul__(self,other):
      if isinstance(other,self.__class__):
         return self.unpack()*other.unpack()
      else:
         return other*self
   def random(self):
      matrix.random(self)
      if self.anti:
         for i in range(self.sshape[0]):
            self[i,i]=0
      return self

         
def header(str):
   print """
   %s
   %s
   """%(str,len(str)*'-')
if __name__ == "__main__":
   #
   # init
   #
   if 0:
      header('init')
      a=matrix((3,3))
   if 0:
      header('fmt')
      print a.fmt
      b=matrix((3,3),fmt="%5.2f")
      print b.fmt
   #
   # str 
   #
   if 0:
      header('str')
      a=matrix(16).random()
      print a
      b=a.reshape((4,4),order=matrix.order)
      print b
      c=a.reshape((2,2,2,2),order=matrix.order)
      print c
      d=matrix((4,4)).random()
      print d
      d=matrix((5,5)).random()
      print d
      d=matrix((6,6)).random()
      print d
   #
   # random (setitem)
   #
   if 1:
      header('random')
      a=matrix((3,3)).random()
      print a
      a=matrix((3,3))
      a.random()
      print a
      #a.debug()
   if 0:
      header('subblocked')
      asb=a.subblocked([2,1],[2,1])
      print asb
   #
   #
   # slice
   #
   if 0:
      header('slice')
      print a[:,:1]
      print a[:,0]
      print a[:1,:]
      print a[0,:].T
   #
   #  neg
   #
   if 0:
      header('neg')
      print -a
   #
   # unit
   #
   if 0:
      header('unit')
      print unit(3)
      print unit(3,3.14)
   #
   #  add
   #
   if 0:
      header('add')
      print a+a
   #
   #   sub
   #
   if 0:
      header('sub')
      print a-a
   #
   #  mul
   #
   if 0:
      header('mul')
      print a*a
   #
   #  rmul
   #
   if 0:
      header('rmul')
      print 2.0*a
      print a*2
   #
   #  inv
   #
   if 1:
      header('inv')
      print a*a.I
   #
   #  div
   #
   if 0:
      header('div')
      print a/a
      print a/2 - 0.5*a
      print (1/a)*a
   #
   #  transpose
   #
   if 0:
      header('transpose')
      #a=matrix((2,2)).random()
      #print a,a.T
      b=matrix((2,2,2,2)).random()
      print b,b.transpose((1,0,2,3))
   #
   #  __and__
   #
   if 1:
      header('__and__')
      a=matrix((3,2)).random()
      print a&a,(a.T*a).tr()

   #
   # eig
   #
   if 0:
      header('eig')
      a2=a.T*a
      print a2,a2.eig()
   #
   # eigvec
   #
   if 0:
      header('eigvec')
      a2=a.T*a
      u,v=a2.eigvec()
      print u,v,v.T*a2*v
   #
   # func
   #
   if 0:
      header('func')
      a=a+a.T
      a2=a*a
      print a2.func(math.sqrt),a

   if 0:
      header('trace')
      print a
      print a.trace()

   if 0:
      header('exp')
      print a.exp()*(-a).exp()

   #
   #gather-scatter
   #
   if 0:
      header('gather columns')
      c=[0,2]
      print a
      ga=a[:,c]
      print  ga
      header('scatter columns')
      sa=matrix((3,3))
      sa[:,c]=ga
      print sa
      #ga.scatter(sa,columns=[0,2])
      header('gather rows')
      ga=a[c,:]
      print ga
      header('scattter rows')
      sa=matrix((3,3))
      sa[c,:]=ga
      print sa

   #determinant
   if 0:
      header('determinant')
      a=matrix((3,3)).random()
      print a
      print a.det()
   #
   # minor
   #
   if 0:
      header('minor')
      a=matrix((3,3)).random()
      print a
      print a.minor(1,1)
   #
   #cofactor
   #
   if 0:
      header('cofactor')
      a=matrix((3,3)).random()
      print a
      print a.cofactor() - a.inv().T*a.det()
   #
   #outer
   #
   if 0:
      header('x(outer)')
      a=matrix((3,2)).random()
      b=matrix((2,1)); b[0,0]=1.; b[1,0]=-1.
      c=a.x(b)
      print a,b,c
   #
   #
   #
   if 0:
      header('clear')
      a=matrix((3,3)).random()
      print a
      a.clear()
      print a

#
# Tensor dot
#     
   if 0:
      header("dot with tensors")
      a=numpy.arange(36).reshape((2,2,3,3),order=matrix.order).view(matrix)
      b=unit(3)
      #
      # a*b gives (2,2,3,:)(:,3)
      # b*a gives (3,:)(2,2,:,3)
      #
      print "a",a
      print "b",b
      c=a*b
      print "a*b",c.shape,c.order,c
      c=b*a
      print "b*a",c.shape,c.order,c
      print "b*a",c.shape,c.order,c.transpose((1,2,0,3))

   if 0:
      header("normalize")
      S=unit(3,2)
      print a,S
      a[:,0].normalize(S)
      print a
   if 0:
      header("GS")
      I=unit(3)
      S=matrix((3,3)).random(); S=S.T*S
      old=I[:,:2]
      new=old.GS(S)
      print "old",old,"old overlap",old.T*S*old
      print "new",new,"new overlap",new.T*S*new
      # new=old*T
      T=old.GST(S)
      print "T",T,"new-old*T",new-old*T
   if 0:
      header('permute')
      print permute([0,1,3,4,5],9)
#
#
# Triangluar tests
#
   if 0:
      header('new')
      t=triangular((4,4)).random()
      print t.anti
      u=triangular((4,4),anti=True).random()
      print u.anti
   if 0:
      header('print')
      print t,u
   if 0:
      header('setitem')
      t[1,1]=1
      t[2,2]=2
      t[3,3]=3
      print t
   if 0:
      header('unpack')
      print t.unpack()
      print u.unpack()
   if 0:
      header('pack')
      print t.unpack().pack()
      print u.unpack().pack(anti=True)
   if 0:
      header('svd')
      a=matrix((5,3)).random()
      u,s,v=a.svd()
      print u,s,v
      print a-u*s*v.T

  
