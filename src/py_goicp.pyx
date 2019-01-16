#cython: c_string_encoding=ascii  # for cython>=0.19
#cython: embedsignature=False
from  libcpp.string  cimport string as libcpp_string
from  libcpp.string  cimport string as libcpp_utf8_string
from  libcpp.string  cimport string as libcpp_utf8_output_string
from  libcpp.set     cimport set as libcpp_set
from  libcpp.vector  cimport vector as libcpp_vector
from  libcpp.pair    cimport pair as libcpp_pair
from  libcpp.map     cimport map  as libcpp_map
from  libcpp cimport bool
from  libc.string cimport const_char
from cython.operator cimport dereference as deref, preincrement as inc, address as address
from  AutowrapRefHolder cimport AutowrapRefHolder
from  AutowrapPtrHolder cimport AutowrapPtrHolder
from  AutowrapConstPtrHolder cimport AutowrapConstPtrHolder
from  smart_ptr cimport shared_ptr
from goicpcc cimport GoICP as _GoICP
from goicpcc cimport POINT3D as _POINT3D
from goicpcc cimport ROTNODE as _ROTNODE
from goicpcc cimport TRANSNODE as _TRANSNODE

cdef extern from "autowrap_tools.hpp":
    char * _cast_const_away(char *) 

cdef class ROTNODE:
    """
    Cython implementation of _ROTNODE
    """

    cdef shared_ptr[_ROTNODE] inst

    def __dealloc__(self):
         self.inst.reset()

    
    property a:
        def __set__(self, float a):
        
            self.inst.get().a = (<float>a)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().a
            py_result = <float>_r
            return py_result
    
    property b:
        def __set__(self, float b):
        
            self.inst.get().b = (<float>b)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().b
            py_result = <float>_r
            return py_result
    
    property c:
        def __set__(self, float c):
        
            self.inst.get().c = (<float>c)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().c
            py_result = <float>_r
            return py_result
    
    property w:
        def __set__(self, float w):
        
            self.inst.get().w = (<float>w)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().w
            py_result = <float>_r
            return py_result
    
    property ub:
        def __set__(self, float ub):
        
            self.inst.get().ub = (<float>ub)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().ub
            py_result = <float>_r
            return py_result
    
    property lb:
        def __set__(self, float lb):
        
            self.inst.get().lb = (<float>lb)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().lb
            py_result = <float>_r
            return py_result
    
    property l:
        def __set__(self,  l):
        
            self.inst.get().l = (<int>l)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().l
            py_result = <int>_r
            return py_result
    
    def __init__(self):
        """Cython signature: void ROTNODE()"""
        self.inst = shared_ptr[_ROTNODE](new _ROTNODE()) 

cdef class TRANSNODE:
    """
    Cython implementation of _TRANSNODE
    """

    cdef shared_ptr[_TRANSNODE] inst

    def __dealloc__(self):
         self.inst.reset()

    
    property x:
        def __set__(self, float x):
        
            self.inst.get().x = (<float>x)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().x
            py_result = <float>_r
            return py_result
    
    property y:
        def __set__(self, float y):
        
            self.inst.get().y = (<float>y)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().y
            py_result = <float>_r
            return py_result
    
    property z:
        def __set__(self, float z):
        
            self.inst.get().z = (<float>z)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().z
            py_result = <float>_r
            return py_result
    
    property w:
        def __set__(self, float w):
        
            self.inst.get().w = (<float>w)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().w
            py_result = <float>_r
            return py_result
    
    property ub:
        def __set__(self, float ub):
        
            self.inst.get().ub = (<float>ub)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().ub
            py_result = <float>_r
            return py_result
    
    property lb:
        def __set__(self, float lb):
        
            self.inst.get().lb = (<float>lb)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().lb
            py_result = <float>_r
            return py_result
    
    def __init__(self):
        """Cython signature: void TRANSNODE()"""
        self.inst = shared_ptr[_TRANSNODE](new _TRANSNODE()) 

cdef class POINT3D:
    """
    Cython implementation of _POINT3D
    """

    cdef shared_ptr[_POINT3D] inst

    def __dealloc__(self):
         self.inst.reset()

    
    property x:
        def __set__(self, float x):
        
            self.inst.get().x = (<float>x)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().x
            py_result = <float>_r
            return py_result
    
    property y:
        def __set__(self, float y):
        
            self.inst.get().y = (<float>y)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().y
            py_result = <float>_r
            return py_result
    
    property z:
        def __set__(self, float z):
        
            self.inst.get().z = (<float>z)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().z
            py_result = <float>_r
            return py_result
    
    def _init_0(self):
        """Cython signature: void POINT3D()"""
        self.inst = shared_ptr[_POINT3D](new _POINT3D())
    
    def _init_1(self, float in_0 , float in_1 , float in_2 ):
        """Cython signature: void POINT3D(float, float, float)"""
        assert isinstance(in_0, float), 'arg in_0 wrong type'
        assert isinstance(in_1, float), 'arg in_1 wrong type'
        assert isinstance(in_2, float), 'arg in_2 wrong type'
    
    
    
        self.inst = shared_ptr[_POINT3D](new _POINT3D((<float>in_0), (<float>in_1), (<float>in_2)))
    
    def __init__(self, *args , **kwargs):
        """
          - Cython signature: void POINT3D()
          - Cython signature: void POINT3D(float, float, float)
"""
        if not args:
             self._init_0(*args)
        elif (len(args)==3) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and (isinstance(args[2], float)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def pointToString(self):
        """Cython signature: void pointToString()"""
        self.inst.get().pointToString() 

cdef class GoICP:
    """
    Cython implementation of _GoICP
    """

    cdef shared_ptr[_GoICP] inst

    def __dealloc__(self):
         self.inst.reset()

    
    property Nm:
        def __set__(self,  Nm):
        
            self.inst.get().Nm = (<int>Nm)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().Nm
            py_result = <int>_r
            return py_result
    
    property Nd:
        def __set__(self,  Nd):
        
            self.inst.get().Nd = (<int>Nd)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().Nd
            py_result = <int>_r
            return py_result
    
    property MSEThresh:
        def __set__(self, float MSEThresh):
        
            self.inst.get().MSEThresh = (<float>MSEThresh)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().MSEThresh
            py_result = <float>_r
            return py_result
    
    property SSEThresh:
        def __set__(self, float SSEThresh):
        
            self.inst.get().SSEThresh = (<float>SSEThresh)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().SSEThresh
            py_result = <float>_r
            return py_result
    
    property icpThresh:
        def __set__(self, float icpThresh):
        
            self.inst.get().icpThresh = (<float>icpThresh)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().icpThresh
            py_result = <float>_r
            return py_result
    
    property trimFraction:
        def __set__(self, float trimFraction):
        
            self.inst.get().trimFraction = (<float>trimFraction)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().trimFraction
            py_result = <float>_r
            return py_result
    
    property inlierNum:
        def __set__(self,  inlierNum):
        
            self.inst.get().inlierNum = (<int>inlierNum)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().inlierNum
            py_result = <int>_r
            return py_result
    
    property doTrim:
        def __set__(self,  doTrim):
        
            self.inst.get().doTrim = (<bool>doTrim)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().doTrim
            py_result = <bool>_r
            return py_result
    
    def setInitNodeTrans(self, TRANSNODE in_0 ):
        """Cython signature: void setInitNodeTrans(TRANSNODE &)"""
        assert isinstance(in_0, TRANSNODE), 'arg in_0 wrong type'
    
        self.inst.get().setInitNodeTrans((deref(in_0.inst.get())))
    
    def optimalTranslation(self):
        """Cython signature: libcpp_vector[double] optimalTranslation()"""
        _r = self.inst.get().optimalTranslation()
        cdef list py_result = _r
        return py_result
    
    def loadModelAndData(self,  in_0 , list in_1 ,  in_2 , list in_3 ):
        """Cython signature: void loadModelAndData(int, libcpp_vector[POINT3D], int, libcpp_vector[POINT3D])"""
        assert isinstance(in_0, (int, long)), 'arg in_0 wrong type'
        assert isinstance(in_1, list) and all(isinstance(elemt_rec, POINT3D) for elemt_rec in in_1), 'arg in_1 wrong type'
        assert isinstance(in_2, (int, long)), 'arg in_2 wrong type'
        assert isinstance(in_3, list) and all(isinstance(elemt_rec, POINT3D) for elemt_rec in in_3), 'arg in_3 wrong type'
    
        cdef libcpp_vector[_POINT3D] * v1 = new libcpp_vector[_POINT3D]()
        cdef POINT3D item1
        for item1 in in_1:
            v1.push_back(deref(item1.inst.get()))
    
        cdef libcpp_vector[_POINT3D] * v3 = new libcpp_vector[_POINT3D]()
        cdef POINT3D item3
        for item3 in in_3:
            v3.push_back(deref(item3.inst.get()))
        self.inst.get().loadModelAndData((<int>in_0), deref(v1), (<int>in_2), deref(v3))
        del v3
        del v1
    
    def optimalRotation(self):
        """Cython signature: libcpp_vector[libcpp_vector[double]] optimalRotation()"""
        _r = self.inst.get().optimalRotation()
        cdef list py_result = _r
        return py_result
    
    def Register(self):
        """Cython signature: float Register()"""
        cdef float _r = self.inst.get().Register()
        py_result = <float>_r
        return py_result
    
    def BuildDT(self):
        """Cython signature: void BuildDT()"""
        self.inst.get().BuildDT()
    
    def setDTSizeAndFactor(self,  in_0 , double in_1 ):
        """Cython signature: void setDTSizeAndFactor(int, double)"""
        assert isinstance(in_0, (int, long)), 'arg in_0 wrong type'
        assert isinstance(in_1, float), 'arg in_1 wrong type'
    
    
        self.inst.get().setDTSizeAndFactor((<int>in_0), (<double>in_1))
    
    def __init__(self):
        """Cython signature: void GoICP()"""
        self.inst = shared_ptr[_GoICP](new _GoICP())
    
    def setInitNodeRot(self, ROTNODE in_0 ):
        """Cython signature: void setInitNodeRot(ROTNODE &)"""
        assert isinstance(in_0, ROTNODE), 'arg in_0 wrong type'
    
        self.inst.get().setInitNodeRot((deref(in_0.inst.get()))) 
