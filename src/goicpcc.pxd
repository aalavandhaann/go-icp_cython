from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "jly_goicp.hpp":
    cdef cppclass POINT3D:
        float x
        float y
        float z
        POINT3D()
        POINT3D(float, float, float)
        void pointToString()
    
    cdef cppclass ROTNODE:
        float a
        float b
        float c
        float w
        float ub
        float lb
        int l
        ROTNODE()
    
    cdef cppclass TRANSNODE:
        float x
        float y
        float z 
        float w
        float ub
        float lb
        TRANSNODE()
        
    cdef cppclass GoICP:
        int Nm
        int Nd
        float MSEThresh
        float SSEThresh
        float icpThresh
        float trimFraction
        int inlierNum
        bool doTrim
    
        GoICP()
        float Register() except +
        void BuildDT() except + 
        libcpp_vector[libcpp_vector[double]] optimalRotation() except +
        libcpp_vector[double] optimalTranslation() except +
        
        void loadModelAndData(int, libcpp_vector[POINT3D], int, libcpp_vector[POINT3D]) except +
        void setInitNodeRot(ROTNODE&) except +
        void setInitNodeTrans(TRANSNODE&) except +
        void setDTSizeAndFactor(int, double) except +