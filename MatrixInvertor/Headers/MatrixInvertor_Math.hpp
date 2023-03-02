#ifndef MatrixInvertor_Math_hpp

#define MatrixInvertor_Math_hpp



#include "MatrixInvertor.hpp"



typedef const bool (*InverseFnc)(const float* _Src, float* _Dest, const size_t _N);
typedef const float (*DeterminantFnc)(const float* _Mat, const size_t _N);

const float Determinant(const float* _Mat, const size_t _N);
const float DeterminantGauss(const float* _Mat, const size_t _N);
const bool Inverse(const float* _Src, float* _Dest, const size_t _N);
const bool InverseGauss(const float* _Src, float* _Dest, const size_t _N);
void Multiply(const float* _A, const float* _B, float* _Rez, const size_t _N);
const bool CheckIdentity(const float* _I, const size_t _N);
void GetRandMat(float* _A, const size_t _N, const time_t _Seed);



#endif
