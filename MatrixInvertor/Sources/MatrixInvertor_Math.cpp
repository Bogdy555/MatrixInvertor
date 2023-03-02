#include "..\Headers\MatrixInvertor.hpp"



const float Determinant(const float* _Mat, const size_t _N)
{
	if (_N == 0)
	{
		return 0.0f;
	}

	if (_N == 1)
	{
		return *_Mat;
	}

	float _Rez = 0.0f;

	float* _SubMatrix = new float[(_N - 1) * (_N - 1)];

	for (size_t _X = 0; _X < _N; _X++)
	{
		size_t _SubI = 0;

		for (size_t _I = 1; _I < _N; _I++)
		{
			size_t _SubJ = 0;

			for (size_t _J = 0; _J < _N; _J++)
			{
				if (_J == _X)
				{
					continue;
				}

				_SubMatrix[_SubJ + _SubI * (_N - 1)] = _Mat[_J + _I * _N];

				_SubJ++;
			}

			_SubI++;
		}

		_Rez = _Rez + powf(-1.0f, (float)(_X)) * _Mat[_X + 0 * _N] * Determinant(_SubMatrix, _N - 1);
	}

	delete[] _SubMatrix;

	return _Rez;
}

const float DeterminantGauss(const float* _Mat, const size_t _N)
{
	if (_N == 0)
	{
		return 0.0f;
	}

	if (_N == 1)
	{
		return *_Mat;
	}

	float _DetRez = 1.0f;

	float* _AuxMat = new float[_N * _N];

	for (size_t _Y = 0; _Y < _N; _Y++)
	{
		for (size_t _X = 0; _X < _N; _X++)
		{
			_AuxMat[_X + _Y * _N] = _Mat[_X + _Y * _N];
		}
	}

	for (size_t _Line = 0; _Line < _N; _Line++)
	{
		{
			size_t _SwapLine = _Line;

			while (_AuxMat[_Line + _SwapLine * _N] == 0.0f)
			{
				_SwapLine++;

				if (_SwapLine >= _N)
				{
					delete[] _AuxMat;
					return 0.0f;
				}
			}

			if (_SwapLine != _Line)
			{
				_DetRez *= -1.0f;

				for (size_t _X = 0; _X < _N; _X++)
				{
					float _Aux = _AuxMat[_X + _Line * _N];
					_AuxMat[_X + _Line * _N] = _AuxMat[_X + _SwapLine * _N];
					_AuxMat[_X + _SwapLine * _N] = _Aux;
				}
			}
		}

		{
			float _Div = _AuxMat[_Line + _Line * _N];

			_DetRez *= _Div;

			for (size_t _X = 0; _X < _N; _X++)
			{
				_AuxMat[_X + _Line * _N] /= _Div;
			}
		}

		for (size_t _NextLine = _Line + 1; _NextLine < _N; _NextLine++)
		{
			float _Mul = _AuxMat[_Line + _NextLine * _N];

			for (size_t _X = 0; _X < _N; _X++)
			{
				_AuxMat[_X + _NextLine * _N] -= _Mul * _AuxMat[_X + _Line * _N];
			}
		}
	}

	delete[] _AuxMat;

	return _DetRez;
}

const bool Inverse(const float* _Src, float* _Dest, const size_t _N)
{
	const float _Det = Determinant(_Src, _N);

	if (!_Det)
	{
		return false;
	}

	if (_N == 1)
	{
		*_Dest = 1.0f / (*_Src);
		return true;
	}

	float _InvDet = 1.0f / _Det;

	float* _AuxMat = new float[(_N - 1) * (_N - 1)];

	for (size_t _Y = 0; _Y < _N; _Y++)
	{
		for (size_t _X = 0; _X < _N; _X++)
		{
			for (size_t _CopyY = 0; _CopyY < _Y; _CopyY++)
			{
				for (size_t _CopyX = 0; _CopyX < _X; _CopyX++)
				{
					_AuxMat[_CopyX + _CopyY * (_N - 1)] = _Src[_CopyY + _CopyX * _N];
				}
			}

			for (size_t _CopyY = 0; _CopyY < _Y; _CopyY++)
			{
				for (size_t _CopyX = _X; _CopyX < _N - 1; _CopyX++)
				{
					_AuxMat[_CopyX + _CopyY * (_N - 1)] = _Src[_CopyY + (_CopyX + 1) * _N];
				}
			}

			for (size_t _CopyY = _Y; _CopyY < _N - 1; _CopyY++)
			{
				for (size_t _CopyX = 0; _CopyX < _X; _CopyX++)
				{
					_AuxMat[_CopyX + _CopyY * (_N - 1)] = _Src[(_CopyY + 1) + _CopyX * _N];
				}
			}

			for (size_t _CopyY = _Y; _CopyY < _N - 1; _CopyY++)
			{
				for (size_t _CopyX = _X; _CopyX < _N - 1; _CopyX++)
				{
					_AuxMat[_CopyX + _CopyY * (_N - 1)] = _Src[(_CopyY + 1) + (_CopyX + 1) * _N];
				}
			}

			_Dest[_X + _Y * _N] = powf(-1.0f, (float)(_X + _Y)) * Determinant(_AuxMat, _N - 1) * _InvDet;
		}
	}

	delete[] _AuxMat;

	return true;
}

const bool InverseGauss(const float* _Src, float* _Dest, const size_t _N)
{
	if (_N == 1)
	{
		if (*_Src == 0.0f)
		{
			return false;
		}

		*_Dest = 1.0f / (*_Src);
		return true;
	}

	float* _AuxMat = new float[_N * _N];

	for (size_t _Y = 0; _Y < _N; _Y++)
	{
		for (size_t _X = 0; _X < _N; _X++)
		{
			_AuxMat[_X + _Y * _N] = _Src[_X + _Y * _N];
			_Dest[_X + _Y * _N] = (float)(_X == _Y);
		}
	}

	for (size_t _Line = 0; _Line < _N; _Line++)
	{
		{
			size_t _SwapLine = _Line;

			while (_AuxMat[_Line + _SwapLine * _N] == 0.0f)
			{
				_SwapLine++;

				if (_SwapLine >= _N)
				{
					delete[] _AuxMat;
					return false;
				}
			}

			if (_SwapLine != _Line)
			{
				for (size_t _X = 0; _X < _N; _X++)
				{
					{
						float _Aux = _AuxMat[_X + _Line * _N];
						_AuxMat[_X + _Line * _N] = _AuxMat[_X + _SwapLine * _N];
						_AuxMat[_X + _SwapLine * _N] = _Aux;
					}

					{
						float _Aux = _Dest[_X + _Line * _N];
						_Dest[_X + _Line * _N] = _Dest[_X + _SwapLine * _N];
						_Dest[_X + _SwapLine * _N] = _Aux;
					}
				}
			}
		}

		{
			float _Div = _AuxMat[_Line + _Line * _N];

			for (size_t _X = 0; _X < _N; _X++)
			{
				_AuxMat[_X + _Line * _N] /= _Div;
				_Dest[_X + _Line * _N] /= _Div;
			}
		}

		for (size_t _NextLine = _Line + 1; _NextLine < _N; _NextLine++)
		{
			float _Mul = _AuxMat[_Line + _NextLine * _N];

			for (size_t _X = 0; _X < _N; _X++)
			{
				_AuxMat[_X + _NextLine * _N] -= _Mul * _AuxMat[_X + _Line * _N];
				_Dest[_X + _NextLine * _N] -= _Mul * _Dest[_X + _Line * _N];
			}
		}
	}

	for (size_t _CurrentLine = 0; _CurrentLine < _N; _CurrentLine++)
	{
		for (size_t _TopLine = _CurrentLine + 1; _TopLine < _N; _TopLine++)
		{
			float _Mult = _AuxMat[(_N - 1 - _CurrentLine) + (_N - 1 - _TopLine) * _N];
			for (size_t _X = 0; _X < _N; _X++)
			{
				_AuxMat[_X + (_N - 1 - _TopLine) * _N] -= _Mult * _AuxMat[_X + (_N - 1 - _CurrentLine) * _N];
				_Dest[_X + (_N - 1 - _TopLine) * _N] -= _Mult * _Dest[_X + (_N - 1 - _CurrentLine) * _N];
			}
		}
	}

	delete[] _AuxMat;

	return true;
}

void Multiply(const float* _A, const float* _B, float* _Rez, const size_t _N)
{
	for (size_t _Index = 0; _Index < _N * _N; _Index++)
	{
		_Rez[_Index] = 0.0f;
	}

	for (size_t i = 0; i < _N; i++)
	{
		for (size_t j = 0; j < _N; j++)
		{
			for (size_t k = 0; k < _N; k++)
			{
				_Rez[i * _N + j] += _A[i * _N + k] * _B[k * _N + j];
			}
		}
	}
}

const bool CheckIdentity(const float* _I, const size_t _N)
{
	for (size_t _Y = 0; _Y < _N; _Y++)
	{
		for (size_t _X = 0; _X < _N; _X++)
		{
			float _ExpVal = 0.0f;
			if (_X == _Y)
			{
				_ExpVal = 1.0f;
			}
			if (abs(_I[_X + _Y * _N] - _ExpVal) > 0.01f)
			{
				return false;
			}
		}
	}

	return true;
}

void GetRandMat(float* _A, const size_t _N, const time_t _Seed)
{
	std::mt19937 _Generator((unsigned int)(_Seed));
	std::uniform_real_distribution<float> _Distribution(-1000.0f, 1000.0f);

	do
	{
		for (size_t _Y = 0; _Y < _N; _Y++)
		{
			for (size_t _X = 0; _X < _N; _X++)
			{
				_A[_X + _Y * _N] = _Distribution(_Generator);
			}
		}
	} while (!Determinant(_A, _N));
}
