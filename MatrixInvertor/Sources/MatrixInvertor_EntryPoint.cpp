#include "..\Headers\MatrixInvertor.hpp"



#define MAX_SIZE 9

#define START_SIZE 1

#define TESTS_COUNT 10



void TestInverse(InverseFnc _Invertor, const wchar_t* _FileName, time_t _Seed)
{
	std::wofstream _FOut(_FileName);

	for (size_t _Size = START_SIZE; _Size < MAX_SIZE + 1; _Size++)
	{
		std::vector<float> _Times;

		_FOut << "Size: " << _Size << L'\n';

		float* _A = new float[_Size * _Size];
		float* _Inv = new float[_Size * _Size];
		float* _Mult = new float[_Size * _Size];

		for (size_t _Test = 0; _Test < TESTS_COUNT; _Test++)
		{
			GetRandMat(_A, _Size, _Seed);
			_Seed++;

			Timer _TestTimer;

			_TestTimer.Start();

			bool _InvertorRez = _Invertor(_A, _Inv, _Size);

			_TestTimer.Stop();

			_Times.push_back(_TestTimer * 1000000.0f);
		}

		delete[] _Mult;
		delete[] _Inv;
		delete[] _A;

		float _Average = 0.0f;

		for (size_t _Index = 0; _Index < _Times.size(); _Index++)
		{
			_Average += _Times[_Index];
		}

		_Average /= (float)(TESTS_COUNT);

		_FOut << _Average << L'\n';
	}

	_FOut.close();
}



void TestDeterminant(DeterminantFnc _DetFnc, const wchar_t* _FileName, time_t _Seed)
{
	std::wofstream _FOut(_FileName);

	for (size_t _Size = START_SIZE; _Size < MAX_SIZE + 1; _Size++)
	{
		std::vector<float> _Times;

		_FOut << "Size: " << _Size << L'\n';

		float* _A = new float[_Size * _Size];

		for (size_t _Test = 0; _Test < TESTS_COUNT; _Test++)
		{
			GetRandMat(_A, _Size, _Seed);
			_Seed++;

			Timer _TestTimer;

			_TestTimer.Start();

			float _DetRez = _DetFnc(_A, _Size);

			_TestTimer.Stop();

			_Times.push_back(_TestTimer * 1000000.0f);
		}

		delete[] _A;

		float _Average = 0.0f;

		for (size_t _Index = 0; _Index < _Times.size(); _Index++)
		{
			_Average += _Times[_Index];
		}

		_Average /= (float)(TESTS_COUNT);

		_FOut << _Average << L'\n';
	}

	_FOut.close();
}



int main()
{
	time_t _Seed = time(nullptr);

	//std::wcout << L"Inverse\n";
	TestInverse(Inverse, L".\\Test Inverse.txt", _Seed);
	//std::wcout << L"InverseGauss\n";
	TestInverse(InverseGauss, L".\\Test InverseGauss.txt", _Seed);
	//std::wcout << L"Determinant\n";
	TestDeterminant(Determinant, L".\\Test Determinant.txt", _Seed);
	//std::wcout << L"DeterminantGauss\n";
	TestDeterminant(DeterminantGauss, L".\\Test DeterminantGauss.txt", _Seed);

	return 0;
}
