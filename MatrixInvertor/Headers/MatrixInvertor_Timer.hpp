#ifndef MatrixInvertor_Timer_hpp

#define MatrixInvertor_Timer_hpp



#include "MatrixInvertor.hpp"



class Timer
{

public:

	Timer();
	Timer(const Timer& _Other);
	Timer(Timer&& _Other) noexcept;
	~Timer();

	void Start();
	void Stop();
	void Restart();

	operator const float() const;

	void operator= (const Timer& _Other);
	void operator= (Timer&& _Other) noexcept;

private:

	std::chrono::time_point<std::chrono::system_clock> Begin;
	std::chrono::time_point<std::chrono::system_clock> End;

};



#endif
