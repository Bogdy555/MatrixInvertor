#include "..\Headers\MatrixInvertor.hpp"



Timer::Timer() : Begin(), End()
{
	Begin = std::chrono::system_clock::now();
	End = Begin;
}

Timer::Timer(const Timer& _Other) : Begin(_Other.Begin), End(_Other.End)
{

}

Timer::Timer(Timer&& _Other) noexcept : Begin(std::move(_Other.Begin)), End(std::move(_Other.End))
{

}

Timer::~Timer()
{

}

void Timer::Start()
{
	Begin = std::chrono::system_clock::now();
}

void Timer::Stop()
{
	End = std::chrono::system_clock::now();
}

void Timer::Restart()
{
	Begin = std::chrono::system_clock::now();
	End = Begin;
}

Timer::operator const float() const
{
	return ((std::chrono::duration<float>)(End - Begin)).count();
}

void Timer::operator= (const Timer& _Other)
{
	Begin = _Other.Begin;
	End = _Other.End;
}

void Timer::operator= (Timer&& _Other) noexcept
{
	Begin = std::move(_Other.Begin);
	End = std::move(_Other.End);
}
