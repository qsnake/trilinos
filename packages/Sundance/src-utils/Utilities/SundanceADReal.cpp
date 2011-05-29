#include "SundanceADReal.hpp"
#include "SundanceExceptions.hpp"


using namespace Sundance;
using namespace Teuchos;

ADReal ADReal::operator-() const
{
	ADReal rtn = *this;
	rtn.value_ = -rtn.value_;
  rtn.gradient_ = -gradient_;

  addFlops(1 + gradient_.dim());
	return rtn;
}

ADReal ADReal::operator+(const ADReal& other) const 
{
	ADReal rtn(*this);
	rtn += other;
	return rtn;
}

ADReal ADReal::operator-(const ADReal& other) const 
{
	ADReal rtn(*this);
	rtn -= other;
	return rtn;
}

ADReal ADReal::operator*(const ADReal& other) const 
{
	ADReal rtn(*this);
	rtn *= other;
	return rtn;
}

ADReal ADReal::operator/(const ADReal& other) const 
{
	ADReal rtn(*this);
	rtn /= other;
	return rtn;
}

ADReal& ADReal::operator+=(const ADReal& other) 
{
	value_ += other.value_;
	gradient_ += other.gradient_;

  addFlops(1 + gradient_.dim());

	return *this;
}

ADReal& ADReal::operator-=(const ADReal& other) 
{
	value_ -= other.value_;
	gradient_ -= other.gradient_;

  addFlops(1 + gradient_.dim());

	return *this;
}

ADReal& ADReal::operator*=(const ADReal& other) 
{
	gradient_ = (other.value_*gradient_ + value_*other.gradient_);
	value_ *= other.value_;

  addFlops(1 + 3*gradient_.dim());

	return *this;
}

ADReal& ADReal::operator/=(const ADReal& other) 
{
  TEST_FOR_EXCEPTION(other.value_ == 0.0, RuntimeError,
                     "ADReal::operator/=() division by zero");

	gradient_ = (gradient_/other.value_ 
							 - value_/(other.value_*other.value_) * other.gradient_);
	value_ /= other.value_;

  addFlops(3 + 3*gradient_.dim());

	return *this;
}


// operations with constant reals 

ADReal& ADReal::operator+=(const double& other) 
{
	value_ += other;
  addFlops(1);
	return *this;
}

ADReal& ADReal::operator-=(const double& other) 
{
	value_ -= other;

  addFlops(1);
	return *this;
}

ADReal& ADReal::operator*=(const double& other) 
{
	gradient_ *= other;
	value_ *= other;

  addFlops(1 + gradient_.dim());

	return *this;
}

ADReal& ADReal::operator/=(const double& other) 
{
  TEST_FOR_EXCEPTION(other == 0.0, RuntimeError,
                     "ADReal::operator/=() division by zero");

  addFlops(2 + gradient_.dim());
	gradient_ *= 1.0/other;
	value_ /= other;
	return *this;
}


ADReal ADReal::operator+(const double& other) const 
{
	ADReal rtn(*this);
	rtn += other;
	return rtn;
}

ADReal ADReal::operator-(const double& other) const 
{
	ADReal rtn(*this);
	rtn -= other;
	return rtn;
}

ADReal ADReal::operator*(const double& other) const 
{
	ADReal rtn(*this);
	rtn *= other;
	return rtn;
}

ADReal ADReal::operator/(const double& other) const 
{
	ADReal rtn(*this);
	rtn /= other;
	return rtn;
}


namespace Sundance
{
ADReal operator+(const double& scalar, const ADReal& a)
{
  return a+scalar;
}

ADReal operator-(const double& scalar, const ADReal& a)
{
  return -a+scalar;
}

ADReal operator*(const double& scalar, const ADReal& a)
{
  return a*scalar;
}

ADReal operator/(const double& scalar, const ADReal& a)
{
  ADReal rtn(a);
  rtn.reciprocate();
  return rtn*scalar;
}
}

void ADReal::reciprocate() 
{
	TEST_FOR_EXCEPTION(value_==0.0, RuntimeError,
                     "div by 0 in ADReal::reciprocate()");

  addFlops(1 + gradient_.dim());

	gradient_ /= value_;
	value_ = 1.0/value_;
}


	





