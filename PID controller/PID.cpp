#include "PID.h"
#include <limits>
using namespace std;

/*
* TODO: Complete the PID class.
*/

PID::PID() {}

PID::~PID() {}

void PID::Init(double Kp, double Ki, double Kd)  {
  this->Kp = Kp; // proportional coefficient
  this->Ki = Ki; // integral coefficient
  this->Kd = Kd; // differential coefficient

  i_error = 0.0; // integral cte
  p_error = numeric_limits<double>::max(); // previous cte
  d_error = 0.0; // differential cte
}

void PID::UpdateError(double cte) {
  // first time through make the previous error equal to cte
  if (p_error == numeric_limits<double>::max())
    p_error = cte;

  d_error = cte - p_error;
  p_error = cte;
  i_error += cte;
}

double PID::TotalError() {
  return -Kp * p_error - Kd * d_error - Ki * i_error;
}

