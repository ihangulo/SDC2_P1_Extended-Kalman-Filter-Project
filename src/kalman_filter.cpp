#include "kalman_filter.h"
#include <iostream> // debug

// https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/3612b91d-9c33-47ad-8067-a572a6c93837/concepts/252f0093-48ac-4122-aaae-f10214d30320

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
  x_ = F_ * x_;

  MatrixXd Ft = F_.transpose();

  P_ = F_ * P_ * Ft + Q_;


}

void KalmanFilter::Update(const VectorXd &z) {
  /**
 	Note. copied frm quiz solution of udacity.
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;

  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_; // H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();

  MatrixXd K = PHt * Si;
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */

  // get values from matrix
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  /**
  * Before and while calculating the Jacobian matrix Hj,
  * make sure your code avoids dividing by zero. For example, both the x and y values might be zero
  * or px*px + py*py might be close to zero.
  */

   // avoid px value divide by zero // and Jacobian calculation
   if (fabs(px) <0.001 && fabs(py) <0.001) { // px == 0 -->Nan Error occured
		   px = (px<0)? -0.001 : 0.001; // remain sign of number
		   py = (py<0)? -0.001 : 0.001;
   }
   float ro_2 = px * px + py * py;
   if( ro_2 < 0.000001)
	   ro_2 = 0.000001;

   float ro = sqrt (ro_2); // min(row_2) = 0.000001 --> min(row) == 0.001 (cannot be zero)


   /**
   * In C++, atan2() returns values between -pi and pi. When calculating phi in y = z - h(x) for
   radar measurements, the resulting angle phi in the y vector should be adjusted so that
   it is between -pi and pi. The Kalman filter is expecting small angle values between the range -pi and pi.
   */

  float phi = atan2(py, px); //  = atan(py/px); // px==py==0 cannot be founded
                             // (if that case, it changed px=(+ or -) 0.001,  py = (+ or -) 0.001
  float ro_dot = (px*vx + py*vy) / ro; // ro cannot be zero

  VectorXd z_pred(3);
  z_pred << ro, phi, ro_dot;


  // z_pred is input value but initialize before.. (change before call)
  // same with Kalman Filter
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_; // H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();

  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}
