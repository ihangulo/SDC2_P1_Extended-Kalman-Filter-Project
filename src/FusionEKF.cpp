#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iomanip>  // std::setprecision
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 
 https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/3612b91d-9c33-47ad-8067-a572a6c93837/concepts/252f0093-48ac-4122-aaae-f10214d30320
 above source // tracking.cpp is helpful.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
 
  //measurement covariance matrix - laser // given value
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar //given value
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

		// ---- given values end
  /**
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  //create a 4D state vector, we don't know yet the values of the x state
  VectorXd x_ = VectorXd(4);
  x_ = VectorXd(4); // will be input
  //x_ << 0.1, 0.1, 0.1, 0.1; // dummy value
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1,0,0,0,
        0,1,0,0;
  
		
  //state covariance matrix P
  MatrixXd P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;
		  
  //the initial transition matrix F_
  MatrixXd F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
    0, 1, 0, 1,
    0, 0, 1, 0,
    0, 0, 0, 1; // it will be modified at prediction

		
  MatrixXd Q_ = MatrixXd(4, 4);
  // Given value : Use noise_ax = 9 and noise_ay = 9 for your Q matrix. // given value

  noise_ax = 9;
  noise_ay = 9;
  

  //void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
   //                     MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {

  ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_); // x_, H_ , R_ are will be changed by type
  
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	  // instead of H, use Hj_
	  float ro = measurement_pack.raw_measurements_[0];
	  float phi = measurement_pack.raw_measurements_[1];
	  float ro_dot = measurement_pack.raw_measurements_[2];
	 
	  
      //set the state with the initial location and real velocity
	  ekf_.x_ << ro * cos(phi), ro * sin(phi), ro_dot* cos(phi), ro_dot*sin(phi); // real velocity
	  // ekf_.x_ << ro * cos(phi), ro * sin(phi), 0.0, 0.0; // start 0 velocity
	 
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	  float px = measurement_pack.raw_measurements_[0];
	  float py = measurement_pack.raw_measurements_[1];
	  
	  /**
	   In Data 2, the starting lidar measurements for x, y are both zero, and this special case can create problems
	   for both the EKF and UKF lidar update states, in particular
	   for the EKF when calculating the Jacobian. One way to catch for this is to observe 
	   when both px, py are zero and instead set them to some small floating value.
	   */
	   if (fabs(px) <0.0001 && fabs(py) <0.0001) { // px == 0 -->Nan Error occured 
		   px = (px<0)? -0.0001 : 0.0001; // remain sign of number
		   py = (py<0)? -0.0001 : 0.0001;
	   }
	   
	   //set the state with the initial location and zero velocity
       ekf_.x_ << px, py, 0.0, 0.0;

    }
	previous_timestamp_ = measurement_pack.timestamp_;
	
	// cout << "init previous_timestamep" << previous_timestamp_ << endl;
	
	// prohibit divide by zero error
	if (noise_ax ==0)
	  noise_ax = 0.0001;
    if(noise_ay == 0)
	  noise_ay= 0.0001;

    // done initializing, no need to predict or update
    is_initialized_ = true;

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
 //cout << "Prediction P_ = " << ekf_.P_ << endl;
  /**
   //TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   
   
  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0f;	//dt - expressed in seconds
  
//   cout << "now dt=" << measurement_pack.timestamp_ << endl << previous_timestamp_ << endl << dt<< endl; // debug
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  /** finally 
  ekf_.F_<< 1, 0, dt, 0,
            0, 1, 0, dt,
			0, 0, 1, 0,
			0, 0, 0, 1; 
  */
  //set the process covariance matrix Q

  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
			   0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
			   dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
			   0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();
  

  /*****************************************************************************
   *  Update
   ****************************************************************************/
 //cout << "Update P_ = " << ekf_.P_ << std::endl;
  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	ekf_.R_ = R_radar_;
	ekf_.H_ = tools.CalculateJacobian(ekf_.x_); // Hj
	ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	
  } else {
	  
    // Laser updates
    ekf_.R_ = R_laser_; // set type
	ekf_.H_ = H_laser_;
	ekf_.Update(measurement_pack.raw_measurements_);
	//ekf_.Update(xx_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
