#include "ukf.h"

using namespace std;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	// initially set to false, set to true in first call of ProcessMeasurement
	is_initialized_ = false;
	
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// state dimention
	n_x_ = 5;

	// initial state vector
	x_ = VectorXd(n_x_);

	// initial covariance matrix
	P_ = MatrixXd(n_x_, n_x_);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 0.5;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.5;

	/**
	* DO NOT MODIFY measurement noise values below.
	* These are provided by the sensor manufacturer.
	*/

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

	/**
	* End DO NOT MODIFY section for measurement noise values
	*/

	// Augmented state dimension
	n_aug_ = n_x_ + 2;

	// Number of Sigma points
	n_sigma_points_ = 2 * n_aug_ + 1;

	// Sigma point spreading parameter
	lambda_ = 3 - n_aug_;

	// Lambda weights
	weights_ = VectorXd(n_sigma_points_);
	weights_.fill(0.5 / (n_aug_ + lambda_));
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	
	//State parameters for Laser and Radar measurements
	n_lidar_ = 2;
	n_radar_ = 3;

	// Predicted sigma points matrix
	Xsig_pred_ = MatrixXd(n_x_, n_sigma_points_);

	// NIS
    //float NIS_lidar_ = 0.0;
    //float NIS_radar_ = 0.0;

	/** Initialize the UKF */
	x_ << 0,0,0,0,0;
	P_ << 1,0,0,0,0,
			0,1,0,0,0,
			0,0,1,0,0,
			0,0,0,1,0,
			0,0,0,0,1;
			
	// Initialize measurement noice covarieance matrix
	R_radar_ = MatrixXd(3, 3);
	R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;

	R_lidar_ = MatrixXd(2, 2);
	R_lidar_ << std_laspx_*std_laspx_,0,
              0,std_laspy_*std_laspy_;

}

UKF::~UKF() {}

void UKF::NormalizeAngle(double& phi)
{
	while (phi > M_PI) phi -= 2. * M_PI;
	while (phi < -M_PI) phi += 2. * M_PI;
}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
	if(!is_initialized_)
	{
		if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
		{
			double px = meas_package.raw_measurements_[0];
			double py = meas_package.raw_measurements_[1];

			x_ << px, py, 0.0, 0.0, 0.0;

		}

		else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
		{
			double rho = meas_package.raw_measurements_[0];
			double psi = meas_package.raw_measurements_[1];
			double rho_dot = meas_package.raw_measurements_[2];
			
			double px = rho * cos(psi);
			double py = rho * sin(psi);

			x_ << px, py, rho_dot, 0.0, 0.0;

		}
		// Initialize the state covariance matrix
		P_ = MatrixXd::Identity(n_x_, n_x_);

		//Update last measurement
		previous_timestamp_ = meas_package.timestamp_;
		is_initialized_ = true;
		
		return;

	}
	//Prediction
	double dt = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;

	// update last measurement time
	previous_timestamp_ = meas_package.timestamp_;
	
	Prediction(dt);
	
	// Update measurement
	if(meas_package.sensor_type_ == MeasurementPackage::LASER)
	{
		UpdateLidar(meas_package);
	}

	if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
		UpdateRadar(meas_package);
	}

}

MatrixXd UKF::GetSigmaPoints(double dt)
{
	// Get Sigma points from current state and covariance matrix
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;
	
	//Augmented state covariance
	MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
	MatrixXd Q = MatrixXd(2,2);
	Q << std_a_*std_a_, 0,
			0, std_yawdd_*std_yawdd_;

	P_aug.topLeftCorner(5, 5) = P_;
	P_aug.bottomRightCorner(2, 2) = Q;

	MatrixXd A = P_aug.llt().matrixL();
	
	// Compute sigma point matrix
	MatrixXd sigma_pts = MatrixXd::Zero(n_aug_, n_sigma_points_);
	sigma_pts.col(0) = x_aug;
	
	for(int i=0;i < n_aug_; i++)
	{
		sigma_pts.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
		sigma_pts.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
	}
	
	// Sigma point prediction from estimates sigma points and dt
	for (int i = 0; i < n_sigma_points_; i++) 
	{

		// Auxiliary variables for readability
		double p_x		= sigma_pts(0, i);
		double p_y		= sigma_pts(1, i);
		double v		= sigma_pts(2, i);
		double yaw		= sigma_pts(3, i);
		double yawd		= sigma_pts(4, i);
		double nu_a		= sigma_pts(5, i);
		double nu_yawdd	= sigma_pts(6, i);

		// Sanity check
		if (fabs(p_x) < 0.001 && fabs(p_y) < 0.001) {
			p_x = 0.1;
			p_y = 0.1;
		}

		// Predicted state values
		double px_p, py_p;
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd * dt) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * dt));
		}
		else {
			px_p = p_x + v * dt * cos(yaw);
			py_p = p_y + v * dt * sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd * dt;
		double yawd_p = yawd;

		// Handle noise
		px_p = px_p + 0.5 * nu_a * dt * dt * cos(yaw);
		py_p = py_p + 0.5 * nu_a * dt * dt * sin(yaw);
		v_p = v_p + nu_a * dt;
		yaw_p = yaw_p + 0.5 * nu_yawdd * dt * dt;
		yawd_p = yawd_p + nu_yawdd * dt;

		// Fill current column of Xsig_pred matrix with sigma point just computed
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}
	return Xsig_pred_;
}


void UKF::Prediction(double delta_t)
{
	// Get predicted sigma points
	MatrixXd pred_sigma_x_ = GetSigmaPoints(delta_t);

	// Calculate Predicted mean
	VectorXd x = VectorXd::Zero(n_x_);

	for(int i=0; i<n_sigma_points_;i++)
	{
		x = x + weights_(i) * pred_sigma_x_.col(i);
	}

	// Calculate predicted covariance
	MatrixXd P = MatrixXd::Zero(n_x_,n_x_);

	for(int i=0; i<n_sigma_points_;i++)
	{
		VectorXd x_diff = pred_sigma_x_.col(i) - x;
		NormalizeAngle(x_diff(3));

		P = P + weights_(i) * x_diff * x_diff.transpose();

	}
	
	//cout<<"Predicted x_ = "<<x<<endl;
	x_ = x;
	P_ = P;

}

void UKF::UpdateLidar(MeasurementPackage meas_package)
{
	/** Predict measurement */
	cout<<"UpdateLidar"<<endl;
	// Project sigma points onto meas space
	MatrixXd Zsig = Xsig_pred_.block(0, 0, n_lidar_, n_sigma_points_);

	// Predicted measurement mean
	VectorXd z_pred = VectorXd::Zero(n_lidar_);

	for(int i=0; i<n_sigma_points_; i++)
	{
		z_pred = z_pred + weights_(i) * Zsig.col(i);

	}
	// Predicted  measurement covariance
	MatrixXd S = MatrixXd::Zero(n_lidar_,n_lidar_);
	for(int i=0;i<n_sigma_points_; i++)
	{
		VectorXd z_sig_diff = Zsig.col(i) - z_pred;
		
		S = S + weights_(i)* z_sig_diff * z_sig_diff.transpose();
		
	}
	// Lidar Noise
	S = S + R_lidar_;

	/** Update state */
	MatrixXd T = MatrixXd::Zero(n_x_, n_lidar_);
	for (int i=0;i<n_sigma_points_;i++)
	{
		// Residual
		VectorXd z_sig_diff = Zsig.col(i) - z_pred;
		//Normalize angle
		NormalizeAngle(z_sig_diff(1));
		
		//State difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//NormalizeAngle
		NormalizeAngle(x_diff[3]);
				
		T = T + weights_(i) * x_diff * z_sig_diff.transpose();

	}

	// Kalman Gain
	MatrixXd K = T * S.inverse();

	// state update
	VectorXd z = VectorXd(n_lidar_);
	z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

	VectorXd z_diff = z - z_pred;
	
	x_ = x_ + K * (z_diff);

	// Coraiance mattrix update
	P_ = P_ - K * S * K.transpose();
	// Compute NIS for laser sensor
	NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;

}

void UKF::UpdateRadar(MeasurementPackage meas_package)
{
	/** Predict measurement */
	cout<<"updateRadar"<<endl;
	// Project sigma points onto meas space
	MatrixXd Zsig = MatrixXd(n_radar_, n_sigma_points_);

	for(int i=0; i<n_sigma_points_;i++)
	{
		float px = Xsig_pred_(0,i);
		float py = Xsig_pred_(1,i);
		float v = Xsig_pred_(2, i);
		float phi = Xsig_pred_(3, i);

		float v1 = cos(phi) * v;
		float v2 = sin(phi) * v;

		Zsig(0,i) = sqrt(px*px + py*py);
		if (fabs(px) > 0.001)
		{
			Zsig(1,i) = atan2(py,px);
		}
		else
			Zsig(1,i) = M_PI/2;

		if (Zsig(0,i) > 0.001)
			Zsig(2,i) = (px * v1 + py * v2) / sqrt(px * px + py * py);
		else
			Zsig(2,i) = 0.0;
	}
	// Predicted measurement mean
	VectorXd z_pred = VectorXd::Zero(n_radar_);

	for(int i=0; i<n_sigma_points_; i++)
	{
		z_pred = z_pred + weights_(i) * Zsig.col(i);

	}
	// Predicted  measurement covariance
	MatrixXd S = MatrixXd::Zero(n_radar_,n_radar_);
	
	for(int i=0;i<n_sigma_points_; i++)
	{
		VectorXd z_sig_diff = Zsig.col(i) - z_pred;
		NormalizeAngle(z_sig_diff[1]);

		S = S + weights_(i) * z_sig_diff * z_sig_diff.transpose();
	}
	// Radar Noise
	S = S + R_radar_;
	/** Update state */
	MatrixXd T = MatrixXd::Zero(n_x_, n_radar_);
	for (int i=0;i<n_sigma_points_;i++)
	{
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		NormalizeAngle(x_diff[3]);

		VectorXd z_sig_diff = Zsig.col(i) - z_pred;
		NormalizeAngle(z_sig_diff[1]);

		T = T + weights_(i) * x_diff * z_sig_diff.transpose();

	}

	// Kalman Gain
	MatrixXd K = T * S.inverse();

	// state update
	VectorXd z = VectorXd(n_radar_);
	z = meas_package.raw_measurements_;

	VectorXd z_diff = z - z_pred;
	NormalizeAngle(z_diff[1]);
	x_ = x_ + K * (z_diff);

	// Coraiance mattrix update
	P_ = P_ - K * S * K.transpose();
	// Compute NIS for radar sensor
	NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
