## Introduction
This projects implements [Unscented Kalman Filter](https://www.seas.harvard.edu/courses/cs281/papers/unscented.pdf) (UKF) in C++, using information from Lidar and Radar noisy sensors.

- Lidar sensor measures the position of the object in cartesian co-ordinates '(x, y)'
- Radar sensor identifies the location in polar co-ordinates and also measures the radial speed '(rho, phi, rho_dot)'

UKF uses **constant turn/yaw rate and velocity magnitude model** (CRTV) for this particular system copared to [Extended Kalman Filter's](https://github.com/ashsne/Extended-Kalman-Filter.git) (EKF) **constant velocity** model. UKF model handles non-linear motion better than EKF in essence UKF is more precise than EKF.

This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases).

This repository includes two files that can be used to set up and intall [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac systems. For windows you can use either Docker, VMware, or even [Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) to install uWebSocketIO. 

INPUT: values provided by the simulator to the c++ program

["sensor_measurement"] => the measurment that the simulator observed (either lidar or radar)


OUTPUT: values provided by the c++ program to the simulator

["estimate_x"] <= kalman filter estimated position x
["estimate_y"] <= kalman filter estimated position y
["rmse_x"]
["rmse_y"]
["rmse_vx"]
["rmse_vy"]

---

## Results
Results of the code can be visualized as below

![Dataset1](https://github.com/ashsne/Unscented-Kalman-Filter/blob/master/img/Dataset1.PNG)

## Other Important Dependencies
* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF`
