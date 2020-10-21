# **Particle Filter Project**


[//]: # (Image References)

[image1]: ./pictures/PF_architecture.png "function blocks"

[image2]: ./pictures/Bicycle.png "motion model"

[image3]: ./pictures/Transform.png "Coord_Transformation"

[image4]: ./pictures/Gaussian.png "Multivariate Gaussian model"

[image5]: ./pictures/Results.png "Simulation results"

[image6]: ./pictures/PF_algorithm.png "Algorithm"

---

## Overview

In this project, Particle filter was implemented to estimate the localization of moving vehicle. As input, it takes a (noisy) GPS estimate of its initial location, a map of landmark locations, and lots of (noisy) sensor and control data at each time step. Output is a posterior probability distribution of where the vehicle is.

Inside the system, each particle has a guess of vehicle's state(location & heading) with weight(probability). Each particle's guess is a time-varying state and updated by 'move and sense' routine of Bayesian. A set of particles represents the output posterior probability distribution of where the vehicle is. 

![alt text][image1]

Steps:

- To initialise, N particles were created and set by GPS pose information with Gaussian noise.

- Make a prediction for each particle on the vehicles pose assuming a bicycle (motion) model.

- Update the weights of each particle based on the distance from map landmark to measured landmark positions,
  assuming Gaussian model.

- Resample N particles assigning a higher likeihood of picking particles with a higher weight.

- Repeat the process.


---

## Implementing particle filter

My project includes the following files:
* src/main.cpp uses for uWebSocketIO in communicating with the simulator.
* src/particle_filter.cpp, src/particle_filter.h contains localization process & measurement mgmt.
* src/map.h contains map data format
* src/helper_functions.h contains control/landmark data format and distance calculation routine, etc.,


![alt text][image6]


### Initialization
GPS is inaccurate because of the noise with 1m error at best and signal reception can be blocked in-door as we want centi-meteter level accuray even in in-door environment. Still, if available, initializing our particle filter using GPS often can help localization process initially converge faster.  

Here, a particle's initial guess on 2d position(x,y) and heading(theta) was initialized by Gaussian distribution using gps x,y,theta as mean values. The initial number of particles to 100. Since the posterior probability distribution is dependent on the number of particles, the number of particles needs to be set higher enough to keep our posterior distribution from sparsity problem.

### (Motion) Prediction
With control signal ut(velocity, yaw_rate, delta_t) and previous state p(x_t-1), we predict p(x_t) using Bicycle model. Namely, we calculate p(x_t| ut, x_t-1). Since the yaw_rate can be close to zero(not exactly zero in continous space), we seperate the case to prevent devide-by-zero errors. Evolving states using motion prediction model also involves Gaussian noise every time step, adding errors, contributing to uncertainty.

![alt text][image2]

### Weight Update
For each particle, we calculate its weight wt, or p(z_t | x_t) using Gaussian model. Namely, p(z_t | x_t) ~ N(observation, map_landmark, std_landmark). As z_t is mult-dimensional(multi-landmarks), we get joint probability by product.

First, we transformed the observations in the vehicle coordinate system into the world coordinate
To do this we applied the homogeneous matrix transformation:

![alt text][image3]

Then, we measured the distance from a particle(a guess of vehicle position in world coordinates) to landmark position on the map, leaving only landmarks within range.
        
After that, for each observation we searched the nearest map landmark using the nearest neighbor method, matching 1:1.
For a match, we filled the observation with the same landmark's id as the map landmark's id.
        
For each particle, it's observations and corresponding map landmarks are paired. Since the particles's weight is a joint probability of all dimensions, we took a product of Gaussian probability of each dimension.

![alt text][image4]

Finally, we updated a particle's weight into normalized joint probability by dividing it by the sum of probablities.
       
### Resampling
We used std::discrete_distribution produces random integers on the interval[0,n], with given each number's weight. 
Here, we set the number of particles constant 100, and resampled allowing replacement. Namely, the number of particles is invariant whereas the diversity of particles can become sparse because particles with higher weights have higher chance to be resampled. This sparsity issue could be a potential problem. We'll discuss more on this issue in the future work.


## Running the Code
Once the install for uWebSocketIO is complete, the main program can be built and ran by doing the following from the project top directory.
```sh
1. ./clean.sh
2. ./build.sh
3. ./run.sh
```

![alt text][image5]

This repository includes two files that can be used to set up and install [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac systems. For windows you can use either Docker, VMware, or even [Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) to install uWebSocketIO. 

To run the code,
- run shell script to install [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for Linux, Mac, or Ubuntu systems.
- check dependencies
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
- install Udacity simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

Tips for setting up your environment can be found [here](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/23d376c7-0195-4276-bdf0-e02f1f3c665d)


