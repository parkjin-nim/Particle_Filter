/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
    num_particles = 1000;  // Set the number of particles
    particles(num_particles);
    
    std::default_random_engine gen;
    // creates a normal (Gaussian) distribution for x
    normal_distribution<double> dist_x(gps_x, std[0]);
    normal_distribution<double> dist_y(gps_y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    for (int i = 0; i < num_particles; ++i) {
      // Sample from these normal distributions like this:
      particles[i].x = dist_x(gen);
      particles[i].y = dist_y(gen);
      particles[i].theta = dist_theta(gen);
      particles[i].weight = 1.;
    }

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    std::default_random_engine gen;
    normal_distribution<double> dist_x(particles.x, std_pos[0]);
    normal_distribution<double> dist_y(particles.y, std_pos[1]);
    normal_distribution<double> dist_theta(particles.theta, std_pos[2]);
    
    for (int i = 0; i < num_particles; ++i) {
        particles[i].x = particles[i].x + (velocity/yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
        particles[i].y = particles[i].y + (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
        particles[i].theta = particles[i].theta + yaw_rate * delta_t;
        
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
    vector<LandmarkObs> closest_to_observed;
    float min_cost;
    struct LandmarkObs nearest;
    
    // leave only nearest observation to predicted
    for (int i = 0; i < predicted.size(); ++i) {
        min_cost = std::numeric_limits<const float>::infinity();
        nearest = observations[0];
        
        // find min-cost observation
        for (int j = 0; j < observations.size(); ++j){
            cost = sqrt((predicted[i].x - observations[j].x)**2 + (predicted[i].y - observations[j].y)**2);
            if (cost < min_cost){
                min_cost = cost;
                nearest = observations[j];
            }
        }
        // remove all except nearest observation
        closest_to_observed.push_back(nearest);
    }
    // remove all except nearest observation
    observations.clear();
    observations = closest_to_observed;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
    for (int i = 0; i < num_particles; ++i) {
        for (int j = 0; j < observations.size(); ++j){
            // vehicle coordinates transformation
            double x_part, y_part, x_obs, y_obs, theta;
            x_part = particles[i].x;
            y_part = particles[i].y;
            x_obs  = observations[j].x;
            y_obs  = observations[j].y;
            theta  = observations[j].theta; // -90 degrees

            // transform to map x,y coordinate
            obsevations[j].x = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);
            obsevations[j].y = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);
        }
        
        // get the one nearest observation to map_landmark
        dataAssociation(map, observations);

        // calculate weight for a particle
        double gauss_norm = 1 / (2 * M_PI * std_landmark[0] * sigma_landmark[1]);
        double exponent 0.0f;
        double weight = 1.0f;
        for (int k = 0; k < observations.size(); ++k){
            // calculate exponent
            exponent = (pow(observations[k] - map_landmarks.x, 2) / (2 * pow(std_landmark[0], 2)))
            + (pow(observations[k] - map_landmarks.y, 2) / (2 * pow(sigma_landmark[1], 2)));

            // calculate weight using normalization terms and exponent
            weight *= gauss_norm * exp(-exponent);

        }
        particles[i].weight = weight;
    }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    std::default_random_engine gen;
    std::uniform_int_distribution<int> dist_idx(0, num_particles);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::vector<Particle> new_particles;
    vector<double> w;
    
    int idx = dist_idx(gen);
    float beta = 0.0f;
    for (i = 0; i < num_particles; ++i){
        w[i] = particles[i].weight;
    }
    double mw = max(w);
    
    for (i = 0; i < num_particles; ++i){
        beta += distribution(gen) * 2.0f * mw;
        while beta > w[index]:
            beta -= w[index];
            index = (index + 1) % num_particles;
        new_particles.push_back(particles[index])
    }
    
    particles.clear();
    particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
