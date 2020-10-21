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
    /* x,y,z is initialized by gps x,y,z w stddev*/
    
    num_particles = 100;  // Set the number of particles
    
    std::default_random_engine gen;
    
    std::normal_distribution<double> dist_x(0., std[0]);
    std::normal_distribution<double> dist_y(0., std[1]);
    std::normal_distribution<double> dist_theta(0., std[2]);

    for (int i = 0; i < num_particles; ++i) {
        
      Particle p;
        
      p.id = i;
      p.x = x + dist_x(gen);
      p.y = y + dist_y(gen);
      p.theta = theta + dist_theta(gen);
      p.weight = 1.0;
        
      particles.push_back(p);
        
      weights.push_back(1.0);
    }
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
    /* With control signal ut: velocity, yaw_rate, delta_t,
       predict p(x_t) using Bicycle model.
       Namely, calc p(x_t| ut, x_t-1)   */
    
    std::default_random_engine gen;
    
    std::normal_distribution<double> noise_x(0., std_pos[0]);
    std::normal_distribution<double> noise_y(0., std_pos[1]);
    std::normal_distribution<double> noise_theta(0., std_pos[2]);
    
    for (int i = 0; i < num_particles; ++i) {
        double x0 = particles[i].x;
        double y0 = particles[i].y;
        double theta0 = particles[i].theta;
        double x_pred,y_pred,theta_pred;
        
        //Bicycle model as motion model,
        if (fabs(yaw_rate) < 0.0001){ //density prob.
            x_pred = x0 + velocity * delta_t * cos(theta0);
            y_pred = y0 + velocity * delta_t * sin(theta0);
            theta_pred = theta0;
        }
        else {
            x_pred = x0 + (velocity / yaw_rate) * (sin(theta0 + yaw_rate * delta_t) - sin(theta0));
            y_pred = y0 + (velocity / yaw_rate) * (cos(theta0) - cos(theta0 + yaw_rate * delta_t));
            theta_pred = theta0 + yaw_rate * delta_t;
        }
        particles[i].x = x_pred + noise_x(gen);
        particles[i].y = y_pred + noise_y(gen);
        particles[i].theta = theta_pred + noise_theta(gen);
    }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
    /* For each observation, find the nearest landmark.
       'predicted' is estimated ground truth.
       'observation' is sensor-measurred landmark estimation. */
    
    for (unsigned i = 0; i < observations.size(); ++i) {
        int nearest = observations[i].id;
        //double min_cost = std::numeric_limits<double>::max();
        double min_cost = std::numeric_limits<const float>::infinity();
        
        for (unsigned j = 0; j < predicted.size(); ++j){
            double cost = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
            if (cost < min_cost){
                min_cost = cost;
                nearest = predicted[j].id;
            }
        }
        observations[i].id = nearest;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
    /* For each particle, calc wt, or p(z_t | x_t) using Gaussian model
       p(z_t | x_t) ~ N(observation, map_landmark, std_landmark)
       As z_t is mult-dimensional(multi-landmarks), get joint prob. by product */
     
    for (int i = 0; i < num_particles; ++i) {
        
        //transform observation into world coordinates
        vector<LandmarkObs> w_observations;
        
        for (unsigned j = 0; j < observations.size(); ++j){
            double xp = particles[i].x;
            double yp = particles[i].y;
            double theta  = particles[i].theta;
            double x_obs  = observations[j].x;
            double y_obs  = observations[j].y;

            LandmarkObs landmk;
            
            landmk.x = xp + (cos(theta) * x_obs) - (sin(theta) * y_obs);
            landmk.y = yp + (sin(theta) * x_obs) + (cos(theta) * y_obs);
            
            w_observations.push_back(landmk);
        }
        
        //From a particle(a guess of vehicle position) in world coordinates,
        //to map-provided landmark position, measure distance
        //Only leave landmarks within range
        vector<LandmarkObs> landmarks_in_range;
        
        for (unsigned k = 0; k < map_landmarks.landmark_list.size(); ++k){
            double x_part = particles[i].x;
            double y_part = particles[i].y;
            int map_id = map_landmarks.landmark_list[k].id_i;
            double map_x = map_landmarks.landmark_list[k].x_f;
            double map_y = map_landmarks.landmark_list[k].y_f;
            
            double distance = dist(x_part, y_part, map_x, map_y);
            if (distance <= sensor_range){
                
                LandmarkObs landmark;
                
                landmark.id = map_id;
                landmark.x = map_x;
                landmark.y = map_y;
                
                landmarks_in_range.push_back(landmark);
            }
        }
        
        //get the nearest landmarks_in_range for each observation.
        //by filling w_observations.id as matching landmark.id
        dataAssociation(landmarks_in_range, w_observations);
        
        //For each particle, observations and map landmarks are paired.
        //it's probability is a joint prob., a product of each dimension Gaussian prob.
        double sig_x = std_landmark[0];
        double sig_y = std_landmark[1];
        double gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);
        double joint_prob = 1.0;
        
        for (unsigned l = 0; l < w_observations.size(); ++l){
            double obs_x = w_observations[l].x;
            double obs_y = w_observations[l].y;
            double obs_id = w_observations[l].id;
            
            for (unsigned m = 0; m < landmarks_in_range.size(); ++m){
                    double ldmk_x = landmarks_in_range[m].x;
                    double ldmk_y = landmarks_in_range[m].y;
                    double ldmk_id = landmarks_in_range[m].id;

                    if ((obs_id - ldmk_id) == 0){
                        
                        double exponent = (pow(obs_x - ldmk_x, 2) / (2 * pow(sig_x, 2))) + (pow(obs_y - ldmk_y, 2) / (2 * pow(sig_y, 2)));

                        joint_prob = joint_prob * gauss_norm * exp(-1. * exponent);
                        break;
                }
            }
        }
        particles[i].weight = joint_prob;
        weights[i] = joint_prob;
    }
    
    //update particle's weight into normalized joint prob.
    double normalizer = std::accumulate(weights.begin(), weights.end(), 0.0f);
    for (int i = 0; i < num_particles; ++i){
        particles[i].weight = particles[i].weight / normalizer;
        weights[i] = particles[i].weight;
    }
}

void ParticleFilter::resample() {
    /* std::discrete_distribution produces random integers on the interval[0,n],
     with given each number's weight. Refer to
     https://en.cppreference.com/w/cpp/numeric/random/discrete_distribution */

    std::default_random_engine gen;
    
    std::discrete_distribution<int> sample(weights.begin(), weights.end());
    
    //particles' number in invariant
    //whereas particles' diversity getting sparse
    //by resampling on weights
    vector<Particle> new_particles;
    for (int i = 0; i < num_particles; ++i)
    {
        int k = sample(gen);
        new_particles.push_back(particles[k]);
    }
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
