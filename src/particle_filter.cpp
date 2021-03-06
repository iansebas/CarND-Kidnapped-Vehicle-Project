/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <cmath>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).



	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[3]);


	this->particles.clear();
	this->weights.clear();
	#pragma omp parallel for
	for (int i = 0; i < this->num_particles; ++i){
		Particle one_particle;
		one_particle.id = i;
		one_particle.x = dist_x(gen);
		one_particle.y = dist_y(gen);
		one_particle.theta = dist_theta(gen);
		one_particle.weight = 1;
		this->weights.push_back(one_particle.weight);
		this->particles.push_back(one_particle);

	}

	this->is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[3]);
	#pragma omp parallel for
	for (int i = 0; i < this->num_particles; ++i){
	    if (fabs(yaw_rate) > 0.000000001) {
	        this->particles[i].x += velocity/yaw_rate * ( sin (this->particles[i].theta + yaw_rate*delta_t) - sin(this->particles[i].theta));
	        this->particles[i].y += + velocity/yaw_rate * ( cos(this->particles[i].theta) - cos(this->particles[i].theta+yaw_rate*delta_t) );
	        this->particles[i].theta += + yaw_rate*delta_t;
	    }
	    else {
	        this->particles[i].x += velocity*delta_t*cos(this->particles[i].theta);
	        this->particles[i].y += velocity*delta_t*sin(this->particles[i].theta);
	    }
	    this->particles[i].x += dist_x(gen);
	    this->particles[i].y += dist_y(gen);
	    this->particles[i].theta += dist_theta(gen);
	}

}

LandmarkObs global_to_local(Particle p, Map::single_landmark_s landmark){
	LandmarkObs one_landmark;
	one_landmark.id = landmark.id_i;
	//Translate
	landmark.x_f -= p.x;
	landmark.y_f -= p.y;
	//Rotate
	one_landmark.x = cos(p.theta)*landmark.x_f + sin(p.theta)*landmark.y_f;
	one_landmark.y = -sin(p.theta)*landmark.x_f + cos(p.theta)*landmark.y_f;
	return one_landmark;
}

Map::single_landmark_s local_to_global(const Particle &p, const LandmarkObs &obs){
	Map::single_landmark_s global_obs;
	global_obs.id_i = obs.id;
	global_obs.x_f = p.x + cos(p.theta)*obs.x - sin(p.theta)*obs.y;
	global_obs.y_f = p.y + sin(p.theta)*obs.x + cos(p.theta)*obs.y;
	return global_obs;
}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	std::vector<LandmarkObs> new_observations;
	int matching_index = 0;
	double min_dist;
	for (int i = 0; i < predicted.size(); ++i){
		matching_index = 0;
		min_dist = numeric_limits<double>::max();
		for (int j = 0; j < observations.size(); ++j){
			double dist = sqrt((predicted[i].x - observations[j].x)*(predicted[i].x - observations[j].x)+(predicted[i].y - observations[j].y)*(predicted[i].y - observations[j].y));
			if (dist<min_dist){
				min_dist = dist;
				matching_index = j;
			}
		}

		LandmarkObs matching_obs = observations[matching_index];
		matching_obs.id = predicted[i].id;
		new_observations.push_back(matching_obs);
		observations.erase(observations.begin()+matching_index);

		if(!(observations.size()>0)){
			break;
		}		
	}

	observations = new_observations;
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {

	double sig_x = std_landmark[0];
	double sig_y = std_landmark[1];
	double gauss_norm= (1/(2 * M_PI * sig_x * sig_y));
	sensor_range = sensor_range*sensor_range;
	#pragma omp parallel for
	for (int i = 0; i < this->num_particles; ++i){
		double mu_x, mu_y, x_obs, y_obs, exponent, multiplier, weight;
		vector<int> associations;
		vector<double> sense_x;
		vector<double> sense_y;
		int matching_index;
		double obs_dist, sense_dist, min_dist;
		vector<Map::single_landmark_s> close_landmarks;
		for (int k = 0; k < map_landmarks.landmark_list.size(); ++k){
			sense_dist = ((map_landmarks.landmark_list[k].x_f - this->particles[i].x)*(map_landmarks.landmark_list[k].x_f - this->particles[i].x)+(map_landmarks.landmark_list[k].y_f - this->particles[i].y)*(map_landmarks.landmark_list[k].y_f - this->particles[i].y));
			if (sense_dist<sensor_range){
				close_landmarks.push_back(map_landmarks.landmark_list[k]);
			}
		}

		weight = 1;
		for (int j = 0; j < observations.size(); ++j){
			Map::single_landmark_s global_observation = local_to_global(this->particles[i],observations[j]);
			matching_index = 0;
			min_dist = numeric_limits<double>::max();
			for (int k = 0; k < close_landmarks.size(); ++k){
				obs_dist = ((close_landmarks[k].x_f - global_observation.x_f)*(close_landmarks[k].x_f - global_observation.x_f)+(close_landmarks[k].y_f - global_observation.y_f)*(close_landmarks[k].y_f - global_observation.y_f));
				if (obs_dist<min_dist){
					min_dist = obs_dist;
					matching_index = k;
				}
	
			}

			x_obs = global_observation.x_f;
			y_obs = global_observation.y_f;
			mu_x = close_landmarks[matching_index].x_f;
			mu_y = close_landmarks[matching_index].y_f;
			
			exponent= ((x_obs - mu_x)*(x_obs - mu_x))/(2 * sig_x*sig_x) + ((y_obs - mu_y)*(y_obs - mu_y))/(2 * sig_y*sig_y);
			multiplier = gauss_norm * exp(-exponent);

			if (multiplier>0){
				weight *= multiplier;
			}

			associations.push_back(close_landmarks[matching_index].id_i);
			sense_x.push_back(x_obs);
			sense_y.push_back(y_obs);



		}

		this->SetAssociations(this->particles[i], associations, sense_x, sense_y);
		this->particles[i].weight = weight;
		this->weights[i] = weight;
	}


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution


	default_random_engine gen;
	discrete_distribution<int> dist(weights.begin(), weights.end());
	std::vector<Particle> new_particles;
	for (int i = 0; i < this->num_particles; ++i){
		new_particles.push_back(this->particles[dist(gen)]);
	}

	this->particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}