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
	std::vector<Particle> new_particles;

	for (int i = 0; i < this->num_particles; ++i){
		Particle one_particle;
		one_particle = this->particles[i];
	    if (fabs(yaw_rate) > 0.000000001) {
	        one_particle.x = this->particles[i].x + velocity/yaw_rate * ( sin (this->particles[i].theta + yaw_rate*delta_t) - sin(this->particles[i].theta));
	        one_particle.y = this->particles[i].y + velocity/yaw_rate * ( cos(this->particles[i].theta) - cos(this->particles[i].theta+yaw_rate*delta_t) );
	        one_particle.theta = this->particles[i].theta + yaw_rate*delta_t;
	    }
	    else {
	        one_particle.x = this->particles[i].x + velocity*delta_t*cos(this->particles[i].theta);
	        one_particle.y = this->particles[i].y + velocity*delta_t*sin(this->particles[i].theta);
	        one_particle.theta = this->particles[i].theta;
	    }
		normal_distribution<double> dist_x(one_particle.x, std_pos[0]);
		normal_distribution<double> dist_y(one_particle.y, std_pos[1]);
		normal_distribution<double> dist_theta(one_particle.theta, std_pos[3]);
	    one_particle.x = dist_x(gen);
	    one_particle.y = dist_y(gen);
	    one_particle.theta = dist_theta(gen);
		new_particles.push_back(one_particle);
	}

	this->particles = new_particles;

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

Map::single_landmark_s local_to_global(Particle p, LandmarkObs obs){
	Map::single_landmark_s global_obs;
	global_obs.id_i = obs.id;
	//Rotate
	global_obs.x_f = cos(p.theta)*obs.x - sin(p.theta)*obs.y;
	global_obs.y_f = sin(p.theta)*obs.x + cos(p.theta)*obs.y;
	//Translate
	global_obs.x_f += p.x;
	global_obs.y_f += p.y;
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
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	std::vector<LandmarkObs> associated_observations;
	double mu_x, mu_y, sig_x, sig_y, x_obs, y_obs, gauss_norm, exponent, multiplier, weight;
	sig_x = std_landmark[0];
	sig_y = std_landmark[1];
	vector<int> associations;
	vector<double> sense_x;
	vector<double> sense_y;

	for (int i = 0; i < this->num_particles; ++i){

		associations.clear();
		sense_x.clear();
		sense_y.clear();

		std::vector<LandmarkObs> predicted_landmark;
		std::vector<Map::single_landmark_s> close_landmarks;
		for (int k = 0; k < map_landmarks.landmark_list.size(); ++k){
			double dist = sqrt((map_landmarks.landmark_list[k].x_f - this->particles[i].x)*(map_landmarks.landmark_list[k].x_f - this->particles[i].x)+(map_landmarks.landmark_list[k].y_f - this->particles[i].y)*(map_landmarks.landmark_list[k].y_f - this->particles[i].y));
			if (dist<=sensor_range){
				LandmarkObs one_landmark = global_to_local(this->particles[i],map_landmarks.landmark_list[k]);
				predicted_landmark.push_back(one_landmark);
				close_landmarks.push_back(map_landmarks.landmark_list[k]);
			}
		}

		if (!(predicted_landmark.size()>0) or !(observations.size()>0)){
			this->SetAssociations(this->particles[i], associations, sense_x, sense_y);
			continue;
		}

		associated_observations = observations;
		dataAssociation(predicted_landmark, associated_observations);

		weight = 1;
		for (int j = 0; j < associated_observations.size(); ++j){
			for (int k = 0; k < close_landmarks.size(); ++k){
				if(associated_observations[j].id == close_landmarks[k].id_i){

					Map::single_landmark_s global_observation = local_to_global(this->particles[i],associated_observations[j]);

					x_obs = global_observation.x_f;
					y_obs = global_observation.y_f;
					mu_x = close_landmarks[k].x_f;
					mu_y = close_landmarks[k].y_f;

					gauss_norm= (1/(2 * M_PI * sig_x * sig_y));
					exponent= ((x_obs - mu_x)*(x_obs - mu_x))/(2 * sig_x*sig_x) + ((y_obs - mu_y)*(y_obs - mu_y))/(2 * sig_y*sig_y);
					multiplier = gauss_norm * exp(-exponent);

					if (multiplier>0){
						weight *= multiplier;
					}

					associations.push_back(global_observation.id_i);
					sense_x.push_back(x_obs);
					sense_y.push_back(y_obs);

				}
			}
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
