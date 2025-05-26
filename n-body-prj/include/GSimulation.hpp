#ifndef _GSIMULATION_HPP
#define _GSIMULATION_HPP

#include <random>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <functional>

#include <omp.h>

#include "Particle.hpp"

class GSimulation
{
public:
  GSimulation();
  ~GSimulation();

  void set_number_of_particles(int N);
  void set_sim_time(real_type N);
  void set_area_size(real_type N);
  void start(std::function<void(void)> cb);

  Particle* get_parts() { return particles; }
  size_t get_part_num() { return _npart; }
  real_type get_area_lim() { return _area_lim; }
  void have_to_drop() { _drop = true; }

private:
  Particle *particles;

  bool _drop;

  real_type _area_lim;
  int _npart;         // number of particles
  int _nsteps;        // number of integration steps
  real_type _tstep;   // time step of the simulation
  real_type _simtime; // total simulation time

  real_type _energy;  // energy of the system
  real_type _impulse; // impulse of the system
  real_type _impulse_ndim[3];

  void init_pos();
  void init_vel();
  void init_mass();

  inline void set_npart(const int &N) { _npart = N; }
  inline int get_npart() const { return _npart; }

  inline void set_simtime(const real_type &time) { _simtime = time; }
  inline real_type get_simtime() const { return _simtime; }

  inline int get_nsteps() const { return _nsteps; }

  inline void init_nsteps() { _nsteps = _simtime / _tstep; };
  inline real_type get_tstep() { return _tstep; };

  real_type compute_impulse();
  real_type compute_k_energy();
  real_type compute_p_energy();
  void compute_f(Particle *y, Particle *f);
  void update_system(Particle *dst, const Particle *src_1, const Particle *f, double dt);

  void print_header();
};

#endif
