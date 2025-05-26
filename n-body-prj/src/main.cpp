#include <iostream>

#include "GSimulation.hpp"
#include <functional>

int main(int argc, char **argv)
{
  int N;        // number of particles
  int sim_time; // simulation time
  real_type area_lim;

  if (argc == 4)
  {
    GSimulation sim;
    N = atoi(argv[1]);
    sim.set_number_of_particles(N);
    sim_time = atoi(argv[2]);
    sim.set_sim_time(sim_time);
    area_lim = atof(argv[3]);
    sim.set_area_size(area_lim);

    sim.start(std::function<void(void)>{});
  }
  else
  {
    std::cout << "Use: n_body_simulation.exe PARTICLES_NUM SIM_TIME AREA\n";
  }

  return 0;
}
