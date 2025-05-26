#include "GSimulation.hpp"
#include "cpu_time.hpp"
#include <functional>

GSimulation::GSimulation()
{
  std::cout << "===============================" << std::endl;
  std::cout << " Initialize Gravity Simulation" << std::endl;
  set_npart(10000);
  set_simtime(1);
  _impulse_ndim[0] = _impulse_ndim[1] = _impulse_ndim[2] = 0;
  _drop = false;
  _tstep = 0.1;
  _area_lim = 100;
  particles = new Particle[get_npart()];
}

GSimulation::~GSimulation() {
  delete[] particles;
}

void GSimulation::set_number_of_particles(int N)
{
  set_npart(N);
}

void GSimulation::set_sim_time(real_type N)
{
  _simtime = N;
}

void GSimulation::set_area_size(real_type N)
{
  _area_lim = N;
}

void GSimulation::init_pos()
{
  std::random_device rd;
  std::mt19937 gen(42);
  real_type l = 0, r = _area_lim;
  std::uniform_real_distribution<real_type> unif_d(l, r);

  for (int i = 0; i < get_npart(); ++i)
  {
    particles[i].pos[0] = unif_d(gen);
    particles[i].pos[1] = unif_d(gen);
    particles[i].pos[2] = unif_d(gen);
  }
}

void GSimulation::init_vel()
{
  std::random_device rd;
  std::mt19937 gen(42);
  std::uniform_real_distribution<real_type> unif_d(-1.0, 1.0);

  for (int i = 0; i < get_npart(); ++i)
  {
    particles[i].vel[0] = unif_d(gen) * 1e-9;
    particles[i].vel[1] = unif_d(gen) * 1e-9;
    particles[i].vel[2] = unif_d(gen) * 1e-9;
  }
}

void GSimulation::init_mass()
{
  real_type n = static_cast<real_type>(get_npart());
  std::random_device rd;
  std::mt19937 gen(42);
  std::uniform_real_distribution<real_type> unif_d(0.0, 1.0);
  real_type mass_sum = 0;

  for (int i = 0; i < get_npart(); ++i)
  {
    particles[i].mass = n * unif_d(gen);
    mass_sum += particles[i].mass;
  }
}

// prevents explosion in the case the particles are really close to each other
static constexpr real_type softeningSquared = 1e-6;
static constexpr real_type G = 6.67430e-11;

template <typename type, typename added_type>
void sum_with_correction(type &sum, const added_type &value_to_add, type &correction)
{
  type corrected = value_to_add - correction;
  type new_sum = sum + corrected;
  correction = (new_sum - sum) - corrected;
  sum = new_sum;
}

real_type GSimulation::compute_impulse()
{
  real_type correction[] = {0, 0, 0};
  _impulse_ndim[0] = _impulse_ndim[1] = _impulse_ndim[2] = 0;
  for (int i = 0; i < _npart; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      real_type curr_impulse = particles[i].mass * particles[i].vel[j];
      sum_with_correction(_impulse_ndim[j], curr_impulse, correction[j]);
    }
  }
  return sqrt(pow(_impulse_ndim[0], 2) + pow(_impulse_ndim[1], 2) + pow(_impulse_ndim[2], 2));
}

real_type GSimulation::compute_k_energy()
{
  double energy = 0.;
  double correction = 0.;
  for (int i = 0; i < _npart; ++i)
  {
    double curr_energy = particles[i].mass * (particles[i].vel[0] * particles[i].vel[0] +
                                              particles[i].vel[1] * particles[i].vel[1] +
                                              particles[i].vel[2] * particles[i].vel[2]);
    sum_with_correction(energy, curr_energy, correction);
  }
  return energy / 2;
}

real_type GSimulation::compute_p_energy()
{
  double p_energy = 0.;
  double correction = 0.;
  for (int i = 0; i < _npart; ++i)
  {
    for (int j = 0; j < _npart; ++j)
    {
      if (i == j)
        continue;
      double dx = particles[i].pos[0] - particles[j].pos[0];
      double dy = particles[i].pos[1] - particles[j].pos[1];
      double dz = particles[i].pos[2] - particles[j].pos[2];

      double distanceSqr = dx * dx + dy * dy + dz * dz + softeningSquared;
      double distanceInv = 1.0 / sqrt(distanceSqr);
      double curr_energy = -G * particles[j].mass * distanceInv * particles[i].mass;
      sum_with_correction(p_energy, curr_energy, correction);
    }
  }
  return p_energy / 2;
}

void GSimulation::compute_f(Particle *y, Particle *f)
{
  for (int i = 0; i < _npart; i++)
  {
    f[i].vel[0] = 0.;
    f[i].vel[1] = 0.;
    f[i].vel[2] = 0.;

    real_type corr[3] = {0, 0, 0};

    for (int j = 0; j < _npart; j++)
    {
      if (i != j)
      {
        real_type dx, dy, dz;
        real_type distanceSqr = 0.0;
        real_type distanceInv = 0.0;

        dx = y[i].pos[0] - y[j].pos[0];
        dy = y[i].pos[1] - y[j].pos[1];
        dz = y[i].pos[2] - y[j].pos[2];

        distanceSqr = dx * dx + dy * dy + dz * dz + softeningSquared;
        distanceInv = 1.0 / sqrt(distanceSqr);

        sum_with_correction(f[i].vel[0], -dx * G * y[j].mass * distanceInv * distanceInv * distanceInv, corr[0]);
        sum_with_correction(f[i].vel[1], -dy * G * y[j].mass * distanceInv * distanceInv * distanceInv, corr[1]);
        sum_with_correction(f[i].vel[2], -dz * G * y[j].mass * distanceInv * distanceInv * distanceInv, corr[2]);    
      }
    }
    f[i].pos[0] = y[i].vel[0];
    f[i].pos[1] = y[i].vel[1];
    f[i].pos[2] = y[i].vel[2];
  }
}

void GSimulation::update_system(Particle *dst, const Particle *src, const Particle *f, double dt)
{
  for (int i = 0; i < _npart; ++i)
  {
    dst[i].pos[0] = src[i].pos[0] + f[i].pos[0] * dt;
    dst[i].pos[1] = src[i].pos[1] + f[i].pos[1] * dt;
    dst[i].pos[2] = src[i].pos[2] + f[i].pos[2] * dt;

    dst[i].vel[0] = src[i].vel[0] + f[i].vel[0] * dt;
    dst[i].vel[1] = src[i].vel[1] + f[i].vel[1] * dt;
    dst[i].vel[2] = src[i].vel[2] + f[i].vel[2] * dt;
  }
}

void GSimulation::start(std::function<void(void)> cb)
{
  init_nsteps();
  init_pos();
  init_vel();
  init_mass();

  real_type energy_k, energy_p;
  real_type dt = get_tstep();

  Particle* f1 = new Particle[get_npart()];
  // Particle* f2 = new Particle[get_npart()];
  // Particle* f3 = new Particle[get_npart()];
  // Particle* f4 = new Particle[get_npart()];

  // Particle* tmp = new Particle[get_npart()];
  // for (int i = 0; i < get_npart(); ++i) {
  //   tmp[i].mass = particles[i].mass;
  // }

  energy_k = compute_k_energy();
  energy_p = compute_p_energy();
  _energy = energy_k + energy_p;
  _impulse = compute_impulse();

  std::cout << "Initial system energy k: " << energy_k << " p:" << energy_p << " Sum: " << _energy << " Impulse: " << _impulse << std::endl;

  print_header();

  double _totTime = 0.;

  CPUTime time;
  double ts0 = 0;
  double ts1 = 0;

  const double t0 = time.start();
  for (int s = 1; s <= get_nsteps(); ++s)
  { // start of the time step loop
    ts0 += time.start();

    // Simple Euler Method
    compute_f(particles, f1);
    update_system(particles, particles, f1, dt);

    // только эта часть кода (выше) должна измениться при добавлении вашего варианта.
    // при добавлении своего варианта нужно будет лишь поменять кол-во вызовов функций выше
    // допускается введение доп. вспомогательных массивов (аналогично tmp_part), то есть
    // их может быть несколько для разных вариантов.

    if (cb) {
      cb();
      if (_drop) break;
    }

    energy_k = compute_k_energy();
    energy_p = compute_p_energy();

    real_type curr_energy = energy_k + energy_p;

    double curr_impulse = compute_impulse();

    ts1 += time.stop();
    {
      std::cout << " "
                << std::left << std::setw(8) << s
                << std::left << std::setprecision(5) << std::setw(8) << s * get_tstep()
                << std::left << std::setprecision(9) << std::setw(16) << fabs(100 * (curr_energy - _energy) / _energy)
                << std::left << std::setprecision(9) << std::setw(16) << fabs(100 * (curr_impulse - _impulse) / _impulse)
                << std::left << std::setprecision(5) << std::setw(16) << (ts1 - ts0)
                << std::endl;
      ts0 = 0;
      ts1 = 0;
    }
  } // end of the time step loop

  const double t1 = time.stop();
  _totTime = (t1 - t0);

  std::cout << std::endl;
  std::cout << "# Total Time (s)     : " << _totTime << std::endl;
  std::cout << "===============================" << std::endl;
}

void GSimulation ::print_header()
{

  std::cout << " nPart = " << get_npart() << "; "
            << "nSteps = " << get_nsteps() << "; "
            << "dt = " << get_tstep() << std::endl;

  std::cout << "------------------------------------------------" << std::endl;
  std::cout << " "
            << std::left << std::setw(8) << "s"
            << std::left << std::setw(8) << "dt"
            << std::left << std::setw(16) << "d_energy"
            << std::left << std::setw(16) << "d_impulse"
            << std::left << std::setw(16) << "time (s)"
            << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
}
