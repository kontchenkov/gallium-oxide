#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <probable.hh>

using namespace probable;
const double density = 5.88 * units::g / pow(units::cm, 3);     // density
const double sound_velocity = 6.8e3 * units::m / units::s;      // sound velocity
const double acoustic_deformation_potential = 16.6 * units::eV; // acoustic deformation potential
const double nonpolar_optical_deformation_potential =
    8.5e8 * units::eV / units::cm; // acoustic deformation potential
const double eps_inf = 4.21;       // high-frequency dielectric constant
const double eps_static = 11.4;    // static dielectric constant;

struct AcousticScattering : public Scattering {
  AcousticScattering(const Material &m, double temperature) : Scattering(m, 0) {
    constant = pow(acoustic_deformation_potential, 2) * consts::kB * temperature * m.mass /
               (math::pi * pow(consts::hbar, 4) * density * pow(sound_velocity, 2));
  }
  double rate(const Vec3 &p) const { return constant * p.length(); }

private:
  double constant;
};

struct ImpurityScattering : public Scattering {
  ImpurityScattering(const Material &m, double temperature) : Scattering(m, 0) {
    const double z = 1;
    const double t = temperature / units::K;

    // аппроксимация радиуса экранирования
    const double a = 2.73e-11;
    const double b = 24.11;
    const double c = 5.69e-15;
    const double d = -1.99e-13;
    const double r02 = (a / (t - b) + c * t + d) * pow(units::cm, 2);

    // концентрация электронов
    const double n = 1e17 * (0.44 * atan(0.021 * t - 2.05) + 0.35) * pow(units::cm, -3);
    // концентрация акцепторов
    const double n_a = 4.2e16 * pow(units::cm, -3);
    // концентрация ионизированных примесей
    const double nci = n + 2 * n_a;

    constant = nci / math::pi * pow(z * r02 / consts::eps0 / eps_static, 2) *
               pow(consts::e / consts::hbar, 4) * m.mass;
    p02inv = r02 * pow(2 / consts::hbar, 2);
  }
  double rate(const Vec3 &p) const { return constant * p.length() / (1 + p.dot(p) * p02inv); }

private:
  double constant;
  double p02inv;
};

struct NonpolarOpticalAbsorptionScattering : public Scattering {
  NonpolarOpticalAbsorptionScattering(const Material &m, double temperature, double energy)
      : Scattering(m, -energy) {
    constant = pow(nonpolar_optical_deformation_potential, 2) * pow(m.mass, 1.5) /
               (sqrt(2) * math::pi * pow(consts::hbar, 2) * density * energy) /
               (exp(energy / consts::kB / temperature) - 1);
  }
  double rate(const Vec3 &p) const { return constant * sqrt(m.energy(p) - energy); }

private:
  double constant;
};

struct NonpolarOpticalEmissionScattering : public Scattering {
  NonpolarOpticalEmissionScattering(const Material &m, double temperature, double energy)
      : Scattering(m, energy) {
    constant = pow(nonpolar_optical_deformation_potential, 2) * pow(m.mass, 1.5) /
               (sqrt(2) * math::pi * pow(consts::hbar, 2) * density * energy) *
               (1 + 1 / (exp(energy / consts::kB / temperature) - 1));
  }
  double rate(const Vec3 &p) const {
    double e = m.energy(p);
    if (e < energy) {
      return 0;
    }
    return constant * sqrt(e - energy);
  }

private:
  double constant;
};

struct PolarOpticalAbsorptionScattering : public Scattering {
  PolarOpticalAbsorptionScattering(const Material &m, double temperature, double energy)
      : Scattering(m, -energy) {
    double n = 1 / (exp(energy / consts::kB / temperature) - 1);
    constant = pow(consts::e, 2) * energy / (2 * math::pi * consts::eps0 * pow(consts::hbar, 2)) *
               (1 / eps_inf - 1 / eps_static) * n;
  }
  PolarOpticalAbsorptionScattering(const Material &m, double temperature, double energy, double c)
      : Scattering(m, -energy) {
    double n = 1 / (exp(energy / consts::kB / temperature) - 1);
    constant = pow(c, 2) / (2 * math::pi * density * energy) * n;
  }
  double rate(const Vec3 &p) const {
    double v = m.velocity(p).length();
    if (v < 1e-20) {
      return constant * sqrt(m.mass / std::abs(2 * energy));
    }
    return constant / v * asinh(sqrt(m.energy(p) / std::abs(energy)));
  }

private:
  double constant;
};

struct PolarOpticalEmissionScattering : public Scattering {
  PolarOpticalEmissionScattering(const Material &m, double temperature, double energy)
      : Scattering(m, energy) {
    double n = 1 + 1 / (exp(energy / consts::kB / temperature) - 1);
    constant = pow(consts::e, 2) * energy / (2 * math::pi * consts::eps0 * pow(consts::hbar, 2)) *
               (1 / eps_inf - 1 / eps_static) * n;
  }
  PolarOpticalEmissionScattering(const Material &m, double temperature, double energy, double c)
      : Scattering(m, energy) {
    double n = 1 + 1 / (exp(energy / consts::kB / temperature) - 1);
    constant = pow(c, 2) / (2 * math::pi * density * energy) * n;
  }
  double rate(const Vec3 &p) const {
    double e = m.energy(p);
    if (e < energy) {
      return 0;
    }
    double v = m.velocity(p).length();
    return constant / v * acosh(sqrt(e / std::abs(energy)));
  }

private:
  double constant;
};

template <typename T> T parse(const std::string &s) {
  std::stringstream ss(s);
  T result;
  ss >> result;
  return result;
}

template <typename T> std::ostream &operator<<(std::ostream &s, std::vector<T> t) {
  s << "[";
  for (std::size_t i = 0; i < t.size(); i++) {
    s << t[i] << (i == t.size() - 1 ? "" : ", ");
  }
  return s << "]" << std::endl;
}

template <typename T> T sum(std::vector<T> t) {
  T s{};
  for (auto &v : t) {
    s += v;
  }
  return s;
}

int main(int argc, char const *argv[]) {
  if (argc != 4) {
    std::cout << "Invalid number of arguments\n";
    std::cout << "Usage: " << argv[0] << " <ansemble size> <temperature> <electric field>\n";
    return 1;
  }

  size_t ansemble_size = parse<int>(argv[1]);
  double temperature = parse<double>(argv[2]) * units::K;
  Vec3 electric_field{parse<double>(argv[3]) * units::V / units::m, 0, 0};
  Vec3 magnetic_field{0, 0, 0};

  Material gallium_oxide{0.29 * consts::me}; // electron effective mass

  std::vector<Scattering *> scattering_mechanisms{
      new AcousticScattering(gallium_oxide, temperature),
      new ImpurityScattering(gallium_oxide, temperature),
      new PolarOpticalAbsorptionScattering(
          gallium_oxide, temperature, 25e-3 * units::eV, 2.0e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalEmissionScattering(
          gallium_oxide, temperature, 25e-3 * units::eV, 2.0e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalAbsorptionScattering(
          gallium_oxide, temperature, 29e-3 * units::eV, 1.6e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalEmissionScattering(
          gallium_oxide, temperature, 29e-3 * units::eV, 1.6e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalAbsorptionScattering(
          gallium_oxide, temperature, 35e-3 * units::eV, 0.5e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalEmissionScattering(
          gallium_oxide, temperature, 35e-3 * units::eV, 0.5e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalAbsorptionScattering(
          gallium_oxide, temperature, 43e-3 * units::eV, 1.15e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalEmissionScattering(
          gallium_oxide, temperature, 43e-3 * units::eV, 1.15e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalAbsorptionScattering(
          gallium_oxide, temperature, 50e-3 * units::eV, 3.2e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalEmissionScattering(
          gallium_oxide, temperature, 50e-3 * units::eV, 3.2e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalAbsorptionScattering(
          gallium_oxide, temperature, 70e-3 * units::eV, 4.5e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalEmissionScattering(
          gallium_oxide, temperature, 70e-3 * units::eV, 4.5e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalAbsorptionScattering(
          gallium_oxide, temperature, 80e-3 * units::eV, 3.1e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalEmissionScattering(
          gallium_oxide, temperature, 80e-3 * units::eV, 3.1e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalAbsorptionScattering(
          gallium_oxide, temperature, 90e-3 * units::eV, 3.2e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalEmissionScattering(
          gallium_oxide, temperature, 90e-3 * units::eV, 3.2e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalAbsorptionScattering(
          gallium_oxide, temperature, 14e-3 * units::eV, 0.15e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalEmissionScattering(
          gallium_oxide, temperature, 14e-3 * units::eV, 0.15e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalAbsorptionScattering(
          gallium_oxide, temperature, 37e-3 * units::eV, 2.0e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalEmissionScattering(
          gallium_oxide, temperature, 37e-3 * units::eV, 2.0e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalAbsorptionScattering(
          gallium_oxide, temperature, 60e-3 * units::eV, 3.8e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalEmissionScattering(
          gallium_oxide, temperature, 60e-3 * units::eV, 3.8e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalAbsorptionScattering(
          gallium_oxide, temperature, 81e-3 * units::eV, 3.0e2 * units::eV / pow(units::nm, 2)),
      new PolarOpticalEmissionScattering(
          gallium_oxide, temperature, 81e-3 * units::eV, 3.0e2 * units::eV / pow(units::nm, 2))};

  double time_step = 1e-16 * units::s;
  double all_time = 1e-11 * units::s;

  // print info on stdout
  std::cout << "Ansemble size:    " << ansemble_size << "\n";
  std::cout << "Time step:        " << time_step / units::s << " s\n";
  std::cout << "Simulation time:  " << all_time / units::s << " s\n";
  std::cout << "Temperature:      " << temperature / units::K << " K\n";
  std::cout << "Electric field:   " << electric_field / units::V * units::m << " V/m\n";
  std::cout << "Magnetic field:   " << magnetic_field / units::T << " T\n";
  std::cout << "Scattering mechanisms:\n";
  for (size_t i = 0; i < scattering_mechanisms.size(); ++i) {
    std::cout << i + 1 << ": " << *scattering_mechanisms[i] << '\n';
  }
  auto results = simulate(gallium_oxide,
                          scattering_mechanisms,
                          temperature,
                          electric_field,
                          magnetic_field,
                          time_step,
                          all_time,
                          ansemble_size,
                          DumpFlags(DumpFlags::none));
  Vec3 average_velocity;
  Vec3 average_velocity2;
  std::vector<double> scattering_rates(scattering_mechanisms.size(), 0);
  for (std::size_t i = 0; i < results.size(); ++i) {
    average_velocity += (results[i].average_velocity - average_velocity) / (i + 1);
    average_velocity2 +=
        (results[i].average_velocity * results[i].average_velocity - average_velocity2) / (i + 1);
    for (std::size_t j = 0; j < scattering_mechanisms.size(); ++j) {
      scattering_rates[j] +=
          (results[i].scattering_count[j] / all_time * units::s - scattering_rates[j]) / (i + 1);
    }
  }
  Vec3 std_velocity = (average_velocity2 - average_velocity * average_velocity).sqrt();
  std::cout << "Average velocity: " << average_velocity << '\n';
  std::cout << "             std: " << std_velocity << '\n';
  std::cout << "Scattering rates: " << scattering_rates << '\n';
  std::cout << "           Total: " << sum(scattering_rates) << '\n';
  return 0;
}