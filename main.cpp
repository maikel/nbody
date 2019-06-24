// Copyright (c) 2019 Maikel Nadolski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "particle_data.hpp"
#include "run_simulation.hpp"

#include <Vc/vector.h>
#include <boost/program_options.hpp>
#include <fmt/format.h>
#include <iostream>

namespace nbody {

namespace po = boost::program_options;

void main(const po::variables_map& vm) {
  const std::size_t num_particles = vm["num_particles"].as<std::size_t>();
  const float time_step_size = vm["time_step_size"].as<float>();
  const std::size_t n_time_steps = vm["n_time_steps"].as<std::size_t>();
  particle_data<float> particles(num_particles);
  auto wtime = std::chrono::steady_clock::now();
  auto start = wtime;
  run_simulation(
      particles, time_step_size,
      [wtime, &start, n_time_steps](std::size_t n, float time_point) {
        auto stop = std::chrono::steady_clock::now();
        std::chrono::duration<double> diff = stop - start;
        std::chrono::duration<double> wdiff = stop - wtime;
        const int progress = int(100.f * float(n) / float(n_time_steps));
        fmt::print("[{:3}%] wall-time: {}s, time-step time: {:e}s\n", progress,
                   wdiff.count(), diff.count());
        start = stop;
        return n < n_time_steps;
      });
}

po::variables_map parse_command_line(int argc, char** argv) {
  po::options_description desc;
  desc.add_options()("help", "Print help message.")(
      "num_particles,p", po::value<std::size_t>()->default_value(100),
      "Number of particles in the simulation")(
      "n_time_steps", po::value<std::size_t>()->default_value(100),
      "Number of time steps taken")("time_step_size,dt",
                                    po::value<float>()->default_value(0.05f),
                                    "Time step size for each time step.");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc;
    exit(0);
  }
  return vm;
}

} // namespace nbody

int main(int argc, char** argv) {
  nbody::main(nbody::parse_command_line(argc, argv));
}
