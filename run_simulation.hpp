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

#ifndef NBODY_RUN_SIMULATION_HPP
#define NBODY_RUN_SIMULATION_HPP

#include "particle_data.hpp"

namespace nbody {

void integrate_in_time(particle_data<float>& particles, float dt) noexcept;

template <typename Callback>
void run_simulation(particle_data<float>& particles, float dt,
                    Callback callback) {
  float time_point{0.f};
  std::size_t counter{0};
  do {
    integrate_in_time(particles, dt);
    counter += 1;
    time_point += dt;
  } while (callback(counter, time_point));
}

} // namespace nbody

#endif // !RUN_SIMULATION_HPP
