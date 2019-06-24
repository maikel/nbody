// Copyright (c) 2018 Maikel Nadolski
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

#include "run_simulation.hpp"

namespace nbody {
namespace {
void drift(float* x, float* y, const float* vx, const float* vy, float dt,
           std::size_t n) noexcept {
  using Vc::float_v;
  const std::size_t pack_size = float_v::size();
  assert(n % pack_size == 0);
  float_v dt_(dt);
  while (n >= pack_size) {
    float_v x_(x, Vc::Aligned);
    float_v y_(y, Vc::Aligned);
    const float_v vy_(vy, Vc::Aligned);
    const float_v vx_(vx, Vc::Aligned);
    x_ += dt_ * vx_;
    y_ += dt_ * vy_;
    x_.store(x, Vc::Aligned);
    y_.store(y, Vc::Aligned);
    n -= pack_size;
    x += pack_size;
    y += pack_size;
    vx += pack_size;
    vy += pack_size;
  }
}

void drift(particle_data<float>& particles, float dt) noexcept {
  drift(particles.x(), particles.y(), particles.velocity_x(),
        particles.velocity_y(), dt, particles.padded_size());
}

template <std::size_t N> std::size_t align_for(std::size_t n) noexcept {
  return n + (N - (n % N)) % N;
}

void accumulate_forces(float* fx, float* fy, const float* x, const float* y,
                       const float* mass, std::size_t n) noexcept {
  using Vc::float_v;
  const std::size_t pack_size = float_v::size();
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t j0 = align_for<pack_size>(i);
    for (std::size_t j = i + 1; j < j0; ++j) {
      const float dx = x[i] - x[j];
      const float dy = y[i] - y[j];
      const float dist2 = dx * dx + dy * dy;
      const float dist = std::sqrt(dist2);
      const float fxij = mass[i] * mass[j] / dist2 * (dx / dist);
      const float fyij = mass[i] * mass[j] / dist2 * (dy / dist);
      fx[i] += fxij;
      fx[j] += fxij;
      fy[i] += fyij;
      fy[j] += fyij;
    }
    const float_v xi(x[i]);
    const float_v yi(y[i]);
    const float_v mi(mass[i]);
    for (std::size_t j = j0; j < n; j += pack_size) {
      const float_v xj(x + j, Vc::Aligned);
      const float_v yj(x + j, Vc::Aligned);
      const float_v dx = xi - xj;
      const float_v dy = yi - yj;
      const float_v dist2 = dx * dx + dy * dy;
      const float_v rdist = rsqrt(dist2);
      const float_v mj(mass + j, Vc::Aligned);
      const float_v fxij = mi * mj / dist2 * (dx * rdist);
      const float_v fyij = mi * mj / dist2 * (dy * rdist);
      fx[i] += fxij.sum();
      fy[i] += fyij.sum();
      fxij.store(fx + j, Vc::Aligned);
      fyij.store(fy + j, Vc::Aligned);
    }
  }
}

void kick(particle_data<float>& particles, float dt) noexcept {
  accumulate_forces(particles.force_x(), particles.force_y(), particles.x(),
                    particles.y(), particles.mass(), particles.size());
  drift(particles.velocity_x(), particles.velocity_y(), particles.force_x(),
        particles.force_y(), dt, particles.padded_size());
}
} // namespace

void integrate_in_time(particle_data<float>& particles, float dt) noexcept {
  kick(particles, dt);
  drift(particles, dt);
}

} // namespace nbody
