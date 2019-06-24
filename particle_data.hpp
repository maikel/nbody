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

#ifndef NBODY_PARTICLE_DATA_HPP
#define NBODY_PARTICLE_DATA_HPP

#include <Vc/vector.h>
#include <type_traits>
#include <vector>

namespace nbody {

template <typename T, typename Allocator = std::allocator<T>>
class particle_data {
public:
  particle_data(std::size_t size, const Allocator& alloc = Allocator());

  template <typename Generator>
  particle_data(std::size_t size, Generator gen,
                const Allocator& alloc = Allocator());

  T* x() noexcept;
  const T* x() const noexcept;

  T* y() noexcept;
  const T* y() const noexcept;

  T* velocity_x() noexcept;
  const T* velocity_x() const noexcept;

  T* velocity_y() noexcept;
  const T* velocity_y() const noexcept;

  T* mass() noexcept;
  const T* mass() const noexcept;

  T* force_x() noexcept;
  const T* force_x() const noexcept;

  T* force_y() noexcept;
  const T* force_y() const noexcept;

  std::size_t size() const noexcept;
  std::size_t padded_size() const noexcept;

private:
  T* at_component(int n) noexcept;
  const T* at_component(int n) const noexcept;

  std::vector<T, Allocator> buffer_;
  std::size_t n_particles_;
  std::size_t offset_;
};

template <typename T>
std::size_t compute_particle_buffer_offset(std::size_t n_particles) {
  const std::size_t pack_size = Vc::native_simd<T>::size();
  const std::size_t missing =
      (pack_size - (n_particles % pack_size)) % pack_size;
  return n_particles + missing;
}

template <typename T, typename Allocator>
particle_data<T, Allocator>::particle_data(std::size_t size,
                                           const Allocator& alloc)
    : buffer_(alloc), n_particles_(size),
      offset_(compute_particle_buffer_offset<T>(n_particles_)) {
  buffer_.resize(7 * offset_);
}

template <typename T, typename Allocator>
T* particle_data<T, Allocator>::at_component(int n) noexcept {
  return buffer_.data() + n * offset_;
}

template <typename T, typename Allocator>
const T* particle_data<T, Allocator>::at_component(int n) const noexcept {
  return buffer_.data() + n * offset_;
}

template <typename T, typename Allocator>
T* particle_data<T, Allocator>::x() noexcept {
  return at_component(0);
}
template <typename T, typename Allocator>
const T* particle_data<T, Allocator>::x() const noexcept {
  return at_component(0);
}

template <typename T, typename Allocator>
T* particle_data<T, Allocator>::y() noexcept {
  return at_component(1);
}
template <typename T, typename Allocator>
const T* particle_data<T, Allocator>::y() const noexcept {
  return at_component(1);
}

template <typename T, typename Allocator>
T* particle_data<T, Allocator>::velocity_x() noexcept {
  return at_component(2);
}
template <typename T, typename Allocator>
const T* particle_data<T, Allocator>::velocity_x() const noexcept {
  return at_component(2);
}

template <typename T, typename Allocator>
T* particle_data<T, Allocator>::velocity_y() noexcept {
  return at_component(3);
}
template <typename T, typename Allocator>
const T* particle_data<T, Allocator>::velocity_y() const noexcept {
  return at_component(3);
}

template <typename T, typename Allocator>
T* particle_data<T, Allocator>::mass() noexcept {
  return at_component(4);
}
template <typename T, typename Allocator>
const T* particle_data<T, Allocator>::mass() const noexcept {
  return at_component(4);
}

template <typename T, typename Allocator>
T* particle_data<T, Allocator>::force_x() noexcept {
  return at_component(5);
}
template <typename T, typename Allocator>
const T* particle_data<T, Allocator>::force_x() const noexcept {
  return at_component(5);
}

template <typename T, typename Allocator>
T* particle_data<T, Allocator>::force_y() noexcept {
  return at_component(6);
}
template <typename T, typename Allocator>
const T* particle_data<T, Allocator>::force_y() const noexcept {
  return at_component(6);
}

template <typename T, typename Allocator>
std::size_t particle_data<T, Allocator>::size() const noexcept {
  return n_particles_;
}

template <typename T, typename Allocator>
std::size_t particle_data<T, Allocator>::padded_size() const noexcept {
  return offset_;
}

} // namespace nbody

#endif // !PARTICLE_DATA_HPP
