// MIT License
//
// Copyright (c) 2024 TriLite https://github.com/MeshLite/TriLite
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

#ifndef TRILITE_GENERATOR_H
#define TRILITE_GENERATOR_H

#include <coroutine>
#include <iostream>
#include <optional>
#include <vector>

namespace TL {
// use Generator class from
// https://en.cppreference.com/w/cpp/coroutine/coroutine_handle while waiting
// for std::generator of c++23
template <std::movable T>
class Generator {
 public:
  struct promise_type {
    Generator<T> get_return_object() {
      return Generator{Handle::from_promise(*this)};
    }
    static std::suspend_always initial_suspend() noexcept { return {}; }
    static std::suspend_always final_suspend() noexcept { return {}; }
    std::suspend_always yield_value(T value) noexcept {
      current_value = std::move(value);
      return {};
    }
    // Disallow co_await in Generator coroutines.
    void await_transform() = delete;
    [[noreturn]]
    static void unhandled_exception() {
      throw;
    }

    std::optional<T> current_value;
  };

  using Handle = std::coroutine_handle<promise_type>;

  explicit Generator(const Handle coroutine) : m_coroutine{coroutine} {}

  Generator() = default;
  ~Generator() {
    if (m_coroutine) m_coroutine.destroy();
  }

  Generator(const Generator&) = delete;
  Generator& operator=(const Generator&) = delete;

  Generator(Generator&& other) noexcept : m_coroutine{other.m_coroutine} {
    other.m_coroutine = {};
  }
  Generator& operator=(Generator&& other) noexcept {
    if (this != &other) {
      if (m_coroutine) m_coroutine.destroy();
      m_coroutine = other.m_coroutine;
      other.m_coroutine = {};
    }
    return *this;
  }

  // Range-based for loop support.
  class Iter {
   public:
    void operator++() { m_coroutine.resume(); }
    const T& operator*() const { return *m_coroutine.promise().current_value; }
    bool operator==(std::default_sentinel_t) const {
      return !m_coroutine || m_coroutine.done();
    }

    explicit Iter(const Handle coroutine) : m_coroutine{coroutine} {}

   private:
    Handle m_coroutine;
  };

  Iter begin() {
    if (m_coroutine) m_coroutine.resume();
    return Iter{m_coroutine};
  }

  std::default_sentinel_t end() { return {}; }

 private:
  Handle m_coroutine;
};

template <std::movable T>
std::vector<T> ToVector(auto view) {
  std::vector<T> out_vector;
  for (T t : view) {
    out_vector.push_back(t);
  }
  return out_vector;
}
}  // namespace TL
#endif  // TRILITE_GENERATOR_H