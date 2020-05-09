#include "parallel.h"

#include <algorithm>
#include <functional>
#include <future>
#include <thread>

#include "ThreadPool.h"

namespace Prl2 {

Parallel::Parallel()
    : pool(ThreadPool(std::max(1U, std::thread::hardware_concurrency()))) {}

void Parallel::parallelFor1D(const std::function<void(unsigned int)>& job,
                             unsigned int nChunks, unsigned int n) {
  std::vector<std::future<void>> results;

  const unsigned int chunkSize = n / nChunks;
  for (unsigned int chunk_id = 0; chunk_id < nChunks; ++chunk_id) {
    results.push_back(pool.enqueue([chunk_id, chunkSize, job] {
      const unsigned int start_i = chunk_id * chunkSize;
      const unsigned int end_i = (chunk_id + 1) * chunkSize;
      for (unsigned int i = start_i; i < end_i; ++i) {
        job(i);
      }
    }));
  }

  for (auto&& result : results) {
    result.get();
  }
}

void Parallel::parallelFor2D(
    const std::function<void(unsigned int, unsigned int)>& job,
    unsigned int nChunks_x, unsigned int nChunks_y, unsigned int nx,
    unsigned int ny) {
  std::vector<std::future<void>> results;

  const unsigned int chunkSize_x = nx / nChunks_x;
  const unsigned int chunkSize_y = ny / nChunks_y;
  for (unsigned int chunk_y = 0; chunk_y < nChunks_y; ++chunk_y) {
    for (unsigned int chunk_x = 0; chunk_x < nChunks_x; ++chunk_x) {
      results.push_back(
          pool.enqueue([chunk_x, chunk_y, chunkSize_x, chunkSize_y, job] {
            const unsigned int start_x = chunk_x * chunkSize_x;
            const unsigned int end_x = (chunk_x + 1) * chunkSize_x;
            const unsigned int start_y = chunk_y * chunkSize_y;
            const unsigned int end_y = (chunk_y + 1) * chunkSize_y;
            for (unsigned int y = start_y; y < end_y; ++y) {
              for (unsigned int x = start_x; x < end_x; ++x) {
                job(x, y);
              }
            }
          }));
    }
  }

  for (auto&& result : results) {
    result.get();
  }
}

}  // namespace Prl2