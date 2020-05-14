#ifndef RANDOM_H
#define RANDOM_H
#include "samplers/rng.h"
#include "samplers/sampler.h"

namespace Prl2 {

// ただの乱数を返すSampler
class RandomSampler : public Sampler {
 public:
  RandomSampler(){};
  RandomSampler(uint64_t seed) : rng(RNG(seed)){};

  void setSeed(uint64_t seed) override { rng.setSeed(seed); };

  Real getNext() override { return rng.uniformReal(); }
  Vec2 getNext2D() override {
    return Vec2(rng.uniformReal(), rng.uniformReal());
  };

  std::unique_ptr<Sampler> clone(uint64_t seed) override {
    return std::unique_ptr<Sampler>(new RandomSampler(seed));
  };

 private:
  RNG rng;
};

}  // namespace Prl2

#endif
