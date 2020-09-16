#include <memory>

#include "lens-system/lens-system.h"

int main() {
  const auto film = std::make_shared<Film>(512, 512, 0.025, 0.025);
  const LensSystem lsys = LensSystem("../data/dgauss100mm.json", film);
  return 0;
}