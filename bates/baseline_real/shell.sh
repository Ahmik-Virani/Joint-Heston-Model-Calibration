set -e

g++ -std=c++17 ./app/main.cpp \
  ./src/params.cpp \
  ./src/surfacefit.cpp \
  ./src/particle_filters.cpp \
  ./src/surface.cpp \
  ./src/path.cpp \
  ./src/latent_path_mcmc.cpp \
  ./src/garch.cpp \
  -I/usr/local/include \
  -I/opt/homebrew/include \
  -I/opt/homebrew/include/eigen3 \
  -Iinclude \
  -L/usr/local/lib \
  -L/opt/homebrew/lib \
  -lQuantLib \
  -Wl,-rpath,/usr/local/lib \
  -Wl,-rpath,/opt/homebrew/lib \
  -o demo

./demo