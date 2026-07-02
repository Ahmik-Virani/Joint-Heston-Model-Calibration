g++ -std=c++17 ./app/main.cpp ./src/*.cpp \
  -I../../eigen-3.4.0 \
  -I../../LBFGSpp/include \
  -I$HOME/local/include \
  -Iinclude \
  -L$HOME/local/lib \
  -lQuantLib \
  -Wl,-rpath,$HOME/local/lib \
  -o app/main

./app/main

# g++ -std=c++17 ./app/garch_calibration.cpp ./src/*.cpp \
#   -I../../eigen-3.4.0 \
#   -I../../LBFGSpp/include \
#   -I$HOME/local/include \
#   -Iinclude \
#   -L$HOME/local/lib \
#   -lQuantLib \
#   -Wl,-rpath,$HOME/local/lib \
#   -o app/garch

# ./app/garch

# g++ -std=c++17 ./app/check_ppath.cpp ./src/garch.cpp ./src/path.cpp ./src/surface.cpp \
#   -I../../eigen-3.4.0 \
#   -I../../LBFGSpp/include \
#   -I$HOME/local/include \
#   -Iinclude \
#   -L$HOME/local/lib \
#   -lQuantLib \
#   -Wl,-rpath,$HOME/local/lib \
#   -o app/debug_ppath
