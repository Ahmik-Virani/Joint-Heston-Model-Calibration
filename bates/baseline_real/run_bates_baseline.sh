g++ -std=c++17 ./app/main.cpp ./src/*.cpp \
  -I../../../eigen-3.4.0 \
  -I../../../LBFGSpp/include \
  -I$HOME/local/include \
  -Iinclude \
  -L$HOME/local/lib \
  -lQuantLib \
  -Wl,-rpath,$HOME/local/lib \
  -o app/main

./app/main
