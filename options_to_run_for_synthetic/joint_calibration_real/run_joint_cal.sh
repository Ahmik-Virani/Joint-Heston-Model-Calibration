g++ -std=c++17 -pthread ./app/main.cpp \
  -I../../../eigen-3.4.0 \
  -I../../../LBFGSpp/include \
  -I$HOME/local/include \
  -Iinclude \
  -L$HOME/local/lib \
  -lQuantLib \
  -Wl,-rpath,$HOME/local/lib \
  -o app/main
