# g++ -std=c++17 test.cpp \
#     -I/usr/local/include \
#     -I/opt/homebrew/include \
#     -I../../eigen-3.4.0 \
#     -L/usr/local/lib \
#     -lQuantLib \
#     -Wl,-rpath,/usr/local/lib \
#     -o test

g++ -std=c++17 test.cpp \
  -I$HOME/local/include \
  -I../../eigen-3.4.0 \
  -Iinclude \
  -L$HOME/local/lib \
  -lQuantLib \
  -Wl,-rpath,$HOME/local/lib \
  -o test



./test 
