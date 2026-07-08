set -e

g++ -std=c++17 ./app/main.cpp \
  -I/usr/local/include \
  -I/opt/homebrew/include \
  -I/opt/homebrew/include/eigen3 \
  -Iinclude \
  -L/usr/local/lib \
  -L/opt/homebrew/lib \
  -lQuantLib \
  -Wl,-rpath,/usr/local/lib \
  -Wl,-rpath,/opt/homebrew/lib \

./a.out

rm ./a.out