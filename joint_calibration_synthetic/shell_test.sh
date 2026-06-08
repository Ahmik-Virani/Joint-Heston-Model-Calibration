g++ -std=c++17 test.cpp \
    -I/usr/local/include \
    -I/opt/homebrew/include \
    -L/usr/local/lib \
    -lQuantLib \
    -Wl,-rpath,/usr/local/lib \
    -o test

./test 