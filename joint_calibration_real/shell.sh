g++ -std=c++17 ./app/main.cpp \
    -I/usr/local/include \
    -I/opt/homebrew/include \
    -Iinclude \
    -L/usr/local/lib \
    -L/opt/homebrew/lib \
    -lQuantLib \
    -Wl,-rpath,/usr/local/lib \
    -o run_main

./run_main