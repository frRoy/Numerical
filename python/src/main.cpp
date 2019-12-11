#include <pybind11/pybind11.h>

int add(int i, int j) {
    return i + j;
}

// c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` example.cpp -o example`python3-config --extension-suffix`

PYBIND11_MODULE(pynumerical, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function which adds two numbers");
}