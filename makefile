

int_test_1: int_test_1.cpp
	g++ -Wall -c int_test_1.cpp $$(root-config --glibs --cflags --libs)
	g++ int_test_1.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o int_test_1.exe
