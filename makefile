

int_test_1: int_test_1.cpp
	g++ -Wall -c int_test_1.cpp
	g++ int_test_1.o -lgsl -lgslcblas -lm -o int_test_1.exe
