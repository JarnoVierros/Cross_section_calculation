all: inclusive_ccbar_cross_section.cpp intQ2_inclusive_ccbar_cross_section.cpp
	make inclusive_ccbar_cross_section
	make intQ2_inclusive_ccbar_cross_section
	make intx_inclusive_ccbar_cross_section

inclusive_ccbar_cross_section: inclusive_ccbar_cross_section.cpp
	g++ -Wall -c inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o inclusive_ccbar_cross_section.exe
	
intQ2_inclusive_ccbar_cross_section: intQ2_inclusive_ccbar_cross_section.cpp
	g++ -Wall -c intQ2_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ intQ2_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o intQ2_inclusive_ccbar_cross_section.exe

intx_inclusive_ccbar_cross_section: intx_inclusive_ccbar_cross_section.cpp
	g++ -Wall -c intx_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ intx_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o intx_inclusive_ccbar_cross_section.exe

clean:
	-rm *.exe
	-rm *.o
