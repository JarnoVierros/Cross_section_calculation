all: inclusive_ccbar_cross_section.cpp intQ2_inclusive_ccbar_cross_section.cpp
	make inclusive_ccbar_cross_section
	make exclusive_ccbar_cross_section
	make QAG_exclusive_ccbar_cross_section
	make intQ2_inclusive_ccbar_cross_section
	make intx_inclusive_ccbar_cross_section
	make W_inclusive_ccbar_cross_section
	make reduced_inclusive_ccbar_cross_section
	make parameter_fitter_inclusive_ccbar
	make fit_data_comparison_inclusive_ccbar
	make analytical_QAG_exclusive_ccbar_cross_section

inclusive_ccbar_cross_section: inclusive_ccbar_cross_section.cpp
	g++ -Wall -c inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o inclusive_ccbar_cross_section.exe
	
exclusive_ccbar_cross_section: exclusive_ccbar_cross_section.cpp
	g++ -Wall -c exclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ exclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o exclusive_ccbar_cross_section.exe

QAG_exclusive_ccbar_cross_section: QAG_exclusive_ccbar_cross_section.cpp
	g++ -Wall -c QAG_exclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ QAG_exclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o QAG_exclusive_ccbar_cross_section.exe

analytical_QAG_exclusive_ccbar_cross_section: analytical_QAG_exclusive_ccbar_cross_section.cpp
	g++ -Wall -c analytical_QAG_exclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ analytical_QAG_exclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o analytical_QAG_exclusive_ccbar_cross_section.exe
	
intQ2_inclusive_ccbar_cross_section: intQ2_inclusive_ccbar_cross_section.cpp
	g++ -Wall -c intQ2_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ intQ2_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o intQ2_inclusive_ccbar_cross_section.exe

intx_inclusive_ccbar_cross_section: intx_inclusive_ccbar_cross_section.cpp
	g++ -Wall -c intx_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ intx_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o intx_inclusive_ccbar_cross_section.exe

W_inclusive_ccbar_cross_section: W_inclusive_ccbar_cross_section.cpp
	g++ -Wall -c W_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ W_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o W_inclusive_ccbar_cross_section.exe

reduced_inclusive_ccbar_cross_section: reduced_inclusive_ccbar_cross_section.cpp
	g++ -Wall -c reduced_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ reduced_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o reduced_inclusive_ccbar_cross_section.exe

parameter_fitter_inclusive_ccbar: parameter_fitter_inclusive_ccbar.cpp
	g++ -Wall -c parameter_fitter_inclusive_ccbar.cpp $$(root-config --glibs --cflags --libs)
	g++ parameter_fitter_inclusive_ccbar.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o parameter_fitter_inclusive_ccbar.exe

fit_data_comparison_inclusive_ccbar: fit_data_comparison_inclusive_ccbar.cpp
	g++ -Wall -c fit_data_comparison_inclusive_ccbar.cpp $$(root-config --glibs --cflags --libs)
	g++ fit_data_comparison_inclusive_ccbar.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o fit_data_comparison_inclusive_ccbar.exe

clean:
	-rm *.exe
	-rm *.o
