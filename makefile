all: inclusive_ccbar_cross_section.cpp intQ2_inclusive_ccbar_cross_section.cpp
	make inclusive_ccbar_cross_section
	make old_inclusive_ccbar_cross_section
	make Jb_inclusive_ccbar_cross_section
	make J_data_comparison_inclusive_ccbar
	make complicated_exclusive_ccbar_cross_section
	make intQ2_inclusive_ccbar_cross_section
	make intx_inclusive_ccbar_cross_section
	make W_inclusive_ccbar_cross_section
	make reduced_inclusive_ccbar_cross_section
	make parameter_fitter_inclusive_ccbar
	make old_parameter_fitter_inclusive_ccbar
	make fit_data_comparison_inclusive_ccbar
	make J_data_comparison_inclusive_ccbar
	make simplified_exclusive_ccbar_cross_section
	make alternative_exclusive_ccbar_cross_section
	make diffractive_inclusive_ratio
	make dipole_amp_reader
	make plot_comparator
	make J_diff_ccbar_cross_section
	make J_differential_diffractive_sigma
	make J_numerical_differential_diffractive_sigma
	make J_simplified_differential_diffractive_sigma
	make multiplotter
	make J_LHC_inclusive_ccbar_cross_section
	make 	make J_LHC_inclusive_ccbar_cross_section

Jb_inclusive_ccbar_cross_section: Jb_inclusive_ccbar_cross_section.cpp
	g++ -Wall -c Jb_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ Jb_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o Jb_inclusive_ccbar_cross_section.exe
	
J_inclusive_ccbar_cross_section: J_inclusive_ccbar_cross_section.cpp
	g++ -Wall -c J_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ J_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o J_inclusive_ccbar_cross_section.exe

inclusive_ccbar_cross_section: inclusive_ccbar_cross_section.cpp
	g++ -Wall -c inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o inclusive_ccbar_cross_section.exe

old_inclusive_ccbar_cross_section: old_inclusive_ccbar_cross_section.cpp
	g++ -Wall -c old_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ old_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o old_inclusive_ccbar_cross_section.exe

complicated_exclusive_ccbar_cross_section: complicated_exclusive_ccbar_cross_section.cpp
	g++ -Wall -c complicated_exclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ complicated_exclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o complicated_exclusive_ccbar_cross_section.exe

simplified_exclusive_ccbar_cross_section: simplified_exclusive_ccbar_cross_section.cpp
	g++ -Wall -c simplified_exclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ simplified_exclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o simplified_exclusive_ccbar_cross_section.exe
	
alternative_exclusive_ccbar_cross_section: alternative_exclusive_ccbar_cross_section.cpp
	g++ -Wall -c alternative_exclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ alternative_exclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o alternative_exclusive_ccbar_cross_section.exe
	
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

old_parameter_fitter_inclusive_ccbar: old_parameter_fitter_inclusive_ccbar.cpp
	g++ -Wall -c old_parameter_fitter_inclusive_ccbar.cpp $$(root-config --glibs --cflags --libs)
	g++ old_parameter_fitter_inclusive_ccbar.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o old_parameter_fitter_inclusive_ccbar.exe

J_data_comparison_inclusive_ccbar: J_data_comparison_inclusive_ccbar.cpp
	g++ -Wall -c J_data_comparison_inclusive_ccbar.cpp $$(root-config --glibs --cflags --libs)
	g++ J_data_comparison_inclusive_ccbar.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_data_comparison_inclusive_ccbar.exe

J_LHC_inclusive_ccbar_cross_section: J_LHC_inclusive_ccbar_cross_section.cpp
	g++ -Wall -c J_LHC_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ J_LHC_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_LHC_inclusive_ccbar_cross_section.exe
 

fit_data_comparison_inclusive_ccbar: fit_data_comparison_inclusive_ccbar.cpp
	g++ -Wall -c fit_data_comparison_inclusive_ccbar.cpp $$(root-config --glibs --cflags --libs)
	g++ fit_data_comparison_inclusive_ccbar.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o fit_data_comparison_inclusive_ccbar.exe

diffractive_inclusive_ratio: diffractive_inclusive_ratio.cpp
	g++ -Wall -c diffractive_inclusive_ratio.cpp $$(root-config --glibs --cflags --libs)
	g++ diffractive_inclusive_ratio.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o diffractive_inclusive_ratio.exe

dipole_amp_reader: dipole_amp_reader.cpp
	g++ dipole_amp_reader.cpp -o dipole_amp_reader.exe

test: test.cpp
	g++ test.cpp $$(root-config --glibs --cflags --libs) -o test.exe

plot_comparator: plot_comparator.cpp
	g++ plot_comparator.cpp $$(root-config --glibs --cflags --libs) -o plot_comparator.exe

J_diff_ccbar_cross_section: J_diff_ccbar_cross_section.cpp
	g++ -Wall -c J_diff_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ J_diff_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_diff_ccbar_cross_section.exe

J_differential_diffractive_sigma: J_differential_diffractive_sigma.cpp
	g++ -Wall -c J_differential_diffractive_sigma.cpp $$(root-config --glibs --cflags --libs)
	g++ J_differential_diffractive_sigma.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_differential_diffractive_sigma.exe

J_numerical_differential_diffractive_sigma: J_numerical_differential_diffractive_sigma.cpp
	g++ -Wall -c J_numerical_differential_diffractive_sigma.cpp $$(root-config --glibs --cflags --libs)
	g++ J_numerical_differential_diffractive_sigma.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_numerical_differential_diffractive_sigma.exe
	
J_simplified_numerical_differential_diffractive_sigma: J_simplified_numerical_differential_diffractive_sigma.cpp
	g++ -Wall -c J_simplified_numerical_differential_diffractive_sigma.cpp $$(root-config --glibs --cflags --libs)
	g++ J_simplified_numerical_differential_diffractive_sigma.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_simplified_numerical_differential_diffractive_sigma.exe

J_data_comparison_differential_diffractive: J_data_comparison_differential_diffractive.cpp
	g++ J_data_comparison_differential_diffractive.cpp $$(root-config --glibs --cflags --libs) -o J_data_comparison_differential_diffractive.exe

multiplotter: multiplotter.cpp
	g++ multiplotter.cpp $$(root-config --glibs --cflags --libs) -o multiplotter.exe

LHC_multiplotter: LHC_multiplotter.cpp
	g++ LHC_multiplotter.cpp $$(root-config --glibs --cflags --libs) -o LHC_multiplotter.exe

clean:
	-rm *.exe
	-rm *.o
