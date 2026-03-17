all:
	make direct_light_quark_sigma_calculator
	make Phi_generator
	make Phi_generator_polar
	make Phi_generator_jv
	make Phi_plotter

clean:
	-rm -f *.exe
	-rm -f *.o

direct_light_quark_sigma_calculator: direct_light_quark_sigma_calculator.cpp
	g++ direct_light_quark_sigma_calculator.cpp $$(root-config --glibs --cflags --libs) -lgsl -o direct_light_quark_sigma_calculator.exe


single_plot: single_plot.cpp
	g++ single_plot.cpp $$(root-config --glibs --cflags --libs) -lgsl -o single_plot.exe

J_data_comparison_differential_diffractive: J_data_comparison_differential_diffractive.cpp
	g++ J_data_comparison_differential_diffractive.cpp $$(root-config --glibs --cflags --libs) -lgsl -o J_data_comparison_differential_diffractive.exe

Phi_generator: Phi_generator.cpp
	g++ Phi_generator.cpp $$(root-config --glibs --cflags --libs) -lgsl -o Phi_generator.exe

Phi_generator_polar: Phi_generator_polar.cpp
	g++ Phi_generator_polar.cpp $$(root-config --glibs --cflags --libs) -lgsl -o Phi_generator_polar.exe

Phi_generator_jv: Phi_generator_jv.cpp
	g++ Phi_generator_jv.cpp $$(root-config --glibs --cflags --libs) -lgsl -o Phi_generator_jv.exe

Phi_plotter: Phi_plotter.cpp
	g++ Phi_plotter.cpp $$(root-config --glibs --cflags --libs) -lgsl -o Phi_plotter.exe
