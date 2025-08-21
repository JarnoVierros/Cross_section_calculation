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
	make ratio_plotter
	make LHC_ratio_plotter
	make nuclear_suppression_ratio_plotter
	make dipole_amp_reader
	make plot_comparator
	make J_diff_ccbar_cross_section
	make J_diff_Q2_ccbar_cross_section
	make J_differential_diffractive_sigma
	make J_numerical_differential_diffractive_sigma
	make J_simplified_differential_diffractive_sigma
	make multiplotter
	make J_LHC_inclusive_ccbar_cross_section
	make J_improvised_inclusive_data_fit
	make J_light_quark_sigma_calculator
	make direct_light_quark_sigma_calculator
	make J_LHC_exclusive_ccbar_cross_section
	make variable_change_test
	make J_LHC_simplified_exclusive_ccbar_cross_section
	make J_LHC_Q2_int_simplified_exclusive_ccbar_cross_section
	make interpolation_test
	make GBW_LHC_diffractive_ccbar_cross_section
	make GBW_var_change_LHC_diffractive_ccbar_cross_section
	make integration_tests
	make dipole_amplitude_generator
	make direct_J_differential_diffractive_sigma
	make F2_light_sigma_calculator
	make direct_integrated_diffractive_sigma
	make double_plot
	make rapidity_LHC_inclusive_D0

Jb_inclusive_ccbar_cross_section: Jb_inclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c Jb_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ Jb_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o Jb_inclusive_ccbar_cross_section.exe
	
J_inclusive_ccbar_cross_section: J_inclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ J_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o J_inclusive_ccbar_cross_section.exe

J_inclusive_Q2_ccbar_cross_section: J_inclusive_Q2_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_inclusive_Q2_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ J_inclusive_Q2_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o J_inclusive_Q2_ccbar_cross_section.exe

inclusive_ccbar_cross_section: inclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o inclusive_ccbar_cross_section.exe

old_inclusive_ccbar_cross_section: old_inclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c old_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ old_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o old_inclusive_ccbar_cross_section.exe

complicated_exclusive_ccbar_cross_section: complicated_exclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c complicated_exclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ complicated_exclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o complicated_exclusive_ccbar_cross_section.exe

simplified_exclusive_ccbar_cross_section: simplified_exclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c simplified_exclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ simplified_exclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o simplified_exclusive_ccbar_cross_section.exe
	
alternative_exclusive_ccbar_cross_section: alternative_exclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c alternative_exclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ alternative_exclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o alternative_exclusive_ccbar_cross_section.exe
	
intQ2_inclusive_ccbar_cross_section: intQ2_inclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c intQ2_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ intQ2_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o intQ2_inclusive_ccbar_cross_section.exe

intx_inclusive_ccbar_cross_section: intx_inclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c intx_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ intx_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o intx_inclusive_ccbar_cross_section.exe

W_inclusive_ccbar_cross_section: W_inclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c W_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ W_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o W_inclusive_ccbar_cross_section.exe

reduced_inclusive_ccbar_cross_section: reduced_inclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c reduced_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ reduced_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lgsl -lgslcblas -lm -o reduced_inclusive_ccbar_cross_section.exe

parameter_fitter_inclusive_ccbar: parameter_fitter_inclusive_ccbar.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c parameter_fitter_inclusive_ccbar.cpp $$(root-config --glibs --cflags --libs)
	g++ parameter_fitter_inclusive_ccbar.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o parameter_fitter_inclusive_ccbar.exe

old_parameter_fitter_inclusive_ccbar: old_parameter_fitter_inclusive_ccbar.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c old_parameter_fitter_inclusive_ccbar.cpp $$(root-config --glibs --cflags --libs)
	g++ old_parameter_fitter_inclusive_ccbar.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o old_parameter_fitter_inclusive_ccbar.exe

J_data_comparison_inclusive_ccbar: J_data_comparison_inclusive_ccbar.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_data_comparison_inclusive_ccbar.cpp $$(root-config --glibs --cflags --libs)
	g++ J_data_comparison_inclusive_ccbar.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_data_comparison_inclusive_ccbar.exe

J_LHC_inclusive_ccbar_cross_section: J_LHC_inclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_LHC_inclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ J_LHC_inclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_LHC_inclusive_ccbar_cross_section.exe
 

fit_data_comparison_inclusive_ccbar: fit_data_comparison_inclusive_ccbar.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c fit_data_comparison_inclusive_ccbar.cpp $$(root-config --glibs --cflags --libs)
	g++ fit_data_comparison_inclusive_ccbar.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o fit_data_comparison_inclusive_ccbar.exe

ratio_plotter: ratio_plotter.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c ratio_plotter.cpp $$(root-config --glibs --cflags --libs)
	g++ ratio_plotter.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o ratio_plotter.exe

LHC_ratio_plotter: LHC_ratio_plotter.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c LHC_ratio_plotter.cpp $$(root-config --glibs --cflags --libs)
	g++ LHC_ratio_plotter.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o LHC_ratio_plotter.exe

nuclear_suppression_ratio_plotter: nuclear_suppression_ratio_plotter.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c nuclear_suppression_ratio_plotter.cpp $$(root-config --glibs --cflags --libs)
	g++ nuclear_suppression_ratio_plotter.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o nuclear_suppression_ratio_plotter.exe

dipole_amp_reader: dipole_amp_reader.cpp
	g++ dipole_amp_reader.cpp -o dipole_amp_reader.exe

test: test.cpp
	g++ test.cpp $$(root-config --glibs --cflags --libs) -o test.exe

plot_comparator: plot_comparator.cpp
	g++ plot_comparator.cpp $$(root-config --glibs --cflags --libs) -o plot_comparator.exe

rapidity_plot_comparator: rapidity_plot_comparator.cpp
	g++ rapidity_plot_comparator.cpp $$(root-config --glibs --cflags --libs) -o rapidity_plot_comparator.exe

J_diff_ccbar_cross_section: J_diff_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_diff_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ J_diff_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_diff_ccbar_cross_section.exe

J_diff_Q2_ccbar_cross_section: J_diff_Q2_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_diff_Q2_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ J_diff_Q2_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_diff_Q2_ccbar_cross_section.exe

J_differential_diffractive_sigma: J_differential_diffractive_sigma.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_differential_diffractive_sigma.cpp $$(root-config --glibs --cflags --libs)
	g++ J_differential_diffractive_sigma.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_differential_diffractive_sigma.exe

J_numerical_differential_diffractive_sigma: J_numerical_differential_diffractive_sigma.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_numerical_differential_diffractive_sigma.cpp $$(root-config --glibs --cflags --libs)
	g++ J_numerical_differential_diffractive_sigma.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_numerical_differential_diffractive_sigma.exe
	
J_simplified_numerical_differential_diffractive_sigma: J_simplified_numerical_differential_diffractive_sigma.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_simplified_numerical_differential_diffractive_sigma.cpp $$(root-config --glibs --cflags --libs)
	g++ J_simplified_numerical_differential_diffractive_sigma.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_simplified_numerical_differential_diffractive_sigma.exe

J_data_comparison_differential_diffractive: J_data_comparison_differential_diffractive.cpp
	g++ J_data_comparison_differential_diffractive.cpp $$(root-config --glibs --cflags --libs) -o J_data_comparison_differential_diffractive.exe

multiplotter: multiplotter.cpp
	g++ multiplotter.cpp $$(root-config --glibs --cflags --libs) -o multiplotter.exe

LHC_multiplotter: LHC_multiplotter.cpp
	g++ LHC_multiplotter.cpp $$(root-config --glibs --cflags --libs) -o LHC_multiplotter.exe

double_plot: double_plot.cpp
	g++ double_plot.cpp $$(root-config --glibs --cflags --libs) -o double_plot.exe

LHC_multiplotter_2: LHC_multiplotter_2.cpp
	g++ LHC_multiplotter_2.cpp $$(root-config --glibs --cflags --libs) -o LHC_multiplotter_2.exe

J_improvised_inclusive_data_fit: J_improvised_inclusive_data_fit.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_improvised_inclusive_data_fit.cpp $$(root-config --glibs --cflags --libs)
	g++ J_improvised_inclusive_data_fit.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_improvised_inclusive_data_fit.exe

custom_plot: custom_plot.cpp
	g++ custom_plot.cpp $$(root-config --glibs --cflags --libs) -o custom_plot.exe

J_light_quark_sigma_calculator: J_light_quark_sigma_calculator.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_light_quark_sigma_calculator.cpp $$(root-config --glibs --cflags --libs)
	g++ J_light_quark_sigma_calculator.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_light_quark_sigma_calculator.exe

direct_light_quark_sigma_calculator: direct_light_quark_sigma_calculator.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c direct_light_quark_sigma_calculator.cpp $$(root-config --glibs --cflags --libs)
	g++ direct_light_quark_sigma_calculator.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o direct_light_quark_sigma_calculator.exe

J_LHC_exclusive_ccbar_cross_section: J_LHC_exclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_LHC_exclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ J_LHC_exclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_LHC_exclusive_ccbar_cross_section.exe

variable_change_test: variable_change_test.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c variable_change_test.cpp $$(root-config --glibs --cflags --libs)
	g++ variable_change_test.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o variable_change_test.exe

J_LHC_simplified_exclusive_ccbar_cross_section: J_LHC_simplified_exclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_LHC_simplified_exclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ J_LHC_simplified_exclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_LHC_simplified_exclusive_ccbar_cross_section.exe

J_LHC_Q2_int_simplified_exclusive_ccbar_cross_section: J_LHC_Q2_int_simplified_exclusive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c J_LHC_Q2_int_simplified_exclusive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ J_LHC_Q2_int_simplified_exclusive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o J_LHC_Q2_int_simplified_exclusive_ccbar_cross_section.exe

interpolation_test: interpolation_test.cpp
	g++ interpolation_test.cpp $$(root-config --glibs --cflags --libs) -o interpolation_test.exe

GBW_LHC_diffractive_ccbar_cross_section: GBW_LHC_diffractive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c GBW_LHC_diffractive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ GBW_LHC_diffractive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o GBW_LHC_diffractive_ccbar_cross_section.exe

GBW_var_change_LHC_diffractive_ccbar_cross_section: GBW_var_change_LHC_diffractive_ccbar_cross_section.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c GBW_var_change_LHC_diffractive_ccbar_cross_section.cpp $$(root-config --glibs --cflags --libs)
	g++ GBW_var_change_LHC_diffractive_ccbar_cross_section.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o GBW_var_change_LHC_diffractive_ccbar_cross_section.exe

integration_tests: integration_tests.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c integration_tests.cpp $$(root-config --glibs --cflags --libs)
	g++ integration_tests.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o integration_tests.exe

dipole_amplitude_generator: dipole_amplitude_generator.cpp
	g++ dipole_amplitude_generator.cpp -o dipole_amplitude_generator.exe

direct_J_differential_diffractive_sigma: direct_J_differential_diffractive_sigma.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c direct_J_differential_diffractive_sigma.cpp $$(root-config --glibs --cflags --libs)
	g++ direct_J_differential_diffractive_sigma.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o direct_J_differential_diffractive_sigma.exe	

rapidity_LHC_inclusive_D0: rapidity_LHC_inclusive_D0.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c rapidity_LHC_inclusive_D0.cpp $$(root-config --glibs --cflags --libs)
	g++ rapidity_LHC_inclusive_D0.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o rapidity_LHC_inclusive_D0.exe
	
rapidity_xyxbar_LHC_inclusive_D0: rapidity_xyxbar_LHC_inclusive_D0.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c rapidity_xyxbar_LHC_inclusive_D0.cpp $$(root-config --glibs --cflags --libs)
	g++ rapidity_xyxbar_LHC_inclusive_D0.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o rapidity_xyxbar_LHC_inclusive_D0.exe

rapidity_polar_xyxbar_LHC_inclusive_D0: rapidity_polar_xyxbar_LHC_inclusive_D0.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c rapidity_polar_xyxbar_LHC_inclusive_D0.cpp $$(root-config --glibs --cflags --libs)
	g++ rapidity_polar_xyxbar_LHC_inclusive_D0.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o rapidity_polar_xyxbar_LHC_inclusive_D0.exe

F2_light_sigma_calculator: F2_light_sigma_calculator.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c F2_light_sigma_calculator.cpp $$(root-config --glibs --cflags --libs)
	g++ F2_light_sigma_calculator.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o F2_light_sigma_calculator.exe

direct_integrated_diffractive_sigma: direct_integrated_diffractive_sigma.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c direct_integrated_diffractive_sigma.cpp $$(root-config --glibs --cflags --libs)
	g++ direct_integrated_diffractive_sigma.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o direct_integrated_diffractive_sigma.exe

qqg_contribution_calculator: qqg_contribution_calculator.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c qqg_contribution_calculator.cpp $$(root-config --glibs --cflags --libs)
	g++ qqg_contribution_calculator.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o qqg_contribution_calculator.exe
	
Jani_qqg_contribution_calculator: Jani_qqg_contribution_calculator.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c Jani_qqg_contribution_calculator.cpp $$(root-config --glibs --cflags --libs)
	g++ Jani_qqg_contribution_calculator.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o Jani_qqg_contribution_calculator.exe
	
sandbox: sandbox.cpp
	g++ -Wall -Wno-sign-compare -Wno-unused-but-set-variable -c sandbox.cpp $$(root-config --glibs --cflags --libs)
	g++ sandbox.o $$(root-config --glibs --cflags --libs) -lMinuit -lgsl -lgslcblas -lm -o sandbox.exe

linterp_test: linterp_test.cpp
	g++ linterp_test.cpp -o linterp_test.exe

clean:
	-rm *.exe
	-rm *.o
