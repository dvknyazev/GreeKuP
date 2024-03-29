version=GreeKuP_v1.0

CC=$(GREEKUP_COMPILER)
#Set the C++ compiler available on your machine. It will be used to compile GreeKuP


ifeq ($(CC),)
$(error The variable CC, that specifies the C++ compiler for GreeKuP is not set. Please set this variable to proceed. You may set the environment variable GREEKUP_COMPILER or set CC manually in the Makefile. Execution of the Makefile is terminated)
endif

makefile_dir:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
#This_command returns absolute path to the current Makefile. This is crucial if this Makefile is run from other directory

MOD:=$(makefile_dir)/modules_c++/
PROG:=$(makefile_dir)/condcalc_c++/
OMEGA:=$(makefile_dir)/omega_generation/


obj_simple:= lattice.o keys.o omegainout.o mymath.o
obj_module:= myinout.o nablainout.o $(obj_simple)

progs_simple:= condmultipleprocess dccalc heatcalc

progs_omega:= omega_generate printomega

programs_condcalc:= GreeKuP printnabla \
		$(progs_simple)
programs_omega:= $(progs_omega)




all: $(addprefix $(PROG), $(programs_condcalc)) \
	$(addprefix $(OMEGA), $(programs_omega)) 

$(programs_condcalc): %: $(PROG)%
$(programs_omega): %: $(OMEGA)%

$(obj_module): %: $(MOD)%



$(PROG)GreeKuP: $(PROG)GreeKuP.o $(MOD)myinout.o $(MOD)keys.o $(MOD)omegainout.o $(MOD)mymath.o
	$(CC) $(PROG)GreeKuP.o $(MOD)myinout.o $(MOD)keys.o $(MOD)omegainout.o $(MOD)mymath.o -o $(PROG)GreeKuP

$(PROG)printnabla: $(PROG)printnabla.o $(MOD)keys.o $(MOD)lattice.o $(MOD)nablainout.o $(MOD)mymath.o
	$(CC) $(PROG)printnabla.o $(MOD)keys.o $(MOD)lattice.o $(MOD)nablainout.o $(MOD)mymath.o -o $(PROG)printnabla

$(addprefix $(PROG), $(progs_simple)): %: %.o $(MOD)myinout.o $(MOD)keys.o
	$(CC) $@.o $(MOD)myinout.o $(MOD)keys.o -o $@



$(addprefix $(OMEGA), $(progs_omega)): %: %.o $(MOD)keys.o $(MOD)omegainout.o
	$(CC) $@.o $(MOD)keys.o $(MOD)omegainout.o -o $@



$(PROG)GreeKuP.o: $(PROG)GreeKuP.C $(MOD)myinout.h $(MOD)keys.h $(MOD)omegainout.h
	$(CC) -c $(PROG)GreeKuP.C -I$(MOD) -o $(PROG)GreeKuP.o

$(PROG)printnabla.o: $(PROG)printnabla.C $(MOD)keys.h $(MOD)lattice.h $(MOD)nablainout.h $(MOD)myinout.h
	$(CC) -c $(PROG)printnabla.C -I$(MOD) -o $(PROG)printnabla.o

$(addsuffix .o, $(addprefix $(PROG), $(progs_simple))): %.o: %.C $(MOD)myinout.h $(MOD)keys.h
	$(CC) -c $*.C -I$(MOD) -o $*.o




$(addsuffix .o, $(addprefix $(OMEGA), $(progs_omega))): %.o: %.C $(MOD)keys.h $(MOD)omegainout.h
	$(CC) -c $*.C -I$(MOD) -o $*.o



$(MOD)myinout.o: $(MOD)myinout.C $(MOD)myinout.h
	$(CC) -c $(MOD)myinout.C -o $(MOD)myinout.o

$(MOD)nablainout.o: $(MOD)nablainout.C $(MOD)nablainout.h $(MOD)myinout.h $(MOD)mymath.h
	$(CC) -c $(MOD)nablainout.C -I$(MOD) -o $(MOD)nablainout.o

$(addprefix $(MOD), $(obj_simple)): %.o: %.C %.h
	$(CC) -c $*.C -o $*.o

clean:
	cd $(MOD); rm -f *.o; cd -
	cd $(PROG); rm -f *.o; rm -f $(programs_condcalc); cd -
	cd $(OMEGA); rm -f *.o; rm -f $(programs_omega); cd -
