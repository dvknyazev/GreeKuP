version=GreeKuP_v1.0

CC=mpiCC

#MYHOME should contain path to the folder, where $(DIRNAME) is located, e.g., MYHOME=/home/knyazev/GreeKuP
#MYHOME should not contain / at the end
MYHOME=$(GREEKUP_HOME)
DIRNAME=GreeKuP_v1.0

MOD=$(MYHOME)/$(DIRNAME)/modules_c++/
PROG=$(MYHOME)/$(DIRNAME)/condcalc_c++/
OMEGA=$(MYHOME)/$(DIRNAME)/omega_generation/

ifeq ($(MYHOME),)
$(error The variable MYHOME is not set. Please set this variable to proceed. Execution of the Makefile is terminated)
endif

obj_simple= lattice.o keys.o omegainout.o mymath.o
obj_module= myinout.o nablainout.o $(obj_simple)

progs_simple= condmultipleprocess dccalc heatcalc

progs_omega= omega_generate printomega

programs_condcalc= GreeKuP printnabla \
		$(progs_simple)
programs_omega= $(progs_omega)




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
