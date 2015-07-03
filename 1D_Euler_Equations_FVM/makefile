EXEC := 1D_Euler_solver
EXEC_DIRECTORY := ./
EXECUTABLES := 1D_Euler_solver

BASE_PATH := ../
BASIC_PATH := 1D_Euler_Equations_FVM/
Plot_PATH:=1D_Codes/
IPATH  := MYLIBRARY/include/
SPATH  := MYLIBRARY/src/
LPATH  := ../
IFLAGS := -I $(IPATH)  -I$(LPATH)
LFLAGS :=  #-lgsl -lgslcblas -lm -l fnccheck -L$(CFDRC)/lib -lDTF
DFLAGS := -ggdb3  -g
OFLAGS :=  -O3  -O2 -O
WFLAGS := -Wall -Wno-long-long
EFLAGS := -pg
AFLAGS := # -m64
CFLAGS := $(EFLAGS) -c $(AFLAGS) $(IFLAGS) $(DFLAGS) $(OFLAGS) $(WFLAGS)
CC     := g++ # i586-mingw32msvc-g++
STRIP  := strip # i586-mingw32msvc-strip

OBJS          =   	$(BASE_PATH)$(BASIC_PATH)Initialize.o\
			$(BASE_PATH)$(BASIC_PATH)Test_Cases.o\
			$(BASE_PATH)$(BASIC_PATH)Limiters.o\
			$(BASE_PATH)$(BASIC_PATH)Primitive.o\
			$(BASE_PATH)$(BASIC_PATH)LLF.o\
			$(BASE_PATH)$(BASIC_PATH)PVU.o\
			$(BASE_PATH)$(BASIC_PATH)PVU_P.o\
			$(BASE_PATH)$(BASIC_PATH)KFVS.o\
			$(BASE_PATH)$(BASIC_PATH)Steger_Warming.o\
			$(BASE_PATH)$(BASIC_PATH)Vanleer.o\
			$(BASE_PATH)$(BASIC_PATH)AUSM.o\
			$(BASE_PATH)$(BASIC_PATH)Roe.o\
			$(BASE_PATH)$(BASIC_PATH)HLL.o\
			$(BASE_PATH)$(BASIC_PATH)HLLC.o\
			$(BASE_PATH)$(BASIC_PATH)Read_Input.o\
			$(BASE_PATH)$(BASIC_PATH)Boundary_Conditions.o\
			$(BASE_PATH)$(BASIC_PATH)Output.o\
			$(BASE_PATH)$(BASIC_PATH)Solver.o\
			$(BASE_PATH)$(BASIC_PATH)Plot.o\
			$(BASE_PATH)$(BASIC_PATH)gnuplot_i.o\
			$(BASE_PATH)$(BASIC_PATH)Euler_Main.o

default:@echo "Type one of the following "
		@echo "  make  pc"
		@echo "  make  pc_debug"
		@echo "  make  pc_static"
		@echo "  make  pc_64"
		@echo "  make  pc_64_debug"
		@echo "  make  pc_64_static"

$(EXEC): $(OBJS)
		$(CC)  $(EFLAGS) $(AFLAGS) $(OBJS)  $(LFLAGS) -o $(EXEC_DIRECTORY)$@

%.o: %.cc
		$(CC) $(CFLAGS) `pwd`/$< -o $@
		@-makedepend -I$(IPATH) -I$(LPATH) -a -fMake_dependencies $< > /dev/null 2>&1

%.o: %.cpp
		$(CC)$(CFLAGS) `pwd`/$< -o $@
		@-makedepend -I$(IPATH) -I$(LPATH) -a -fMake_dependencies $< > /dev/null 2>&1

clean:
		rm -f $(OBJS)
		rm -f $(EXEC_DIRECTORY)$(EXEC)
		@-touch Make_dependencies
		@-makedepend -fMake_dependencies __arbitrary_file > /dev/null 2>&1
archive:
		ar rcs $(EXEC_DIRECTORY)$(EXEC).a $(OBJS)

#linux 32 bit

pc:
	@for i in ${EXECUTABLES}; do make  EXEC="$$i" "$$i" \
	"AFLAGS=-m32" "OFLAGS=-O3"; done
pc_debug:
		@for i in ${EXECUTABLES}; do make  EXEC="$$i" "$$i" \
		"AFLAGS=-m32" "DFLAGS=-ggdb3"; done
pc_static:
		@for i in ${EXECUTABLES}; do make  EXEC="$$i" "$$i" \
		"AFLAGS=-m32" "OFLAGS=-O2" "LFLAGS=-static"; done
pc64:
		@for i in ${EXECUTABLES}; do make  EXEC="$$i" "$$i" \
		"AFLAGS=-m64" "OFLAGS=-O3"; done
pc64_debug:
		@for i in ${EXECUTABLES}; do make  EXEC="$$i" "$$i" \
		"AFLAGS=-m64" "DFLAGS=-ggdb3"; done
pc64_static:
		@for i in ${EXECUTABLES}; do make  EXEC="$$i" "$$i" \
		"AFLAGS=-m64" "OFLAGS=-O2" "LFLAGS=-static"; done
-include Make_dependencies