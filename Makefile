CC = g++

CLASSDIR = /home/cristinel/power_systems/pronerds
INCDIRS = $(CLASSDIR)/include

LIB_DIR = -L/usr/lib/X11
LIB = -lX11 -lm
X11_INCLUDE = -I/usr/X11R6/include

WARN_FLAGS = -Wall -Wpointer-arith -Wcast-qual -Wstrict-prototypes -O -D__USE_FIXED_PROTOTYPES__ -ansi -pedantic -Wmissing-prototypes -Wshadow -Wcast-align -D_POSIX_SOURCE
#WARN_FLAGS =
#DEBUG_FLAGS = -g
DEBUG_FLAGS =
OPT_FLAGS = -O3
#OPT_FLAGS =
# LDFLAGS = -lpthread
LDFLAGS =

FLAGS = $(OPT_FLAGS)
FLAGS += $(DEBUG_FLAGS)
FLAGS += $(addprefix -I, $(INCDIRS))

EXE = nerds

OBJ = nerds_utils.o nerds.o nerds_main.o nerds_gui.o nerds_mcmf.o 
SRC = nerds_utils.cpp nerds.cpp nerds_main.cpp nerds_gui.cpp nerds_mcmf.cpp
H = include/nerds_utils.h include/nerds.h include/nerds_gui.h include/nerds_mcmf.h 

$(EXE): $(OBJ)
	$(CC) $(FLAGS) $(OBJ) -o $(EXE) $(LIB_DIR) $(LIB) $(LDFLAGS)

nerds_utils.o: nerds_utils.cpp $(H)
	$(CC) -c $(FLAGS) nerds_utils.cpp

nerds.o: nerds.cpp $(H)
	$(CC) -c $(FLAGS) nerds.cpp

nerds_main.o: nerds_main.cpp $(H)
	$(CC) -c $(FLAGS) nerds_main.cpp

nerds_gui.o: nerds_gui.cpp $(H)
	$(CC) -c $(FLAGS) $(X11_INCLUDE) nerds_gui.cpp

nerds_mcmf.o: nerds_mcmf.cpp $(H)
	$(CC) -c $(FLAGS) nerds_mcmf.cpp
