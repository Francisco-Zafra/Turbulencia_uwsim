HDR = src/settings.h src/simulation.h src/simulation_file.h src/results.h src/photon_block.h src/utilities.h src/status.h src/urand.h src/processing.h src/configuration.h
BDY = src/settings.c src/simulation.c src/simulation_file.c src/results.c src/photon_block.c src/utilities.c src/status.c src/urand.c src/processing.c src/main.c
SRC = $(BDY) $(HDR)

FLAGS = -lm -fopenmp -lrt `pkg-config --cflags --libs json-glib-1.0`

main : uwsim

all : uwsim uwsim_dbg

install : /usr/local/bin/uwsim

/usr/local/bin/uwsim : uwsim
	@if [ $(USER) = "root" ]; then\
		cp uwsim /usr/local/bin/;\
	else\
		echo "You must be root";\
	fi

uninstall :
	@if [ $(USER) = "root" ]; then\
		rm -fv /usr/local/bin/uwsim;\
	else\
		echo "You must be root";\
	fi

uwsim : $(SRC)
	gcc -Wall -O1 -o uwsim $(BDY) $(FLAGS)

debug : uwsim_dbg

uwsim_dbg : $(SRC)
	gcc -Wall -g -o uwsim_dbg $(BDY) $(FLAGS)

clean :
	@rm -fv uwsim uwsim_dbg
