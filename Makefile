FLAGS = -lm -fopenmp -lrt `pkg-config --cflags --libs json-glib-1.0`

BIN		:= bin
SRC		:= src
INCLUDE	:= include

EXECUTABLE	:= uwsim

main : $(EXECUTABLE)

all : $(EXECUTABLE) $(EXECUTABLE)_dbg

install : /usr/local/bin/$(EXECUTABLE)

/usr/local/bin/$(EXECUTABLE) : $(EXECUTABLE)
	@if [ $(USER) = "root" ]; then\
		cp $(EXECUTABLE) /usr/local/bin/;\
	else\
		echo "You must be root";\
	fi

uninstall :
	@if [ $(USER) = "root" ]; then\
		rm -fv /usr/local/bin/$(EXECUTABLE);\
	else\
		echo "You must be root";\
	fi

$(EXECUTABLE) : $(SRC)/*.c
	gcc -Wall -O1 $^ -o $(BIN)/$@ $(FLAGS) -I$(INCLUDE)

debug : $(EXECUTABLE)_dbg

$(EXECUTABLE)_dbg : $(SRC)/*.c
	gcc -Wall -g $^ -o $(BIN)/$@ $(FLAGS) -I$(INCLUDE)

clean :
	@rm -fv $(BIN)/$(EXECUTABLE) $(BIN)/$(EXECUTABLE)_dbg
