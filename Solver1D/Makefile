CC = gcc
CFLAGS = -Wall -Werror
LIBS = -lm

SRCS = intMain.c EOS_GGL.c read_par_ic.c read_par_eos.c auxFun.c hugoniot_curves.c integral_curves.c mixed_curves.c principalLoops.c
OBJS = $(SRCS:.c=.o)
EXECUTABLE = exactRP

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(EXECUTABLE) $(OBJS)
