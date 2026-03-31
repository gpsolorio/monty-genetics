# Makefile for het-vs-hom

LIB = -lm
CC = gcc
CFLAGS = -O2 -Wall -Werror

APP = hvh
SRC = hvh.c
OBJ = hvh.o

# Targets

default:
	make $(APP)

$(APP): $(OBJ)
	$(CC) -o $(APP) $(OBJ) $(LIB)

clean:
	rm -f *.o $(APP)

# Inference Rules

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

