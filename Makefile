CC=gcc
CFLAGS=-Wall -Wextra -Werror -std=c99 -pedantic

all: 
	$(CC) $(CFLAGS) -o linalg util.c matrix.c main.c

clean:
	rm -f linalg