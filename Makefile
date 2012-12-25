CC=gcc
CFLAGS=-Wall -Wextra -Werror -std=c99 -pedantic

linalg: 
	$(CC) $(CFLAGS) -o linalg util.c matrix.c main.c