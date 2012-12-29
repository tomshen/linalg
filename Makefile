CC=gcc
CFLAGS=-Wall -Wextra -Werror -std=c99 -pedantic

all: 
	$(CC) $(CFLAGS) -o linalg util.c matrix_classification.c matrix_creation.c matrix_util.c vectors.c matrix_basic.c matrix_applications.c main.c

clean:
	rm -f linalg