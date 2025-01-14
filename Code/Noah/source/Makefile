
# CC legt den verwendeten Compiler fest.
# CFLAGS erlaubt die Verwendung von Compiler-Flags

CC     = gcc
CFLAGS = -pedantic -Wall 

# Die Module werden zusammen kompiliert, in "made" unter "out" abgelegt und direkt ausgeführt.

made: main.c my_RK4.c solve_DGLs.c DGLs.c tools.c header.h
	$(CC) -o ./out/made main.c my_RK4.c solve_DGLs.c DGLs.c tools.c header.h -lm		
	./out/made 
	
clean:
	rm -v ./out/made