
# ------------------------------------------------------------------------------------------------------
# This is a makefile. The Code can be compiled and run by using the command "make" inside the directory.
# ------------------------------------------------------------------------------------------------------

CC     = gcc				# CC 		is the compiler used
CFLAGS = -pedantic -Wall 	# CFLAGS 	sets the compiler flags to compile 

made: main.c RK_procedure.c solve_simp_pend.c solve_doub_pend.c solve_trip_pend.c tools.c header.h							# the files to compile
	$(CC) -o ./output main.c RK_procedure.c solve_simp_pend.c solve_doub_pend.c solve_trip_pend.c tools.c header.h -lm
	./output 

	
clean:
	rm -v ./output
	rm -v ../data/*