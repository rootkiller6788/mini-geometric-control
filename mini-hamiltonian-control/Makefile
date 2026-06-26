# Makefile for mini-hamiltonian-control
CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -Wno-misleading-indentation -O2 -g -Iinclude
LDFLAGS  = -lm

SRC = src/hamiltonian_control.c       src/symplectic_geometry.c       src/pontryagin_maximum.c       src/port_hamiltonian.c       src/energy_shaping.c       src/symplectic_integrator.c       src/hjb_equation.c       src/poisson_geometry.c

OBJ = $(SRC:.c=.o)
TEST_SRC = tests/test_hamiltonian.c
TEST_BIN = test_hamiltonian
EXAMPLES = examples/example_harmonic_oscillator            examples/example_rigid_body            examples/example_dc_motor_ph
LIB = libhamiltonian_control.a

.PHONY: all clean test examples

all: $(LIB) test examples

$(LIB): $(OBJ)
	ar rcs $@ $^

src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(TEST_SRC) $(LIB)
	$(CC) $(CFLAGS) $< -L. -lhamiltonian_control $(LDFLAGS) -o $@

examples: $(EXAMPLES)

examples/example_harmonic_oscillator: examples/example_harmonic_oscillator.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lhamiltonian_control $(LDFLAGS) -o $@

examples/example_rigid_body: examples/example_rigid_body.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lhamiltonian_control $(LDFLAGS) -o $@

examples/example_dc_motor_ph: examples/example_dc_motor_ph.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lhamiltonian_control $(LDFLAGS) -o $@

clean:
	rm -f $(OBJ) $(LIB) $(TEST_BIN) $(EXAMPLES)
