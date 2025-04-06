OBJ_DIR = build
SRC_DIR = src

CXX = mpic++
CXXFLAGS = -O3 -MMD -MP -std=c++17 -Wall -Wextra

SRC_CPP      := $(shell find $(SRC_DIR) -type f -name '*.cpp')
OBJ_CPP      := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_CPP))

OUT = kp

-include $(patsubst %.o, %.d, $(OBJ_CPP))

.PHONY: clean

all: $(OBJ_CPP)
	$(CXX) $(OBJ_CPP) -o $(OUT)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp Makefile
	@mkdir -p $(shell dirname $@)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	rm -rf $(OBJ_DIR)
	rm $(OUT)

# Cibles pour les tests "faciles" avec les différents solvers
ez1: all
	mpirun -n 4 $(OUT) 1 ./instances/kp_100_1.in
	mpirun -n 4 $(OUT) 1 ./instances/kp_100_2.in
	mpirun -n 4 $(OUT) 1 ./instances/kp_1000_1.in
	mpirun -n 4 $(OUT) 1 ./instances/kp_1000_2.in
	mpirun -n 4 $(OUT) 1 ./instances/kp_10000_1_0.01.in
	mpirun -n 4 $(OUT) 1 ./instances/kp_10000_2_0.01.in

ez2: all
	mpirun -n 4 $(OUT) 2 ./instances/kp_100_1.in
	mpirun -n 4 $(OUT) 2 ./instances/kp_100_2.in
	mpirun -n 4 $(OUT) 2 ./instances/kp_1000_1.in
	mpirun -n 4 $(OUT) 2 ./instances/kp_1000_2.in
	mpirun -n 4 $(OUT) 2 ./instances/kp_10000_1_0.01.in
	mpirun -n 4 $(OUT) 2 ./instances/kp_10000_2_0.01.in

ez3: all
	mpirun -n 4 $(OUT) 3 ./instances/kp_100_1.in
	mpirun -n 4 $(OUT) 3 ./instances/kp_100_2.in
	mpirun -n 4 $(OUT) 3 ./instances/kp_1000_1.in
	mpirun -n 4 $(OUT) 3 ./instances/kp_1000_2.in
	mpirun -n 4 $(OUT) 3 ./instances/kp_10000_1_0.01.in
	mpirun -n 4 $(OUT) 3 ./instances/kp_10000_2_0.01.in

# Cibles pour les 4ests "difficiles" avec les différents solvers
hard1: all
	mpirun -n 4 $(OUT) 1 ./instances/kp_10000_1_0.04.in
	mpirun -n 4 $(OUT) 1 ./instances/kp_10000_1_0.05.in
	mpirun -n 4 $(OUT) 1 ./instances/kp_10000_1_0.1.in
	mpirun -n 4 $(OUT) 1 ./instances/kp_10000_2_0.04.in
	mpirun -n 4 $(OUT) 1 ./instances/kp_10000_2_0.05.in
	mpirun -n 4 $(OUT) 1 ./instances/kp_10000_2_0.1.in

hard2: all
	mpirun -n 4 $(OUT) 2 ./instances/kp_10000_1_0.04.in
	mpirun -n 4 $(OUT) 2 ./instances/kp_10000_1_0.05.in
	mpirun -n 4 $(OUT) 2 ./instances/kp_10000_1_0.1.in
	mpirun -n 4 $(OUT) 2 ./instances/kp_10000_2_0.04.in
	mpirun -n 4 $(OUT) 2 ./instances/kp_10000_2_0.05.in
	mpirun -n 4 $(OUT) 2 ./instances/kp_10000_2_0.1.in

hard3: all
	mpirun -n 4 $(OUT) 3 ./instances/kp_10000_1_0.04.in
	mpirun -n 4 $(OUT) 3 ./instances/kp_10000_1_0.05.in
	mpirun -n 4 $(OUT) 3 ./instances/kp_10000_1_0.1.in
	mpirun -n 4 $(OUT) 3 ./instances/kp_10000_2_0.04.in
	mpirun -n 4 $(OUT) 3 ./instances/kp_10000_2_0.05.in
	mpirun -n 4 $(OUT) 3 ./instances/kp_10000_2_0.1.in

