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

ez: all
	mpirun -n 4 $(OUT) ./instances/kp_100_1.in
	mpirun -n 4 $(OUT) ./instances/kp_100_2.in
	mpirun -n 4 $(OUT) ./instances/kp_1000_1.in
	mpirun -n 4 $(OUT) ./instances/kp_1000_2.in
	mpirun -n 4 $(OUT) ./instances/kp_10000_1_0.01.in
	mpirun -n 4 $(OUT) ./instances/kp_10000_2_0.01.in

hard: all
	mpirun -n 4 $(OUT) ./instances/kp_10000_1_0.04.in
	mpirun -n 4 $(OUT) ./instances/kp_10000_1_0.05.in
	mpirun -n 4 $(OUT) ./instances/kp_10000_1_0.1.in
	mpirun -n 4 $(OUT) ./instances/kp_10000_2_0.04.in
	mpirun -n 4 $(OUT) ./instances/kp_10000_2_0.05.in
	mpirun -n 4 $(OUT) ./instances/kp_10000_2_0.1.in

