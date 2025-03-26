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
