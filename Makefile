OBJ_DIR = build
SRC_DIR = src

CXX = clang++
CXXFLAGS = -O3 -MMD -MP

.PHONY: clean

all: 1

1: $(OBJ_DIR)/kp.o
	$(CXX) $(OBJ_DIR)/kp.o -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp Makefile
	@mkdir -p $(shell dirname $@)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	rm -rf $(OBJ_DIR)
	rm 1
