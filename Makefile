CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
INCLUDES := -Iinclude

SRC := src/gauss.cpp src/solve.cpp src/inverse_norm_estimate.cpp src/qr.cpp src/iterative.cpp
OBJ_DIR := build/obj
OBJ := $(SRC:src/%.cpp=$(OBJ_DIR)/%.o)

TEST_SRC := test/test.cpp
TIME_SRC := test/timeConsumingTest.cpp
HW2_SRC := test/homework02.cpp
HW3_SRC := test/homework03.cpp
HW5_SRC := test/homework05.cpp
HW8_SRC := test/homework08.cpp
HW10_SRC := test/homework10.cpp
HW11_SRC := test/homework11.cpp
NORM_SRC := test/test_inverse_norm.cpp
QR_SRC := test/test_qr.cpp

TEST_BIN := build/test
TIME_BIN := build/timeConsumingTest
HW2_BIN := build/homework02
HW3_BIN := build/homework03
HW5_BIN := build/homework05
HW8_BIN := build/homework08
HW10_BIN := build/homework10
HW11_BIN := build/homework11
NORM_BIN := build/test_inverse_norm
QR_BIN := build/test_qr

.PHONY: all test timeConsumingTest hw2 hw3 hw5 hw8 hw10 hw11 inverseNorm qr run-test run-timeConsumingTest run-hw2 run-hw3 run-hw5 run-hw8 run-hw10 run-hw11 run-inverseNorm run-qr clean

all: test timeConsumingTest hw2 hw3 hw5 hw8 hw10 hw11 inverseNorm qr

test: $(TEST_BIN)
timeConsumingTest: $(TIME_BIN)
hw2: $(HW2_BIN)
hw3: $(HW3_BIN)
hw5: $(HW5_BIN)
hw8: $(HW8_BIN)
hw10: $(HW10_BIN)
hw11: $(HW11_BIN)
inverseNorm: $(NORM_BIN)
qr: $(QR_BIN)

run-test: $(TEST_BIN)
	./$(TEST_BIN)

run-timeConsumingTest: $(TIME_BIN)
	./$(TIME_BIN) $(START_N) $(END_N) $(STEP)

run-hw2: $(HW2_BIN)
	./$(HW2_BIN)

run-hw3: $(HW3_BIN)
	./$(HW3_BIN)

run-hw5: $(HW5_BIN)
	./$(HW5_BIN)

run-hw8: $(HW8_BIN)
	./$(HW8_BIN)

run-hw10: $(HW10_BIN)
	./$(HW10_BIN)

run-hw11: $(HW11_BIN)
	./$(HW11_BIN)

run-inverseNorm: $(NORM_BIN)
	./$(NORM_BIN)

run-qr: $(QR_BIN)
	./$(QR_BIN)

# Defaults for timeConsumingTest arguments (override via make START_N=... END_N=... STEP=...)
START_N ?= 10
END_N ?= 500
STEP ?= 10

$(TEST_BIN): $(OBJ) $(TEST_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(TEST_SRC) -o $@

$(TIME_BIN): $(OBJ) $(TIME_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(TIME_SRC) -o $@

$(HW2_BIN): $(OBJ) $(HW2_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(HW2_SRC) -o $@

$(HW3_BIN): $(OBJ) $(HW3_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(HW3_SRC) -o $@

$(HW5_BIN): $(OBJ) $(HW5_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(HW5_SRC) -o $@

$(HW8_BIN): $(OBJ) $(HW8_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(HW8_SRC) -o $@

$(HW10_BIN): $(OBJ) $(HW10_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(HW10_SRC) -o $@

$(HW11_BIN): $(OBJ) $(HW11_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(HW11_SRC) -o $@

$(NORM_BIN): $(OBJ) $(NORM_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(NORM_SRC) -o $@

$(QR_BIN): $(OBJ) $(QR_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(QR_SRC) -o $@

$(OBJ_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf build
