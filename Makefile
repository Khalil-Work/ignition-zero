CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic

BUILD := build

.PHONY: all test clean

all: $(BUILD)/ignition_demo $(BUILD)/test_ignition

$(BUILD):
	@if not exist $(BUILD) mkdir $(BUILD)

$(BUILD)/ignition_demo: src/ignition_demo.cpp src/ignition_zero.cpp src/ignition_zero.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) src/ignition_demo.cpp src/ignition_zero.cpp -o $@

$(BUILD)/test_ignition: tests/test_ignition.cpp src/ignition_zero.cpp src/ignition_zero.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) tests/test_ignition.cpp src/ignition_zero.cpp -o $@

test: $(BUILD)/test_ignition
	$(BUILD)\test_ignition.exe

clean:
	@if exist $(BUILD) rmdir /s /q $(BUILD)
