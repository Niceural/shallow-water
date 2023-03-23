# https://www.partow.net/programming/makefile/index.html

CXX      := g++ # -c++
CXXFLAGS := -std=c++11 -pedantic-errors -Wall -Wextra -Wshadow -fopenmp # -Werror
LDFLAGS  := -L/usr/lib -lstdc++ -lm -lblas -lboost_program_options -lboost_timer -fopenmp
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/apps
TARGET   := shallowWater
INCLUDE  := -Iinclude/
SRC      :=                      \
   $(wildcard src/matrices/*.cpp) \
   $(wildcard src/*.cpp)         \

TESTS := $(wildcard tests/*.cpp)

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES := $(OBJECTS:.o=.d)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -MMD -o $@

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

-include $(DEPENDENCIES)

.PHONY: all build clean debug release info test test1 test2 test3 test4

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O3
release: all

analyze: CXXFLAGS += -O3 -DDEBUG -g
analyze: all
	rm -rf ./analyze.er
	OMP_NUM_THREADS=8 collect -o analyze.er ./build/apps/shallowWater --T 100.0 --dt 0.1 --Nx 500 --Ny 500 --ic 3

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

info:
	@echo "[*] Application dir: ${APP_DIR}     "
	@echo "[*] Object dir:      ${OBJ_DIR}     "
	@echo "[*] Sources:         ${SRC}         "
	@echo "[*] Objects:         ${OBJECTS}     "
	@echo "[*] Dependencies:    ${DEPENDENCIES}"

test: release
	./build/apps/shallowWater --dt 0.1 --T 87.0 --Nx 359 --Ny 100 --ic 3

test1: release
	./build/apps/shallowWater --dt 0.1 --T 80.0 --Nx 100 --Ny 100 --ic 1

test2: release
	./build/apps/shallowWater --dt 0.1 --T 80.0 --Nx 100 --Ny 100 --ic 2

test3: release
	./build/apps/shallowWater --dt 0.1 --T 80.0 --Nx 100 --Ny 100 --ic 3

test4: release
	./build/apps/shallowWater --dt 0.1 --T 80.0 --Nx 100 --Ny 100 --ic 4
