default: build

build: ShallowWater

run:

test: test_ic

# src

ShallowWater.cpp:
	g++ -Wall -Wshadow -O2 -o target/ShallowWater.o -c src/ShallowWater.cpp

# tests


test_ic: ShallowWater.cpp initialCondition.cpp
	g++ -o target/initialCondition target/ShallowWater.o target/initialCondition.o -lblas
	./target/initialCondition

initialCondition.cpp:
	g++ -Wall -Wshadow -O2 -o target/initialCondition.o -c tests/initialCondition.cpp
	