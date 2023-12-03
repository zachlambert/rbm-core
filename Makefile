.PHONY: build
build:
	mkdir -p build
	cmake -E chdir build cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DUSE_SVIZ=0 ..
	cmake --build build

.PHONY: build_sviz
build_sviz:
	mkdir -p build
	cmake -E chdir build cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DUSE_SVIZ=1 ..
	cmake --build build

.PHONY: test
test:
	cd build && ctest --output-on-failure

.PHONY: install
install:
	cmake --build build --target install
