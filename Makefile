.PHONY: all build test clean run bench native web-serve

all: build

build:
	./gradlew build -x test

test:
	./gradlew test

clean:
	./gradlew clean
	$(MAKE) -C c clean 2>/dev/null || true

run:
	./gradlew run --console=plain -q

bench:
	./gradlew run --args="--bench" -q

native:
	$(MAKE) -C c all

web-serve:
	@echo "Starting local HTTP server for Calabi-Yau visualizer on http://localhost:8080"
	python3 -m http.server 8080 --directory web
