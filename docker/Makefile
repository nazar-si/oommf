build: Dockerfile Makefile
	docker build --no-cache -t field/oommf:latest .

run:
	docker run --rm -ti -v `pwd`:/io field/oommf:latest bash

test:
	docker run --rm -v `pwd`:/io field/oommf:latest oommf boxsi +fg /usr/local/oommf/oommf/app/oxs/examples/stdprob3.mif -exitondone 1


