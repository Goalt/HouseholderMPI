procs = 16 32 64 128 256 512 1024
# procs = 1 2 4 8

all:
	mpicxx main.cpp functions.cpp -o res -O3
clean:
	rm -rf res
test1:
	mpisubmit.bg -n 1 -w 00:30:00 res 1024
	mpisubmit.bg -n 2 -w 00:20:00 res 1024
	mpisubmit.bg -n 4 -w 00:15:00 res 1024
	mpisubmit.bg -n 8 -w 00:15:00 res 1024

	for v in $(procs) ; \
	do \
		mpisubmit.bg -n $$v -w 00:02:00 res 1024; \
	done
test2:
	mpisubmit.bg -n 1 -w 00:30:00 res 2048
	mpisubmit.bg -n 2 -w 00:20:00 res 2048
	mpisubmit.bg -n 4 -w 00:15:00 res 2048
	mpisubmit.bg -n 8 -w 00:15:00 res 2048

	for v in $(procs) ; \
	do \
		mpisubmit.bg -n $$v -w 00:02:00 res 2048; \
	done
test3:
	mpisubmit.bg -n 1 -w 00:30:00 res 3072
	mpisubmit.bg -n 2 -w 00:20:00 res 3072
	mpisubmit.bg -n 4 -w 00:15:00 res 3072
	mpisubmit.bg -n 8 -w 00:15:00 res 3072

	for v in $(procs) ; \
	do \
		mpisubmit.bg -n $$v -w 00:02:00 res 3072; \
	done