HTSLIB=$(HOME)/dxh/workspace/htslib-sw
GPTL=/home/export/online1/swmore/release/lib/
INCLUDES=-I$(HTSLIB) -I$(GPTL) -I./lwpf2


all: swmapper serial_index mpi_index debug_find_pos

serial_index:
	mpicc -g  -O2 -c serial_index.c -Wl,--wrap=malloc -L/usr/sw-mpp/mpi2/lib -L/usr/sz/mpi2/lib -L/usr/sz/lib -lmpich $(INCLUDES) -L$(GPTL) -lgptl $(LINK_SPC) -o  serial_index.o
	sw5cc -slave -g -c -O2 cpe_index.c -Wl,--wrap=malloc -L/usr/sw-mpp/mpi2/lib -L/usr/sz/mpi2/lib -L/usr/sz/lib -lmpich $(INCLUDES) -L$(GPTL) -lgptl $(LINK_SPC) -o cpe_index.o 
	mpicc -hybrid -g -O2 serial_index.o cpe_index.o -Wl,--wrap=malloc -L/home/export/online1/swmore/opensource/swlu/lib -lswlu_mpi -L/usr/sw-mpp/mpi2/lib -L/usr/sz/mpi2/lib -L/usr/sz/lib -lmpich $(INCLUDES) -L$(GPTL) -lgptl -lz $(LINK_SPC) -o serial_index.out

swmapper:
	mpicc -g -O2 -c mapper_mpe.c -L/usr/sw-mpp/mpi2/lib -L/usr/sz/mpi2/lib -L/usr/sz/lib -lmpich $(INCLUDES) -L$(GPTL) -lgptl $(LINK_SPC) -o mapper_mpe.o
	sw5cc -slave -g -O2 -msimd -c mapper_cpe.c -L/usr/sw-mpp/mpi2/lib -L/usr/sz/mpi2/lib -L/usr/sz/lib -lmpich $(INCLUDES) -L$(GPTL) -lgptl $(LINK_SPC) -o mapper_cpe.o
	mpicc -hybrid -g -O2 mapper_mpe.o mapper_cpe.o -Wl,--wrap=malloc -L/home/export/online1/swmore/opensource/swlu/lib -lswlu_mpi -L/usr/sw-mpp/mpi2/lib -L/usr/sz/mpi2/lib -L/usr/sz/lib -lmpich $(INCLUDES) -L$(GPTL) -lgptl -lz $(LINK_SPC) -o swmapper.out

mpi_index:
	mpicc -g  -O2 -c mpi_index.c -Wl,--wrap=malloc -L/usr/sw-mpp/mpi2/lib -L/usr/sz/mpi2/lib -L/usr/sz/lib -lmpich $(INCLUDES) -L$(GPTL) -lgptl $(LINK_SPC) -o  mpi_index.o
	mpicc -hybrid -g -O2 mpi_index.o -Wl,--wrap=malloc -L/home/export/online1/swmore/opensource/swlu/lib -lswlu_mpi -L/usr/sw-mpp/mpi2/lib -L/usr/sz/mpi2/lib -L/usr/sz/lib -lmpich $(INCLUDES) -L$(GPTL) -lgptl -lz $(LINK_SPC) -o mpi_index.out


debug_find_pos:
	mpicc -g -O2 find_pos.c -Wl,--wrap=malloc -L/home/export/online1/swmore/opensource/swlu/lib -lswlu_mpi -L/usr/sw-mpp/mpi2/lib -L/usr/sz/mpi2/lib -L/usr/sz/lib -lmpich $(INCLUDES) -L$(GPTL) -lgptl $(LINK_SPC) -o debug_find_pos.out


clean:
	rm -f *.o
	rm -f swmapper.out serial_index.out debug_find_pos.out mpi_index.out
