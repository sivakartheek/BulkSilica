
main:main_silica.cpp ewald.cpp
	/cluster/mpi/openmpi/1.6.5-gcc4.8.2/bin/mpic++ main_silica.cpp ewald.cpp -o main 

clean:
	rm -f main
	rm -f output/data/*

