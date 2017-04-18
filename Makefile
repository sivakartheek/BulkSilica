
main:main_silica.cpp ewald.cpp angular.cpp
	mpic++ main_silica.cpp ewald.cpp angular.cpp -o main 

clean:
	rm -f main
	rm -f output/graph/*
	rm -f output/data/*

