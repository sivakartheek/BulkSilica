
main:main_silica.cpp ewald.cpp
	mpic++ main_silica.cpp ewald.cpp -o main 

clean:
	rm -f main
	rm -f output/data/*

