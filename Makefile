a.out:
	clang++ -O3 -pipe -lpthread -march=native -Wall -flto -DNDEBUG main.cpp -o a.out
clean:
	rm *.o a.out
