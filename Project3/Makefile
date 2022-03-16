OBJS	= boundary.o debugger.o driver.o helper.o multigrid.o timeintegration.o updater.o
SOURCE	= boundary.c debugger.c driver.c helper.c multigrid.c timeintegration.c updater.c
HEADER	= headers.h
OUT	= Project3
LFLAGS	 = -lm -g -O3 -march=native

all: Project3

Project3: $(OBJS)
	$(CC) $(LFLAGS) -o $@ $^

%.o: %.c $(HEADER)
	$(CC) $(LFLAGS) -c -o $@ $< 

clean:
	rm -f $(OBJS) $(OUT) uVelocity.csv vVelocity.csv residuals.csv pressure.csv pressure2.csv

run: $(OUT)
	./$(OUT)

debug: $(OUT)
	valgrind $(OUT)

valgrind: $(OUT)
	valgrind $(OUT)

valgrind_leakcheck: $(OUT)
	valgrind --leak-check=full $(OUT)

valgrind_extreme: $(OUT)
	valgrind --leak-check=full --show-leak-kinds=all --leak-resolution=high --track-origins=yes --vgdb=yes $(OUT)
