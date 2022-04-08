OBJS	= driver.o flux.o helper.o mesh.o mesh_refinement.o pressure.o stdatmos.o timemarching.o
SOURCE	= driver.c flux.c helper.c mesh.c mesh_refinement.c pressure.c stdatmos.c timemarching.c
HEADER	= header.h
OUT	= CompEuler
LFLAGS	 = -lpthread -lm

all: CompEuler

CompEuler: $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

%.o: %.c $(HEADER)
	$(CC) -c -o $@ $< $(LFLAGS)

clean:
	rm -f $(OBJS) $(OUT)

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