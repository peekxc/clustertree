include $(MkInclude)

## Library name
LIB = libDT.a

# Dual tree source files components
DT_SRC_FILES = $(wildcard *.cpp ANNdt/*.cpp KNN/*.cpp MST/*.cpp)
OBJS = $(DT_SRC_FILES:.cpp=.o)

# Compile the static, callable library
library:
	$(MAKE) $(LIB)

## Library depends on the dual tree object files
$(LIB): $(OBJS)
	$(AR) -rucs $(LIB) $(OBJS)

## Rule to compile dual tree source files
%.o: %.cpp $(DT_SRC_FILES)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I$(SRC_DIR) -I$(DT_PATH) -c $< -o $@

# Remove object files, but keep the compiled programs and library archives
clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o
	@-rm -f $(OBJS)
