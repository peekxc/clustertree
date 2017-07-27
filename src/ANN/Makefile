include $(MkInclude)

## Library name
# LIB = libANN.a

# Anything with a cpp extension is a source file
ANN_SRC_FILES = $(wildcard *.cpp **/*.cpp)
OBJS = $(ANN_SRC_FILES:.cpp=.o)

# Compile the static, callable library
library:
	$(MAKE) $(LIB)

#$(LIB): $(OBJS)
#	$(AR) -rucs $(LIB) $(OBJS)

%.o: %.cpp $(ANN_SRC_FILES)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I$(SRC_DIR) -I$(ANN_PATH) -c $< -o $@

# Remove object files, but keep the compiled programs and library archives
clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o
	@-rm -f $(OBJS)
