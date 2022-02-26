
TARGET = mc

SRC  = $(shell find ./src -type f -name *.c)
HEADS = $(shell find ./include -type f -name *.h)
OBJS = $(SRC:.c=.o)
DEPS = Makefile.depend

INCLUDE = -I./include
CFLAGS = -O2 -Wall $(INCLUDE)
LDFLAGS = -lm

all: $(TARGET)

$(TARGET): $(OBJS) $(HEADS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS)

run: all
	@./$(TARGET)

.PHONY: depend clean
depend:
	$(CFLAGS) $(INCLUDE) -MM $(SRCS) > $(DEPS)
	@sed -i -E "s/^(.+?).o: ([^ ]+?)\1/\2\1.o: \2\1/g" $(DEPS)

clean:
	$(RM) $(OBJS) $(TARGET)

-include $(DEPS)