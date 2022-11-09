TARGET ?= GraphCounting.out
SRC_DIRS ?= ./

SRCS := $(shell find $(SRC_DIRS) -name *.cpp -or -name *.c -or -name *.s)
OBJS := $(addsuffix .o,$(basename $(SRCS)))
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CCOMPILE=mpic++
CPPFLAGS ?= $(INC_FLAGS) -fopenmp -g0 -O3 -march=native -std=c++14 -Wall -Wextra
CPPFLAGS += -g

LDLIBS = -ltbb -pthread

$(OBJS): 
	$(CCOMPILE) $(CPPFLAGS) $(LDFLAGS) -c -o $@ $(SRCS)

all: $(OBJS)
	$(CCOMPILE) $(CPPFLAGS) $(LDFLAGS) $(OBJS) -o $(TARGET) $(LDLIBS)

.PHONY: clean
clean:
	$(RM) $(TARGET) $(OBJS) $(DEPS)

-include $(DEPS)
