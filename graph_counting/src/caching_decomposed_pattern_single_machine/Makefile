TARGET ?= GraphCounting.out
SRC_DIRS ?= ./

SRCS := $(shell find $(SRC_DIRS) -name *.cpp -or -name *.c -or -name *.s)
OBJS := $(addsuffix .o,$(basename $(SRCS)))
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -g0 -O3 -march=native -std=c++14 -Wall -Wextra
CPPFLAGS += -g

LDLIBS = -ltbb -pthread

$(TARGET): $(OBJS)
	g++ $(CPPFLAGS) $(LDFLAGS) $(OBJS) -o $@ $(LDLIBS)

.PHONY: clean
clean:
	$(RM) $(TARGET) $(OBJS) $(DEPS)

-include $(DEPS)
