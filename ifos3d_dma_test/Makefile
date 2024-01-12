MCC := mpicc


LFLAGS=-lm -lmpi -lcseife
CFLAGS=-Wall -O2
# CFLAGS=-Wall -O4
SFLAGS=-L./../libcseife 
IFLAGS=-I./../libcseife

BIN_NAME := ifos3d

BIN_PATH := ./bin
BUILD_PATH := ./build
SRC_PATH := ./src


SRC_C := $(wildcard $(SRC_PATH)/*.c)
DIR_C := $(notdir $(SRC_C))
OBJS := $(patsubst %.c,$(BUILD_PATH)/%.o,$(DIR_C))




all: $(OBJS)
	$(MCC) $(OBJS) $(SFLAGS) $(LFLAGS) -o $(BIN_PATH)/$(BIN_NAME)
$(BUILD_PATH)/%.o : $(SRC_PATH)/%.c
	$(MCC) $(CFLAGS) $(IFLAGS) -c $< -o $@


.PHONY: clean
clean:
	-rm -f $(BIN_PATH)/$(BIN_NAME)
	-rm -f $(BUILD_PATH)/*.o



