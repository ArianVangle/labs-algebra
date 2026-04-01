CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall -I./include

TARGET = lab1

SRC_DIR = src
INC_DIR = include

SRCS = $(SRC_DIR)/main.cpp \
       $(SRC_DIR)/utils.cpp \
       $(SRC_DIR)/matrix_ops.cpp \
       $(SRC_DIR)/generators.cpp \
       $(SRC_DIR)/gauss.cpp \
       $(SRC_DIR)/lu.cpp \
       $(SRC_DIR)/experiments.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS)

clean:
	rm -f $(TARGET) exp4_*.csv

.PHONY: all run clean help