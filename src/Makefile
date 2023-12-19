VPATH = ./Equation:./Lib:./Post_Precessing:./Pre_Processing:./SPH
# 正则表达式表示目录下所有.c文件
SRC = $(notdir $(wildcard ./*/*.cpp))
# 头文件路径
INCLUDE = -I./Headers
# 生成文SRC
BUILD_DIR = ./build
OBJ = $(addprefix $(BUILD_DIR)/,$(patsubst %.cpp,%.o,$(SRC))) 


TARGET = run

CC = g++

.PHONY:all clean

all: $(TARGET)

# 第一行依赖关系：冒号后面为依赖的文件，相当于Hello: main.o a.o b.o
# 第二行规则：$@表示目标文件，$^表示所有依赖文件，$<表示第一个依赖文件
$(TARGET): $(OBJ)
	$(CC) $^ -o $@ -fopenmp

# 上一句目标文件依赖一大堆.o文件，这句表示所有.o都由相应名字的.cpp文件自动生成
$(BUILD_DIR)/%.o:%.cpp
	if [ ! -d $(BUILD_DIR) ]; then mkdir -p $(BUILD_DIR); fi;\
	$(CC) -c -fopenmp $(INCLUDE) -o $@ $^ 


# make clean删除所有.o和目标文件
clean:
	rm -f $(OBJ) $(TARGET)

debug:
	echo $(SRC) $(OBJ)
