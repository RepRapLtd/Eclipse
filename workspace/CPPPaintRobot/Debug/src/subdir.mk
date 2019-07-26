################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Painter.cpp \
../src/SingleLinkage.cpp \
../src/SingleLinkageOptimized.cpp \
../src/Slink.cpp \
../src/Utils.cpp 

OBJS += \
./src/Painter.o \
./src/SingleLinkage.o \
./src/SingleLinkageOptimized.o \
./src/Slink.o \
./src/Utils.o 

CPP_DEPS += \
./src/Painter.d \
./src/SingleLinkage.d \
./src/SingleLinkageOptimized.d \
./src/Slink.d \
./src/Utils.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/usr/include/opencv2 -I/usr/include/opencv -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


