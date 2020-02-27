#Compiler and Linker
ifdef GCC_8
CC          :=g++-8
else
CC			:=g++
endif

NVCC        :=/usr/local/cuda-9.0/bin/nvcc

#The Target Binary Program
TARGET      := namics

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := src
INCDIR      := inc
BUILDDIR    := obj
TARGETDIR   := bin
RESDIR      := res
SRCEXT      := cpp
CUDAEXT		:= cu
DEPEXT      := d
OBJEXT      := o

#Flags, Libraries and Includes
CFLAGS      := -Wall -Ofast -std=c++14 -march=native
LIB         := -lm -lpthread
INC         := -I/usr/local/cuda-9.0/include -I/usr/local/include -I/usr/include
#INCDEP      := -I$(INCDIR)
ifdef CLENG_EXPERIMENTAL
CFLAGS      := -Wall -Ofast -std=c++14 -DCLENG_EXPERIMENTAL
LIB         += -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5
INC         += -I/usr/include/hdf5/serial
endif
ifdef CUDA
	LIB        += -L/usr/local/cuda/lib64 -lcuda -lcudart -lcurand
	CFLAGS     += -DCUDA
	NVCCFLAGS  := -g -ccbin gcc-5 -arch=sm_61 -std=c++14 -DCUDA
	ifdef PAR_MESODYN
		CFLAGS += -DPAR_MESODYN
		NVCCFLAGS += --expt-relaxed-constexpr --expt-extended-lambda -DPAR_MESODYN
	endif
endif

# %.o: %.cu $(NVCC) $(NVCCFLAGS) -c $< -o $@
#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------

SOURCES     := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
CUDASOURCES := $(shell find $(SRCDIR) -type f -name *.$(CUDAEXT))
OBJECTS     := $(filter $(BUILDDIR)/%, $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT))) \
	$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(CUDASOURCES:.$(CUDAEXT)=.$(OBJEXT))))

#Defauilt Make
all: resources $(TARGET)

#Remake
remake: cleaner all

#Copy Resources from Resources Directory to Target Directory
resources: directories
#	@cp $(RESDIR)/* $(TARGETDIR)/

#Make the Directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)

#Clean only Objects
clean:
	@$(RM) -rf $(BUILDDIR)

#Full Clean, Objects and Binaries
cleaner: clean
	@$(RM) -rf $(TARGETDIR)

#Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)

#Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
	@$(CC) $(CFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp

ifdef CUDA
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(CUDAEXT)
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) $(INC) -c -o $@ $<
endif


#Non-File Targets
.PHONY: all remake clean cleaner resources
