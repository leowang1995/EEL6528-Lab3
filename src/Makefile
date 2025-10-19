# University of Florida EEL6528
# Tan F. Wong
# Aug 31, 2023

CXX = g++
CPPFLAGS = -O3
LDFLAGS = -luhd -lboost_program_options
OBJDIR := ../build

common_hdrs = lab3.hpp 

all: lab2

lab2: $(common_hdrs) $(OBJDIR) $(OBJDIR)/lab3.o
	$(CXX) -o lab3 $(OBJDIR)/lab3.o $(LDFLAGS)

$(OBJDIR):
	mkdir -p $@

$(OBJDIR)/%.o: %.cpp $(common_hdrs)
	$(CXX) $(CPPFLAGS) -c $< -o $@

clean: 
	rm -rf lab3 $(OBJDIR)
