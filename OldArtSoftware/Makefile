 
IMAGE = .
 
# The name of the C++ compiler
CC = g++
 
# Flags to send the compiler (change the last two with care...)
CFLAGS = -g
 
# Object, Runnable and Include directories
 
ODIR = $(IMAGE)/obj
IDIR = $(IMAGE)/include
RDIR = $(IMAGE)/bin
SDIR = $(IMAGE)/src

CART = /home/ensab/Robots/Cartesian
 
 
 
FL = -g -Wno-deprecated -I$(IDIR) -O2 -I/usr/include/freetype2 -D_REENTRANT -D_FILE_OFFSET_BITS=64 -I/usr/X11R6/include -I/usr/X11R6/include/X11/magick -I$(CART)/include 

GB = -L/usr/X11R6/lib -lMagick++ -lMagick -ltiff -ljpeg -lpng -ldpstk -ldps -lXext -lXt -lSM -lICE -lX11 -lbz2 -lz -lpthread -lm
 
art_image:  $(ODIR)/art_image.o
			 $(CC) -o $(RDIR)/art_image $(ODIR)/art_image.o $(GB)
 
# Make the objects
 
$(ODIR)/art_image.o:  $(SDIR)/art_image.cxx
			 $(CC) -c $(FL) -o $(ODIR)/art_image.o $(SDIR)/art_image.cxx	

segment:  $(ODIR)/segment.o
			 $(CC) -o $(RDIR)/segment $(ODIR)/segment.o $(GB)
 
# Make the objects
 
$(ODIR)/segment.o:  $(SDIR)/segment.cxx
			 $(CC) -c $(FL) -o $(ODIR)/segment.o $(SDIR)/segment.cxx

draw:  $(ODIR)/draw.o
			 $(CC) -o $(RDIR)/draw $(ODIR)/draw.o $(CART)/obj/cartesian.o $(GB)
 
# Make the objects
 
$(ODIR)/draw.o:  $(SDIR)/draw.cxx $(IDIR)/draw.h
			 $(CC) -c $(FL) -o $(ODIR)/draw.o $(SDIR)/draw.cxx				    
