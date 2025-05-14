CC = gcc
CFLAGS = -Wall -Wextra -O2 -DHAVE_INLINE -std=gnu2y
CFLAGS += $(EXTRA)
LDFLAGS = -lgsl -lgslcblas -lm
SRCS = main.c wavefunction.c constants.c lse.c
HEADERS = wavefunction.h constants.h lse.h
OBJS = $(SRCS:.c=.o)
TARGET = GaussianVariation

# Shared library settings
WFSRCS = wavefunction.c constants.c
WFHEADERS = $(WFSRCS:.c=.h)
WFOBJS = $(WFSRCS:.c=.pic.o)
WFTARGET = libwavefunction.so

# LSE shared library - add constants.c as a dependency
LSESRCS = lse.c constants.c wavefunction.c
LSEHEADERS = $(LSESRCS:.c=.h)
LSEOBJS = $(LSESRCS:.c=.pic.o)
LSETARGET = liblse.so

all: $(TARGET) $(WFTARGET) $(LSETARGET)

# Original executable target
$(TARGET): $(OBJS) 
	$(CC) $(CFLAGS) $(EXTRA) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# New shared library target
wf: $(WFTARGET)

$(WFTARGET): $(WFOBJS) 
	$(CC) -shared -o $@ $^ $(LDFLAGS)

%.pic.o: %.c
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

lseso: $(LSETARGET)

$(LSETARGET): $(LSEOBJS) 
	$(CC) -shared -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(OBJS) $(WFOBJS) $(LSEOBJS) $(TARGET) $(WFTARGET) $(LSETARGET)
