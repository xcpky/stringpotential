CC = gcc
CFLAGS = -Wall -Wextra -O2
LDFLAGS = -lgsl -lgslcblas -lm
SRCS = main.c wavefunction.c constants.c lse.c
OBJS = $(SRCS:.c=.o)
TARGET = GaussianVariation

# Shared library settings
WFSRCS = wavefunction.c constants.c
WFOBJS = $(WFSRCS:.c=.pic.o)
WFTARGET = libwavefunction.so

# LSE shared library - add constants.c as a dependency
LSESRCS = lse.c constants.c wavefunction.c
LSEOBJS = $(LSESRCS:.c=.pic.o)
LSETARGET = liblse.so

all: $(TARGET)

# Original executable target
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

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
