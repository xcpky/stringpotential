CC = gcc
CFLAGS = -Wall -Wextra -O2
LDFLAGS = -lgsl -lgslcblas -lm
SRCS = main.c wavefunction.c
OBJS = $(SRCS:.c=.o)
TARGET = GaussianVariation

# Shared library settings
LIBSRCS = wavefunction.c
LIBOBJS = $(LIBSRCS:.c=.pic.o)
LIBTARGET = libwavefunction.so

all: $(TARGET)

# Original executable target
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# New shared library target
shared: $(LIBTARGET)

$(LIBTARGET): $(LIBOBJS)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

%.pic.o: %.c
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

clean:
	rm -f $(OBJS) $(LIBOBJS) $(TARGET) $(LIBTARGET)
