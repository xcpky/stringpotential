#include "utils.h"
#include <stdio.h>
char *formatC(double complex x) {
    static char buf[64];
    sprintf(buf, "%.9f%+.9fim", creal(x), cimag(x));
    return buf;
}

void writec(const char *filename, double complex *data, uint64_t n) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        perror("Failed to open file for writing\n");
        return;
    }
    fwrite(data, sizeof(double complex), n, fp);
    fclose(fp);
}
void readc(const char *filename, double complex *data, uint64_t n) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("Failed to open file for reading\n");
        return;
    }
    fread(data, sizeof(double complex), n, fp);
    fclose(fp);
}

void writef(const char *filename, double *data, uint64_t n) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        perror("Failed to open file for writing\n");
        return;
    }
    fwrite(data, sizeof(double), n, fp);
    fclose(fp);
}
void readf(const char *filename, double *data, uint64_t n) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("Failed to open file for reading\n");
        return;
    }
    fread(data, sizeof(double), n, fp);
    fclose(fp);
}
