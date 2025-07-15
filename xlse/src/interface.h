#ifndef INTERFACE_H
#define INTERFACE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
      void *ose00; // Using void* to avoid complex.h in C++ code
      void *ose01;
      void *ose10;
      void *ose11;
} onshellElements;

// Declare the function with void* to avoid complex.h dependency in C++
onshellElements onshellT(void *E, size_t len, size_t pNgauss);

void ose_free(onshellElements ose);

#ifdef __cplusplus
}
#endif

#endif // !INTERFACE_H
