#include "ldPreloder.h"

int inited = 0;

void init(void) {
  int i = 0;

  if (inited) {
    return;
  }
  printf("init\n");
  while (bindings[i].name) {
    *bindings[i].real_address = dlsym(RTLD_NEXT, bindings[i].name);
    if (*bindings[i].real_address == NULL) {
      fprintf(stderr, "Error in `dlsym` of %s: %s\n",
              bindings[i].name, dlerror());
    }
    i++;
  }

  inited = 1;
}
