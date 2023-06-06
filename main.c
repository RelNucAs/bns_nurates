#include <stdio.h>
#include "src/constants.h"
#include "src/bns_nurates.h"
#include "src/integration/integration.h"
#include <unistd.h>
#include <limits.h>
int main() {
  char cwd[PATH_MAX];
  if (getcwd(cwd, sizeof(cwd)) != NULL) {
       printf("Current working dir: %s\n", cwd);
   } else {
       perror("getcwd() error");
       return 1;
   }
  printf("Hello, World!\n");
  return 0;
}
