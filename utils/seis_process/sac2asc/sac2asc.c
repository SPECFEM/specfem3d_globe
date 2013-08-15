
#include <stdlib.h>
#include <stdio.h>
#include <sacio.h>

#define FALSE 0
#define TRUE 1

int
main(int argc, char *argv[]) {
    SACHEAD sachead;
    float *data;
    char* filename;
    int i, swap_bytes;
    float time;

/* DK DK August 2013: by default, do not swap bytes */
   swap_bytes = FALSE;

    if(argc < 2) {
        fprintf(stderr, "sac2asc: sacfile\n");
        exit(-1);
    }
    filename = argv[1];
    if((data = read_sac(filename, &sachead, swap_bytes) ) == 0) {
        fprintf(stderr, "Error reading sacfile: %s\n", filename);
        exit(-1);
    }
    time = sachead.b;
    for(i = 0; i < sachead.npts; i++) {
        fprintf(stdout, "%f %e\n", time, data[i]);
        time = time + sachead.delta;
    }
    free(data);

    return(0);
}

