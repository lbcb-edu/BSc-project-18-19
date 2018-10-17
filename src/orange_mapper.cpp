#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>


void help() {
    printf(
        "This is help menu");
}

int main(int argc, char** argv) {

    static struct option options[] = {
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {NULL, no_argument, NULL, 0}
    };

    char optchr;

    while ((optchr = getopt_long(argc, argv, "hv", options, NULL)) != -1) {
        switch (optchr) {
        case 'h': {
            help();
            exit(0);
        }
        case 'v': {
            printf("Version: 1\n");
            exit(0);
        }
        default: {
            printf("Option not recognised.\n", argv[0]);
            exit(1);
        }
        }
    }
}