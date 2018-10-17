#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "OrangeConfig.h"


void help() {
    printf(
        "\n"
        "The program accepts two files as floating arguments and outputts number of sequences, average length, minimal and maximal length.\n"
        "Supported formats are .fasta, .fa, .fastq, .fq, .fasta.gz, .fa.gz, .fastq.gz, .fq.gz\n"
        "\n"
        "usage: orange_mapper <file_name> <file_name>\n"
        "\n"
        "options: \n"
        "      -h, --help  -  prints help menu (currently displaying)\n"
        "      -v, --version  -  prints help menu (currently displaying)\n");
}

void version() {
    printf("Version: %d\n",
        orange_mapper_VERSION_MAJOR
    );
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
            version();
            exit(0);
        }
        default: {
            printf("Option not recognised.\n");
            exit(1);
        }
        }
    }
}