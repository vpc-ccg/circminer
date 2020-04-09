#ifndef __COMMANDLINEPARSER_H__
#define __COMMANDLINEPARSER_H__

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>

#include <string>

#include "common.h"

int parse_command(int argc, char *argv[]);
void printHELP(void);
void printHELPShort(void);

#endif    //__COMMANDLINEPARSER_H__
