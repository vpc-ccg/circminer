#include "common.h"

FILE* open_file(char* filename, char* mode) {
	FILE* fp;
	fp = fopen(filename, mode);
	if (fp == NULL) {
		fprintf(stderr, "Error: Could not open file %s\n", filename);
		exit(1);
	}
	return fp;
}

void close_file(FILE* fp) {
	if (fp != NULL)
		fclose(fp);
}

// verbose-aware fprintf
void vafprintf(int verbosity, FILE *stream, const char *format, ...) {
	if (verbosity > verboseMode)	return;

	va_list args;
	va_start (args, format);
	vfprintf (stream, format, args);
	va_end (args);
}
