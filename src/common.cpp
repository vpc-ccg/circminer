#include <sys/time.h>
#include <sys/resource.h>

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

double get_cpu_time() {
	struct rusage t;
	getrusage(RUSAGE_SELF, &t);
	return t.ru_utime.tv_sec + t.ru_utime.tv_usec / 1000000.0 + t.ru_stime.tv_sec + t.ru_stime.tv_usec / 1000000.0;
}

double get_real_time() {
	struct timeval t;
	struct timezone tz;
	gettimeofday(&t, &tz);
	return t.tv_sec + t.tv_usec / 1000000.0;
}
