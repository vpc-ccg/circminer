#include "common.h"

// verbose-aware fprintf
void vafprintf(int verbosity, FILE *stream, const char *format, ...) {
	if (verbosity > verboseMode)	return;

	va_list args;
	va_start (args, format);
	vfprintf (stream, format, args);
	va_end (args);
}
