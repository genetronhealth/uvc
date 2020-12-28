#ifndef version_INCLUDED
#define version_INCLUDED
#include <string.h>

#define VERSION123 "0.1.10"

#ifndef COMMIT_VERSION
#define COMMIT_VERSION "NotVersionControlled"
#endif

#ifndef COMMIT_DIFF_SH
#define COMMIT_DIFF_SH "NoVersion"
#endif

#define VERSION_CLEAN VERSION123 "." COMMIT_VERSION
#define VERSION_DIRTY VERSION123 "." COMMIT_VERSION "-dirty"
#define VERSION ((strlen(COMMIT_DIFF_SH) > 0) ? (VERSION_DIRTY) : (VERSION_CLEAN))

#define VERSION_DETAIL_CLEAN (VERSION_CLEAN " (" COMMIT_DIFF_SH ")")
#define VERSION_DETAIL_DIRTY (VERSION_DIRTY " (" COMMIT_DIFF_SH ")")
#define VERSION_DETAIL ((strlen(COMMIT_DIFF_SH) > 0) ? (VERSION_DETAIL_DIRTY) : (VERSION_DETAIL_CLEAN))

extern const char *GIT_DIFF_FULL;

#endif
