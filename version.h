#ifndef version_INCLUDED
#define version_INCLUDED
#include <string.h>
#include <string>

#define VERSION123 "0.0.3"

#ifndef COMMIT_VERSION
#define COMMIT_VERSION "NotVersionControlled"
#endif

#ifndef COMMIT_DIFF_SH
#define COMMIT_DIFF_SH "NoVersion"
#endif

const std::string VERSION = std::string() + VERSION123 + "." + COMMIT_VERSION + (strlen(COMMIT_DIFF_SH) > 0 ? "-dirty" : "");
const std::string VERSION_DETAIL = VERSION + " (" + COMMIT_DIFF_SH + ")";

#endif
