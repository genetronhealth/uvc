#include "logging.hpp"

#define LOG(level) \
if (level > Log::ReportingLevel()) ; \
else Log().Get(level)

// Defined in header: enum TLogLevel{ logCRITICAL ,  logERROR ,  logWARNING ,  logINFO ,  logINFO2,   logDEBUG ,  logDEBUG1 ,  logDEBUG2 ,  logDEBUG3 ,  logDEBUG4 };
const char *TLogLevelToString[10] = {"logCRITICAL", "logERROR", "logWARNING", "logINFO", "logINFO2", "logDEBUG", "logDEBUG1", "logDEBUG2", "logDEBUG3", "logDEBUG4"};

static TLogLevel globalMessageLevel = logINFO2;

char *
nowtime(char *buffer) {
    time_t rawtime;
    time(&rawtime);
    struct tm t;
    localtime_r(&rawtime, &t);
    strftime(buffer, 128, "%F %T %z", &t);
    return buffer;
}

TLogLevel& Log::ReportingLevel() { return globalMessageLevel; };

std::ostringstream& 
Log::Get(TLogLevel level) {
    char buffer[128];
    os << "- " << nowtime(buffer);
    os << " " << TLogLevelToString[level] << ": ";
    messageLevel = level;
    return os;
}

Log::~Log() {
    if (messageLevel <= Log::ReportingLevel()) {
        os << std::endl;
        fprintf(stderr, "%s", os.str().c_str());
        fflush(stderr);
    }
}

