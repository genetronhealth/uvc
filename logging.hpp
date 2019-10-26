// Log, version 0.1: a simple logging class
#ifndef logging_hpp_INCLUDED
#define logging_hpp_INCLUDED

#include <sstream>
#include <time.h>

#define LOG(level) \
if (level > Log::ReportingLevel()) ; \
else Log().Get(level)

enum              TLogLevel         { logCRITICAL ,  logERROR ,  logWARNING ,  logINFO ,  logINFO2,   logDEBUG ,  logDEBUG1 ,  logDEBUG2 ,  logDEBUG3 ,  logDEBUG4 };

static TLogLevel globalMessageLevel = logINFO2;

class Log {
public:
    Log() {};
    virtual ~Log();
    std::ostringstream& 
    Get(TLogLevel level = logINFO2);
public:
    static TLogLevel& 
    ReportingLevel() { return globalMessageLevel; };
protected:
    std::ostringstream os;
private:
    Log(const Log&);
    Log& operator =(const Log&);
private:
    TLogLevel messageLevel;
};

char *
NowTime(char *buffer);

#endif
