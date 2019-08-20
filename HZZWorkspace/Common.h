// Keep common defined types and head files
#ifndef _HZZWS_COMMON_H
#define _HZZWS_COMMON_H

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <TString.h>

// new type defines..
typedef std::vector<std::string> strvec;
typedef std::vector<TString> tstrvec;
typedef std::map<std::string, std::string> strmap;
typedef std::map<std::string, double> dblmap;
typedef std::map<std::string, strmap> strdic;
typedef std::map<std::string, dblmap> dbldic;
static strmap DEFAULT_STRMAP;

// define help function for ERRORs, Warnings and Infos
// and avoid variadic preprocessor macros if we can
template <typename... inArgs> void log_write (const std::string & kind, inArgs... theArgs){
    std::cout << "["<<kind<<"]"<< "("<<__FILE__<<":"<<__LINE__<<") "<<Form(theArgs...)<<std::endl;
}
template <typename... inArgs> void log_err (inArgs... theArgs) {
    log_write("ERROR", theArgs...);
}
template <typename... inArgs> void log_warn (inArgs... theArgs) {
    log_write("Warning", theArgs...);
}
template <typename... inArgs> void log_info (inArgs... theArgs) {
    log_write("INFO", theArgs...);
}

#endif
