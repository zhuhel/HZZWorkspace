// Keep common defined types and head files
#ifndef _HZZWS_COMMON_H
#define _HZZWS_COMMON_H

#include <map>
#include <string>
#include <vector>
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
#ifndef log_err
#define log_err(M, ...) fprintf(stdout, "[ERROR] (%s:%d) " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#define log_warn(M, ...) fprintf(stdout, "[Warning] (%s:%d) " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#define log_info(M, ...) fprintf(stdout, "[INFO] (%s:%d) " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#endif

#endif
