#ifndef __HZZWS_COEFFICIENTFACTORY_H__
#define __HZZWS_COEFFICIENTFACTORY_H__


#include "HZZWorkspace/Coefficient.h"
#include "HZZWorkspace/Helper.h"
#include <map>


namespace CoefficientFactory{

  Coefficient* CreateCoefficient(const std::string& coefline);

}

#endif
