#include "Hzzws/Helper.h"

int main(){

  std::vector<double> vec;

  //try anything here
  Helper::readAcceptancePoly(vec, "VBF", "ggF_2mu2e");

  for (int i(0);i<vec.size();++i)
    std::cout<<vec[i]<<std::endl;
  vec.clear();
  
  //try anything here
  Helper::readAcceptancePoly(vec, "ggH", "ggF_4e");

  for (int i(0);i<vec.size();++i)
    std::cout<<vec[i]<<std::endl;
  vec.clear();
  
}

  
