//////////////////////////////////////////////////////////////////
// little helper class to manage selection of electron candidates
#include <vector>

class SelectElectron {
  std::vector<int> allowed;
 public:
  void add(int);
  bool doesElePass(int);
};
