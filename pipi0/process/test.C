#include "Settings.C"

int test()
{
  Settings s = Settings();
  s.readCard("testcard.txt");
  s.loadSettings();
  
  cout << s.ymax() << endl;
  return 0;
}
