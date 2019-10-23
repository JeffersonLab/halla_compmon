#include <string>
#include <vector>

void exampleMacro3(){
  double x = 10;
  double y = 1.3;
  double z = TMath::Power(x,y);
  std::cout << "z = " << z << "\n";
  std::cout << "z = " << z << "\n";
  std::vector<std::string> s; s.push_back("Some dumb bullshit");
  s.push_back("Some dumber bullshit.");
  s.push_back("Only the dumbest bullshit.");
  for(int i = 0; i < s.size(); i++){
    std::cout << s[i] << std::endl;
  }

  TH1F *h1 = new TH1F("h1", "Histo", 200, -10, 10);
  h1->FillRandom("gaus", 10000);
  h1->Draw();
}
