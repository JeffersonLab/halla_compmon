TTree *quartetwise = 0;

void chainHelper()
{
  TFile *rootfile = gROOT->GetFile();
  quartetwise = (TTree*)rootfile->Get("quartetwise");
  quartetwise->SetAlias("HP0","PosHelAcc0/PosHelNSamples0");
  quartetwise->SetAlias("HN0","NegHelAcc0/NegHelNSamples0");
  quartetwise->SetAlias("HS0","HP0+HN0");
  quartetwise->SetAlias("HD0","HP0-HN0");
  quartetwise->SetAlias("HA0","HD0/HS0");
  for(Int_t a = 1; a < 6; a++) {
    quartetwise->SetAlias(Form("HP%d",a),Form("PosHelAcc%d",a));
    quartetwise->SetAlias(Form("HN%d",a),Form("NegHelAcc%d",a));
    quartetwise->SetAlias(Form("HS%d",a),Form("HP%d+HN%d",a,a));
    quartetwise->SetAlias(Form("HD%d",a),Form("HP%d-HN%d",a,a));
    quartetwise->SetAlias(Form("HA%d",a),Form("HD%d-HS%d",a,a));
  }

}
