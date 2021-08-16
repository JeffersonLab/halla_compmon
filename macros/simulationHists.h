#include "../online/utils.h"

using namespace std;

Float_t findComptonEdgeSim(TH1F *hSum, Int_t nBins){
  printf("  Finding Compton Edge in simulated histogram...\n");
  Float_t max = -1e16;
  Int_t maxBin = 0;
  for(Int_t bin = 1; bin <= nBins; bin++){
    Float_t center = hSum->GetBinCenter(bin);
    if(center < 50) continue;
    Bool_t replaceMax = (hSum->GetBinContent(bin) > max);
    max = replaceMax ? hSum->GetBinContent(bin) : max;
    maxBin = replaceMax ? bin : maxBin;
  }

  printf("  Maximum and bin found! Maximum: %.2f, Max Bin: %i Max SRAU: %.2f\n", max, maxBin, hSum->GetBinCenter(maxBin));

  Float_t ce;
  for(Int_t bin = maxBin; bin <= nBins; bin++){
    if(hSum->GetBinContent(bin) < 0.5*max){
      ce = hSum->GetBinCenter(bin);
      break;
    }
  }

  printf("  Found compton edge: %.2f\n", ce);

  return ce;
}

Float_t getNormFactorPile(TH1F *hPile, Int_t nBins, Float_t match){
  Int_t matchBinNum = 1;
  while(hPile->GetBinCenter(matchBinNum) < 1.5 && matchBinNum <= nBins){
    matchBinNum++;
  }
  return match*1.0/hPile->GetBinContent(matchBinNum);
}

Float_t getNormFactorSim(TH1F *hSpec, Int_t nBins){
  Float_t max = -1e16;
  Int_t maxBin = 0;
  for(Int_t bin = 1; bin <= nBins; bin++){
    Float_t center = hSpec->GetBinCenter(bin);
    if(center < 0.5) continue;
    Bool_t replaceMax = (hSpec->GetBinContent(bin) > max);
    max = replaceMax ? hSpec->GetBinContent(bin) : max;
    maxBin = replaceMax ? bin : maxBin;
  }

  return max;
}

void fillSimHist_6mm_raw(TH1F *h){
   h->SetBinContent(1,6304);
   h->SetBinContent(2,10891);
   h->SetBinContent(3,17129);
   h->SetBinContent(4,19810);
   h->SetBinContent(5,21139);
   h->SetBinContent(6,22667);
   h->SetBinContent(7,23298);
   h->SetBinContent(8,24292);
   h->SetBinContent(9,25254);
   h->SetBinContent(10,25813);
   h->SetBinContent(11,26697);
   h->SetBinContent(12,27511);
   h->SetBinContent(13,28499);
   h->SetBinContent(14,29259);
   h->SetBinContent(15,29928);
   h->SetBinContent(16,30931);
   h->SetBinContent(17,31611);
   h->SetBinContent(18,32598);
   h->SetBinContent(19,33844);
   h->SetBinContent(20,34302);
   h->SetBinContent(21,35292);
   h->SetBinContent(22,35404);
   h->SetBinContent(23,34813);
   h->SetBinContent(24,34718);
   h->SetBinContent(25,34383);
   h->SetBinContent(26,34158);
   h->SetBinContent(27,33690);
   h->SetBinContent(28,33313);
   h->SetBinContent(29,33290);
   h->SetBinContent(30,32837);
   h->SetBinContent(31,32649);
   h->SetBinContent(32,32582);
   h->SetBinContent(33,31750);
   h->SetBinContent(34,31857);
   h->SetBinContent(35,31445);
   h->SetBinContent(36,31270);
   h->SetBinContent(37,30977);
   h->SetBinContent(38,30303);
   h->SetBinContent(39,30195);
   h->SetBinContent(40,29526);
   h->SetBinContent(41,29721);
   h->SetBinContent(42,29424);
   h->SetBinContent(43,28894);
   h->SetBinContent(44,29011);
   h->SetBinContent(45,28566);
   h->SetBinContent(46,28091);
   h->SetBinContent(47,27913);
   h->SetBinContent(48,27794);
   h->SetBinContent(49,27546);
   h->SetBinContent(50,27303);
   h->SetBinContent(51,26931);
   h->SetBinContent(52,26776);
   h->SetBinContent(53,26701);
   h->SetBinContent(54,26473);
   h->SetBinContent(55,26359);
   h->SetBinContent(56,26085);
   h->SetBinContent(57,26196);
   h->SetBinContent(58,25648);
   h->SetBinContent(59,25641);
   h->SetBinContent(60,25543);
   h->SetBinContent(61,25339);
   h->SetBinContent(62,25140);
   h->SetBinContent(63,24896);
   h->SetBinContent(64,24702);
   h->SetBinContent(65,24438);
   h->SetBinContent(66,24752);
   h->SetBinContent(67,24434);
   h->SetBinContent(68,24578);
   h->SetBinContent(69,24265);
   h->SetBinContent(70,24162);
   h->SetBinContent(71,23840);
   h->SetBinContent(72,23969);
   h->SetBinContent(73,23760);
   h->SetBinContent(74,23873);
   h->SetBinContent(75,23786);
   h->SetBinContent(76,23730);
   h->SetBinContent(77,23869);
   h->SetBinContent(78,23424);
   h->SetBinContent(79,23548);
   h->SetBinContent(80,23885);
   h->SetBinContent(81,23573);
   h->SetBinContent(82,23718);
   h->SetBinContent(83,23612);
   h->SetBinContent(84,23686);
   h->SetBinContent(85,23652);
   h->SetBinContent(86,23693);
   h->SetBinContent(87,23642);
   h->SetBinContent(88,23708);
   h->SetBinContent(89,23516);
   h->SetBinContent(90,23831);
   h->SetBinContent(91,23910);
   h->SetBinContent(92,23940);
   h->SetBinContent(93,24107);
   h->SetBinContent(94,23894);
   h->SetBinContent(95,23890);
   h->SetBinContent(96,24256);
   h->SetBinContent(97,24147);
   h->SetBinContent(98,24406);
   h->SetBinContent(99,24298);
   h->SetBinContent(100,24642);
   h->SetBinContent(101,24878);
   h->SetBinContent(102,24610);
   h->SetBinContent(103,24984);
   h->SetBinContent(104,25273);
   h->SetBinContent(105,24875);
   h->SetBinContent(106,25294);
   h->SetBinContent(107,25299);
   h->SetBinContent(108,25739);
   h->SetBinContent(109,25904);
   h->SetBinContent(110,25805);
   h->SetBinContent(111,26085);
   h->SetBinContent(112,25883);
   h->SetBinContent(113,26699);
   h->SetBinContent(114,26593);
   h->SetBinContent(115,26796);
   h->SetBinContent(116,26998);
   h->SetBinContent(117,26933);
   h->SetBinContent(118,27569);
   h->SetBinContent(119,27326);
   h->SetBinContent(120,27631);
   h->SetBinContent(121,27541);
   h->SetBinContent(122,28191);
   h->SetBinContent(123,28298);
   h->SetBinContent(124,28361);
   h->SetBinContent(125,28290);
   h->SetBinContent(126,28663);
   h->SetBinContent(127,28523);
   h->SetBinContent(128,28835);
   h->SetBinContent(129,29063);
   h->SetBinContent(130,29360);
   h->SetBinContent(131,29467);
   h->SetBinContent(132,29289);
   h->SetBinContent(133,29450);
   h->SetBinContent(134,29784);
   h->SetBinContent(135,30457);
   h->SetBinContent(136,30017);
   h->SetBinContent(137,30257);
   h->SetBinContent(138,30538);
   h->SetBinContent(139,30437);
   h->SetBinContent(140,30773);
   h->SetBinContent(141,30565);
   h->SetBinContent(142,30593);
   h->SetBinContent(143,30692);
   h->SetBinContent(144,30388);
   h->SetBinContent(145,30170);
   h->SetBinContent(146,30772);
   h->SetBinContent(147,30024);
   h->SetBinContent(148,30295);
   h->SetBinContent(149,29978);
   h->SetBinContent(150,29543);
   h->SetBinContent(151,29558);
   h->SetBinContent(152,29124);
   h->SetBinContent(153,28821);
   h->SetBinContent(154,27943);
   h->SetBinContent(155,27337);
   h->SetBinContent(156,27159);
   h->SetBinContent(157,26752);
   h->SetBinContent(158,25648);
   h->SetBinContent(159,24788);
   h->SetBinContent(160,23766);
   h->SetBinContent(161,23374);
   h->SetBinContent(162,22605);
   h->SetBinContent(163,21313);
   h->SetBinContent(164,20309);
   h->SetBinContent(165,19242);
   h->SetBinContent(166,18230);
   h->SetBinContent(167,17204);
   h->SetBinContent(168,15985);
   h->SetBinContent(169,14903);
   h->SetBinContent(170,13582);
   h->SetBinContent(171,12393);
   h->SetBinContent(172,11071);
   h->SetBinContent(173,10018);
   h->SetBinContent(174,8905);
   h->SetBinContent(175,7754);
   h->SetBinContent(176,6872);
   h->SetBinContent(177,5797);
   h->SetBinContent(178,4954);
   h->SetBinContent(179,4157);
   h->SetBinContent(180,3399);
   h->SetBinContent(181,2766);
   h->SetBinContent(182,2272);
   h->SetBinContent(183,1848);
   h->SetBinContent(184,1370);
   h->SetBinContent(185,1043);
   h->SetBinContent(186,809);
   h->SetBinContent(187,621);
   h->SetBinContent(188,431);
   h->SetBinContent(189,298);
   h->SetBinContent(190,208);
   h->SetBinContent(191,144);
   h->SetBinContent(192,112);
   h->SetBinContent(193,78);
   h->SetBinContent(194,35);
   h->SetBinContent(195,28);
   h->SetBinContent(196,14);
   h->SetBinContent(197,11);
   h->SetBinContent(198,10);
   h->SetBinContent(199,1);
   h->SetBinContent(200,6);
   h->SetBinContent(201,3);
}

void fillPileHist_6mm_raw(TH1F *hpile6){
   hpile6->SetBinContent(2,3.32061e+10);
   hpile6->SetBinContent(3,1.353877e+11);
   hpile6->SetBinContent(4,3.233289e+11);
   hpile6->SetBinContent(5,5.153087e+11);
   hpile6->SetBinContent(6,7.747058e+11);
   hpile6->SetBinContent(7,1.074594e+12);
   hpile6->SetBinContent(8,1.422394e+12);
   hpile6->SetBinContent(9,1.786299e+12);
   hpile6->SetBinContent(10,2.101183e+12);
   hpile6->SetBinContent(11,2.473524e+12);
   hpile6->SetBinContent(12,2.992068e+12);
   hpile6->SetBinContent(13,3.372812e+12);
   hpile6->SetBinContent(14,3.840925e+12);
   hpile6->SetBinContent(15,4.471071e+12);
   hpile6->SetBinContent(16,4.763684e+12);
   hpile6->SetBinContent(17,5.481582e+12);
   hpile6->SetBinContent(18,5.659574e+12);
   hpile6->SetBinContent(19,6.184112e+12);
   hpile6->SetBinContent(20,6.519637e+12);
   hpile6->SetBinContent(21,7.083127e+12);
   hpile6->SetBinContent(22,7.084223e+12);
   hpile6->SetBinContent(23,8.00838e+12);
   hpile6->SetBinContent(24,8.092882e+12);
   hpile6->SetBinContent(25,8.237216e+12);
   hpile6->SetBinContent(26,8.596731e+12);
   hpile6->SetBinContent(27,8.997363e+12);
   hpile6->SetBinContent(28,9.430463e+12);
   hpile6->SetBinContent(29,9.789213e+12);
   hpile6->SetBinContent(30,1.037643e+13);
   hpile6->SetBinContent(31,1.046791e+13);
   hpile6->SetBinContent(32,1.05557e+13);
   hpile6->SetBinContent(33,1.085701e+13);
   hpile6->SetBinContent(34,1.059112e+13);
   hpile6->SetBinContent(35,1.135137e+13);
   hpile6->SetBinContent(36,1.124423e+13);
   hpile6->SetBinContent(37,1.159906e+13);
   hpile6->SetBinContent(38,1.185982e+13);
   hpile6->SetBinContent(39,1.240347e+13);
   hpile6->SetBinContent(40,1.227658e+13);
   hpile6->SetBinContent(41,1.242497e+13);
   hpile6->SetBinContent(42,1.259935e+13);
   hpile6->SetBinContent(43,1.291336e+13);
   hpile6->SetBinContent(44,1.327274e+13);
   hpile6->SetBinContent(45,1.346249e+13);
   hpile6->SetBinContent(46,1.381326e+13);
   hpile6->SetBinContent(47,1.389908e+13);
   hpile6->SetBinContent(48,1.432294e+13);
   hpile6->SetBinContent(49,1.41195e+13);
   hpile6->SetBinContent(50,1.440296e+13);
   hpile6->SetBinContent(51,1.481504e+13);
   hpile6->SetBinContent(52,1.48076e+13);
   hpile6->SetBinContent(53,1.573524e+13);
   hpile6->SetBinContent(54,1.539071e+13);
   hpile6->SetBinContent(55,1.568738e+13);
   hpile6->SetBinContent(56,1.605377e+13);
   hpile6->SetBinContent(57,1.649786e+13);
   hpile6->SetBinContent(58,1.685123e+13);
   hpile6->SetBinContent(59,1.698378e+13);
   hpile6->SetBinContent(60,1.746162e+13);
   hpile6->SetBinContent(61,1.77953e+13);
   hpile6->SetBinContent(62,1.855595e+13);
   hpile6->SetBinContent(63,1.873638e+13);
   hpile6->SetBinContent(64,1.857268e+13);
   hpile6->SetBinContent(65,1.915779e+13);
   hpile6->SetBinContent(66,1.929335e+13);
   hpile6->SetBinContent(67,2.006321e+13);
   hpile6->SetBinContent(68,2.025582e+13);
   hpile6->SetBinContent(69,2.101979e+13);
   hpile6->SetBinContent(70,2.110328e+13);
   hpile6->SetBinContent(71,2.144166e+13);
   hpile6->SetBinContent(72,2.195509e+13);
   hpile6->SetBinContent(73,2.3048e+13);
   hpile6->SetBinContent(74,2.259514e+13);
   hpile6->SetBinContent(75,2.347694e+13);
   hpile6->SetBinContent(76,2.37515e+13);
   hpile6->SetBinContent(77,2.426597e+13);
   hpile6->SetBinContent(78,2.432384e+13);
   hpile6->SetBinContent(79,2.481406e+13);
   hpile6->SetBinContent(80,2.483657e+13);
   hpile6->SetBinContent(81,2.518641e+13);
   hpile6->SetBinContent(82,2.543427e+13);
   hpile6->SetBinContent(83,2.518356e+13);
   hpile6->SetBinContent(84,2.605613e+13);
   hpile6->SetBinContent(85,2.573869e+13);
   hpile6->SetBinContent(86,2.60914e+13);
   hpile6->SetBinContent(87,2.514113e+13);
   hpile6->SetBinContent(88,2.505322e+13);
   hpile6->SetBinContent(89,2.493685e+13);
   hpile6->SetBinContent(90,2.397308e+13);
   hpile6->SetBinContent(91,2.458824e+13);
   hpile6->SetBinContent(92,2.414158e+13);
   hpile6->SetBinContent(93,2.325187e+13);
   hpile6->SetBinContent(94,2.273504e+13);
   hpile6->SetBinContent(95,2.260052e+13);
   hpile6->SetBinContent(96,2.17065e+13);
   hpile6->SetBinContent(97,2.148578e+13);
   hpile6->SetBinContent(98,2.095296e+13);
   hpile6->SetBinContent(99,2.069695e+13);
   hpile6->SetBinContent(100,2.061577e+13);
   hpile6->SetBinContent(101,1.994671e+13);
   hpile6->SetBinContent(102,1.922501e+13);
   hpile6->SetBinContent(103,1.860626e+13);
   hpile6->SetBinContent(104,1.833033e+13);
   hpile6->SetBinContent(105,1.790257e+13);
   hpile6->SetBinContent(106,1.766996e+13);
   hpile6->SetBinContent(107,1.774769e+13);
   hpile6->SetBinContent(108,1.786128e+13);
   hpile6->SetBinContent(109,1.671002e+13);
   hpile6->SetBinContent(110,1.666938e+13);
   hpile6->SetBinContent(111,1.616491e+13);
   hpile6->SetBinContent(112,1.628043e+13);
   hpile6->SetBinContent(113,1.606833e+13);
   hpile6->SetBinContent(114,1.556144e+13);
   hpile6->SetBinContent(115,1.546141e+13);
   hpile6->SetBinContent(116,1.510769e+13);
   hpile6->SetBinContent(117,1.481468e+13);
   hpile6->SetBinContent(118,1.456828e+13);
   hpile6->SetBinContent(119,1.438687e+13);
   hpile6->SetBinContent(120,1.38397e+13);
   hpile6->SetBinContent(121,1.387284e+13);
   hpile6->SetBinContent(122,1.427458e+13);
   hpile6->SetBinContent(123,1.364349e+13);
   hpile6->SetBinContent(124,1.347135e+13);
   hpile6->SetBinContent(125,1.303305e+13);
   hpile6->SetBinContent(126,1.305059e+13);
   hpile6->SetBinContent(127,1.237641e+13);
   hpile6->SetBinContent(128,1.263388e+13);
   hpile6->SetBinContent(129,1.220069e+13);
   hpile6->SetBinContent(130,1.203169e+13);
   hpile6->SetBinContent(131,1.15933e+13);
   hpile6->SetBinContent(132,1.162997e+13);
   hpile6->SetBinContent(133,1.098271e+13);
   hpile6->SetBinContent(134,1.077927e+13);
   hpile6->SetBinContent(135,1.02592e+13);
   hpile6->SetBinContent(136,1.022483e+13);
   hpile6->SetBinContent(137,9.881253e+12);
   hpile6->SetBinContent(138,9.820153e+12);
   hpile6->SetBinContent(139,9.375477e+12);
   hpile6->SetBinContent(140,9.080573e+12);
   hpile6->SetBinContent(141,8.660575e+12);
   hpile6->SetBinContent(142,8.191987e+12);
   hpile6->SetBinContent(143,8.029066e+12);
   hpile6->SetBinContent(144,7.645037e+12);
   hpile6->SetBinContent(145,7.520915e+12);
   hpile6->SetBinContent(146,6.926998e+12);
   hpile6->SetBinContent(147,6.945928e+12);
   hpile6->SetBinContent(148,6.153391e+12);
   hpile6->SetBinContent(149,5.752735e+12);
   hpile6->SetBinContent(150,5.488233e+12);
   hpile6->SetBinContent(151,5.392067e+12);
   hpile6->SetBinContent(152,4.678201e+12);
   hpile6->SetBinContent(153,4.548931e+12);
   hpile6->SetBinContent(154,3.9761e+12);
   hpile6->SetBinContent(155,3.599476e+12);
   hpile6->SetBinContent(156,3.342955e+12);
   hpile6->SetBinContent(157,2.996012e+12);
   hpile6->SetBinContent(158,2.571277e+12);
   hpile6->SetBinContent(159,2.262419e+12);
   hpile6->SetBinContent(160,1.938952e+12);
   hpile6->SetBinContent(161,1.570958e+12);
   hpile6->SetBinContent(162,1.376935e+12);
   hpile6->SetBinContent(163,1.099553e+12);
   hpile6->SetBinContent(164,8.851618e+11);
   hpile6->SetBinContent(165,7.276967e+11);
   hpile6->SetBinContent(166,5.779174e+11);
   hpile6->SetBinContent(167,4.126839e+11);
   hpile6->SetBinContent(168,3.266407e+11);
   hpile6->SetBinContent(169,2.087774e+11);
   hpile6->SetBinContent(170,1.525711e+11);
   hpile6->SetBinContent(171,9.35151e+10);
   hpile6->SetBinContent(172,5.338599e+10);
   hpile6->SetBinContent(173,2.86137e+10);
   hpile6->SetBinContent(174,1.841288e+10);
   hpile6->SetBinContent(175,9.650183e+09);
   hpile6->SetBinContent(176,3.523847e+09);
   hpile6->SetBinContent(177,1.445073e+09);
   hpile6->SetBinContent(178,4.695789e+08);
   hpile6->SetBinContent(179,1.04723e+08);
   hpile6->SetBinContent(180,2.156298e+07);
   hpile6->SetBinContent(181,2257322);
   hpile6->SetBinContent(182,146537);
   hpile6->SetBinContent(183,7006);
   hpile6->SetBinContent(184,86);
}

void interpPileHist(TH1F *hPile, TH1F *hInt, Int_t nBins, Int_t nBinsCorr){
  for(Int_t i = 1; i <= nBinsCorr; i++){
    Float_t binCenter = hInt->GetBinCenter(i);
    Int_t pileBinNum = 1;
    Float_t pileBinCenter = 0.0;
    while(hPile->GetBinCenter(pileBinNum) < binCenter && pileBinNum <= nBins){
      pileBinNum++;
    }
    Float_t normBin = hPile->GetBinContent(pileBinNum);
    Float_t lowBin = hPile->GetBinContent(pileBinNum - 1);
    Float_t normBinCenter = hPile->GetBinCenter(pileBinNum);
    Float_t lowBinCenter = hPile->GetBinCenter(pileBinNum - 1);
    Float_t slope = (normBin - lowBin)/(normBinCenter - lowBinCenter);
    Float_t xDiff = normBinCenter - binCenter;
    hInt->SetBinContent(i, normBin - slope*xDiff);
  }
}

void fillSimHist_6mm(THStack *hs, TLegend *leg, Int_t color, Float_t matchFit=0.0){
  Float_t upperBound = 170;
  Float_t pileUpperBound = 350;
  Int_t nBins = 200;

  TH1F *h6 = new TH1F("h6", "Basic Simulation Histogram", nBins, 0, upperBound*1.01);
  TH1F *hPile6 = new TH1F("hPile6", "Basic Simulation Histogram", nBins, 0, pileUpperBound);
  fillSimHist_6mm_raw(h6);
  fillPileHist_6mm_raw(hPile6);

  Float_t ce = findComptonEdgeSim(h6, nBins);
  Int_t nBinsCorr = (Int_t)nBins*pileUpperBound/upperBound;
  TH1F *h6_corr = new TH1F("h6_corr", "Scaled Simulation Histogram", nBinsCorr, 0, pileUpperBound*1.0/ce);
  TH1F *hPile6_corr = new TH1F("hPile6_corr", "Scaled Simulation Histogram", nBins, 0, pileUpperBound*1.0/ce);
  TH1F *hPile6_int = new TH1F("hPile6_int", "Scaled Simulation Histogram", nBinsCorr, 0, pileUpperBound*1.0/ce);
  fillSimHist_6mm_raw(h6_corr);
  fillPileHist_6mm_raw(hPile6_corr);
  interpPileHist(hPile6_corr, hPile6_int, nBins, nBinsCorr);
  
  Float_t pileScale = getNormFactorPile(hPile6_int, nBinsCorr, matchFit);
  hPile6_int->Scale(pileScale);

  Float_t scale = getNormFactorSim(h6_corr, nBinsCorr);
  h6_corr->Scale(1.0/scale);
  h6_corr->SetLineColor(color);
  h6_corr->Add(hPile6_int);
  scale = getNormFactorSim(h6_corr, nBinsCorr);
  h6_corr->Scale(1.0/scale);

  TCanvas *c = new TCanvas("cPileupTest", "Pileup Text Canvas", 1200, 800);
  c->cd();
  h6_corr->Draw();
  
  hs->Add(h6_corr);
  leg->AddEntry(h6_corr, "6 mm offset sim");
}

void plotPileHist_6mm(){
  TCanvas *c = new TCanvas("cPile_6mm", "Pileup Canvas 6mm", 1200, 800);
  TH1F *hPile6 = new TH1F("hPile6", "CREX 6mm Offset Pileup Canvas", 200, 0, 350);
  fillPileHist_6mm_raw(hPile6);

  c->cd();
  hPile6->Draw();
}
