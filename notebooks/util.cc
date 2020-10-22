// ----------------------------------------------------------------------------
// Some simple histogramming utilities
// 2012 ESHEP
// Harrison B. Prosper
// ----------------------------------------------------------------------------
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <stdlib.h>

#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "TStyle.h"
#include "TLatex.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TMinuit.h"
#include "util.h"

#include "Math/Interpolator.h"
#include "Math/Chebyshev.h"
#include "Math/IFunction.h"
#include "Math/WrappedFunction.h"
#include "Math/Integrator.h"
#include "Math/RootFinder.h"
// ----------------------------------------------------------------------------
using namespace std;
using namespace ROOT::Math;

namespace {
  const int LINEWIDTH= 1;
  const int TEXTFONT = 42;
  const int NDIV     = 510;
  const int WIDTH    = 500;
  const int HEIGHT   = 500;
  const double TEXTSIZE = 0.04;
  const double MARKERSIZE = 0.8;
}
// ----------------------------------------------------------------------------

void working()
{
  cout << ".";
  cout.flush();
}

string 
strip(string line)
{
  int l = line.size();
  if ( l == 0 ) return string("");
  int n = 0;
  while (((line[n] == 0)    ||
          (line[n] == ' ' ) ||
          (line[n] == '\n') ||
          (line[n] == '\t')) && n < l) n++;
  
  int m = l-1;
  while (((line[m] == 0)    ||
          (line[m] == ' ')  ||
          (line[m] == '\n') ||
          (line[m] == '\t')) && m > 0) m--;
  return line.substr(n,m-n+1);
}

string 
shell(string cmd)
{
  FILE* f = popen(cmd.c_str(),"r");
  int buffsize=8192;
  char s[8192];
  int n = fread(s,1,buffsize,f);
  pclose(f);
  string result = strip(string(s).substr(0,n));
  return result;
}

void 
split(string str, vector<string>& vstr)
{
  vstr.clear();
  istringstream stream(str);
  while ( stream )
    {
      string stringy;
      stream >> stringy;
      if ( stream ) vstr.push_back(stringy);
    }
}

string
replace(string& str, string oldstr, string newstr)
{
  return string(TString(str).ReplaceAll(oldstr, newstr).Data());
}

void error(string message)
{
  cout << "**ERR** " << message << endl;
  exit(0);
}


//----------------------------------------------------------------------------
// Style Utilities
//----------------------------------------------------------------------------

void 
setStyle()
{
  TStyle* style = new TStyle("ESHEP12","ESHEP12");
   
  // For the canvas:
  style->SetCanvasBorderMode(0);
  style->SetCanvasColor(kWhite);
  style->SetCanvasDefH(500); //Height of canvas
  style->SetCanvasDefW(500); //Width of canvas
  style->SetCanvasDefX(0);   //Position on screen
  style->SetCanvasDefY(0);

  // For the Pad:
  style->SetPadBorderMode(0);
  style->SetPadColor(kWhite);
  style->SetPadGridX(kFALSE);
  style->SetPadGridY(kFALSE);
  style->SetGridColor(kGreen);
  style->SetGridStyle(3);
  style->SetGridWidth(1);
    
  // For the frame:
  style->SetFrameBorderMode(0);
  style->SetFrameBorderSize(1);
  style->SetFrameFillColor(0);
  style->SetFrameFillStyle(0);
  style->SetFrameLineColor(1);
  style->SetFrameLineStyle(1);
  style->SetFrameLineWidth(1);
    
  // For the histo:
  style->SetHistLineColor(1);
  style->SetHistLineStyle(0);
  style->SetHistLineWidth(2);

  style->SetEndErrorSize(2);
  style->SetErrorX(0.);
    
  style->SetMarkerSize(0.5);
  style->SetMarkerStyle(20);

  //For the fit/function:
  style->SetOptFit(1);
  style->SetFitFormat("5.4g");
  style->SetFuncColor(2);
  style->SetFuncStyle(1);
  style->SetFuncWidth(1);

  //For the date:
  style->SetOptDate(0);

  // For the statistics box:
  style->SetOptFile(0);
  style->SetOptStat("");
  // To display the mean and RMS:
  //style->SetOptStat("mr"); 
  style->SetStatColor(kWhite);
  style->SetStatFont(TEXTFONT);
  style->SetStatFontSize(0.03);
  style->SetStatTextColor(1);
  style->SetStatFormat("6.4g");
  style->SetStatBorderSize(1);
  style->SetStatH(0.2);
  style->SetStatW(0.3);
    
  // Margins:
  style->SetPadTopMargin(0.05);
  style->SetPadBottomMargin(0.16);
  style->SetPadLeftMargin(0.20);
  style->SetPadRightMargin(0.10);

  // For the Global title:
  style->SetOptTitle(0); 
  style->SetTitleFont(TEXTFONT);
  style->SetTitleColor(1);
  style->SetTitleTextColor(1);
  style->SetTitleFillColor(10);
  style->SetTitleFontSize(0.05);

  // For the axis titles:
  style->SetTitleColor(1, "XYZ");
  style->SetTitleFont(TEXTFONT, "XYZ");
  style->SetTitleSize(0.05, "XYZ");
  style->SetTitleXOffset(1.25);    //(0.9);
  style->SetTitleYOffset(1.60);    //(1.25);

  // For the axis labels:
  style->SetLabelColor(1, "XYZ");
  style->SetLabelFont(TEXTFONT, "XYZ");
  style->SetLabelOffset(0.007, "XYZ");
  style->SetLabelSize(0.05, "XYZ");

  // For the axis:
  style->SetAxisColor(1, "XYZ");
  style->SetStripDecimals(kTRUE);
  style->SetTickLength(0.03, "XYZ");
  style->SetNdivisions(504, "XYZ");
  // To get tick marks on the opposite side of the frame
  style->SetPadTickX(1);  
  style->SetPadTickY(1);

  // Change for log plots:
  style->SetOptLogx(0);
  style->SetOptLogy(0);
  style->SetOptLogz(0);

  // Postscript options:
  style->SetPaperSize(20.,20.);
  style->cd();
}


TH1F* mkhist1(string hname, 
              string xtitle, string ytitle, 
              int nbins, 
              double* bins,
              int color)
{
  TH1F* h = new TH1F(hname.c_str(), "", nbins, bins);

  // Set some reasonable defaults

  h->SetLineColor(color);
  h->SetMarkerSize(1);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(20);

  //h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitle(xtitle.c_str());
  h->GetXaxis()->SetTitleOffset(1.3);
  h->SetNdivisions(504, "X");
  h->SetMarkerSize(1.0);

  //h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitle(ytitle.c_str());
  h->GetYaxis()->SetTitleOffset(1.4); // 1.8
  h->SetNdivisions(504, "Y");

  return h;
}

TH1F* mkhist1(string hname, 
              string xtitle, string ytitle, 
              vector<double>& bins,
              int color)
{
  return mkhist1(hname, xtitle, ytitle, bins.size()-1, &bins[0], color);
}

TH1F* mkhist1(string hname, 
              string xtitle, string ytitle, 
              int nbins, double xmin, double xmax,
              int color)
{
  if ( nbins <= 0 )
    {
      cout << "*** hey! mbins needs to be > 0! " << endl;
      exit(0);
    }
  vector<double> bins(nbins+1);
  double step = (xmax-xmin)/nbins;
  for(int i=0; i <= nbins; i++) bins[i] = xmin + i*step;
  return mkhist1(hname, xtitle, ytitle, nbins, &bins[0], color);
}

Scribe::Scribe()
  : xpos(0),
    ypos(0),
    textsize(0),
    scale(0),
    linewidth(0),
    line(0),
    t(0)
  {}

Scribe::Scribe(float xxpos,
               float yypos,
               float size,
               int font) 
  : xpos(xxpos),
    ypos(yypos),
    textsize(size),
    scale(1.5),
    linewidth(scale*textsize),
    line(0),
    t(new TLatex())
{    
  t->SetNDC();
  t->SetTextSize(textsize);
  t->SetTextFont(font);
  t->SetTextAlign(12);
}

Scribe::~Scribe() { if ( t != 0 ) delete t; }

void Scribe::write(std::string text, float xoffset)
{
  float y = ypos;
  if ( y < 0 ) return;
  t->DrawLatex(xpos+xoffset, y, text.c_str());
  ypos -= linewidth;
  line++;
}

void Scribe::vspace(float f)
{
  float y = ypos;
  if ( y < 0 ) return;
  t->DrawLatex(xpos, y, " ");
  ypos -= linewidth * f;
  line++;
}

vector<double> contents(TH1* hist)
{
  vector<double> c(hist->GetNbinsX());
  for(int i=0; i < hist->GetNbinsX(); i++)
    c[i] = hist->GetBinContent(i+1);
  return c;
}

vector<double> binlowedges(TH1* hist)
{
  vector<double> c(hist->GetNbinsX());
  for(int i=0; i < hist->GetNbinsX(); i++)
    c[i] = hist->GetBinLowEdge(i+1);
  return c;
}

vector<double> binhighedges(TH1* hist)
{
  vector<double> c(hist->GetNbinsX());
  for(int i=0; i < hist->GetNbinsX(); i++)
    c[i] = hist->GetBinLowEdge(i+1) + 0.5*hist->GetBinWidth(i+1);
  return c;
}

vector<double> bincenters(TH1* hist)
{
  vector<double> c(hist->GetNbinsX());
  for(int i=0; i < hist->GetNbinsX(); i++)
    c[i] = hist->GetBinLowEdge(i+1) + 0.5*hist->GetBinWidth(i+1);
  return c;
}

vector<double> errors(TH1* hist)
{
  vector<double> c(hist->GetNbinsX());
  for(int i=0; i < hist->GetNbinsX(); i++)
    c[i] = hist->GetBinError(i+1);
  return c;
}

vector<double> cdf(TH1* hist)
{
  vector<double> c(hist->GetNbinsX());
  c[0] = hist->GetBinContent(1);
  for(int i=1; i < hist->GetNbinsX(); i++)
    c[i] = c[i-1]+hist->GetBinContent(i+1);
  return c;
}

void setContents(TH1* hist, vector<double>& c)
{
  int n = hist->GetNbinsX();
  int nbin = n < (int)c.size() ? n : c.size();
  for(int i=0; i < nbin; i++) hist->SetBinContent(i+1, c[i]);
}

void setErrors(TH1* hist, vector<double>& err)
{
  int n = hist->GetNbinsX();
  int nbin = n < (int)err.size() ? n : err.size();
  for(int i=0; i < nbin; i++) hist->SetBinError(i+1, err[i]);
}

