#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>

#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"
#include "TLatex.h"
#include "TH1.h"
#include "TH1F.h"


std::string strip(std::string line);
std::string shell(std::string cmd);
std::string replace(std::string& str, std::string oldstr, std::string newstr);
void split(std::string str, std::vector<std::string>& vstr);
void error(std::string message);
void setStyle();


TH1F* mkhist1(std::string hname, 
              std::string xtitle, std::string ytitle, 
              int nbins,             
              double* bins,
              int color=kBlue);

TH1F* mkhist1(std::string hname, 
              std::string xtitle, std::string ytitle, 
              std::vector<double>& bins,
              int color);

TH1F* mkhist1(std::string hname, 
              std::string xtitle, std::string ytitle, 
              int nbins, double xmin, double xmax,
              int color);



std::vector<double> contents(TH1* hist);
std::vector<double> binlowedges(TH1* hist);
std::vector<double> binhighedges(TH1* hist);
std::vector<double> bincenters(TH1* hist);
std::vector<double> contents(TH1* hist);
std::vector<double> errors(TH1* hist);
std::vector<double> cdf(TH1* hist);
void setContents(TH1* hist, std::vector<double>& c);
void setErrors(TH1* hist, std::vector<double>& err);

/// Simple wrapper around TLatex that uses NDC coordinates
class Scribe
{
public:
  Scribe();

  Scribe(float xpos,
         float ypos,
         float size=0.04, 
         int font=42);
  ~Scribe();

  void write(std::string text, float xoffset=0);
  void vspace(float f=0.5);

private:
  float xpos;
  float ypos;
  float ymin;
  float ymax;
  bool  logy;
  float textsize;
  float scale;
  float linewidth;
  int   line;
  TLatex* t;
};

#endif
