//-----------------------------------------------------------------------------
// Top Quark Discovery (1995), D0 Results
// Statistical analysis of results using Root
// See 2012 ESHEP Lectures
// Harrison B. Prosper
//-----------------------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "Math/WrappedFunction.h"
#include "Math/Integrator.h"
#include "Math/RootFinder.h"

#include "TApplication.h"
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "util.h"

using namespace std;
//--------------------------------------------------------------------------
// Hide in anonymous namespace
//--------------------------------------------------------------------------
const  float DEFCL     = 0.68;     // Default confidence level
const  double ABSTOL   = 1.e-9;
const  double RELTOL   = 1.e-6;
const  int    SIZE     = 100;
const  int    RULE     = 3;


struct Study
{  
  int    D;
  double B;
  double dB;
  double CL;                 // Confidence Level

  double Q;
  double K;
  double S;
  
  double gamma;
  double mu;
  double beta;
  

  double XMIN;
  double XMAX;
  double NORM;
 
  ROOT::Math::WrappedMemFunction<Study, double (Study::*)(double)> wfm;
  ROOT::Math::IntegratorOneDim ifn;

  ROOT::Math::WrappedMemFunction<Study, double (Study::*)(double)> wfp;
  ROOT::Math::IntegratorOneDim ifp;

  // D0 Results
  // D = 17 events
  // b = 3.8 +/- 0.6 events

  Study(int d=17, double b=3.8, double db=0.6, double cl=0.683)
    : D(d),
      B(b),
      dB(db),
      CL(cl),
      Q(pow(B/dB,2)),
      K(Q/B),
      S(0),
      
      XMIN(0),
      XMAX(20+d+5*sqrt(d)),
      NORM(1),

      gamma(Q+1),
      mu(0),
      beta(1.0/K),
      
      wfm(*this, &Study::fm),
      ifn(wfm, 
          ROOT::Math::IntegrationOneDim::kADAPTIVE,
          ABSTOL,
          RELTOL,
          SIZE,
          RULE),

      wfp(*this, &Study::marginalExact),
      ifp(wfp, 
          ROOT::Math::IntegrationOneDim::kADAPTIVE,
          ABSTOL,
          RELTOL,
          SIZE,
          RULE)
    {
      NORM = 1;
      NORM = ifp.Integral(XMIN, XMAX);
    }
  
  ~Study() {}
  
  //----------------------------------------------------------------------
  // Define likelihood for D0 results
  //----------------------------------------------------------------------
  double likelihood(double s, double b)
  {
    return TMath::Poisson(D, s + b) * TMath::GammaDist(b, gamma, mu, beta);
  }
  
  //----------------------------------------------------------------------
  // Compute exact marginal likelihood p(D|s) found by integrating 
  // likelihood with respect to b
  //----------------------------------------------------------------------
  double marginalExact(double s)
  {
    double y = pow(K/(1+K), Q+1);
    double z = 1.0/(1+K);
    double sum = y * TMath::Poisson(D, s);
    for(int r=1; r <= D; ++r)
      {
        y *= z * (Q+r)/r;
        sum += y * TMath::Poisson(D-r, s);
      }
    return sum;
  }

  //----------------------------------------------------------------------
  // Compute marginal likelihood p(D|s) found by numerically integrating 
  // likelihood with respect to b
  //----------------------------------------------------------------------
  double marginal(double s)
  {
    S = s;
    return ifn.Integral(XMIN, XMAX);
  }
  
  //----------------------------------------------------------------------
  // Compute profile likelihood p(D|s, bestf(s)) 
  //----------------------------------------------------------------------
  double profile(double s)
  {
    return likelihood(s, bestf(s));
  }

  bool limits(double& xmin, double& xmax)
  {
    // function whose root is to be found
    ROOT::Math::WrappedMemFunction<Study,
      double (Study::*)(double)> fn(*this, &Study::fchisq);
    ROOT::Math::RootFinder rootfinder;
    
    double lo = 0;
    double hi = D - B;
    rootfinder.SetFunction(fn, lo, hi);
    int status = rootfinder.Solve();
    if ( status != 1 )
      {
        cout << "*** Post *** RootFinder failed"
             << endl;
        return false;
      }
    xmin = rootfinder.Root();

    lo = D - B;
    hi = XMAX;
    rootfinder.SetFunction(fn, lo, hi);
    status = rootfinder.Solve();
    if ( status != 1 )
      {
        cout << "*** Post *** RootFinder failed"
             << endl;
        return false;
      }
    xmax = rootfinder.Root();
    return true;
  }
  
  bool blimits(double& xmin, double& xmax)
  {
    double oldCL = CL;
    double alpha = (1-CL)/2;
    CL = alpha;

    // function whose root is to be found
    ROOT::Math::WrappedMemFunction<Study,
      double (Study::*)(double)> fn(*this, &Study::fcdf);
    ROOT::Math::RootFinder rootfinder;
    
    double lo = 0;
    double hi = D - B;
    rootfinder.SetFunction(fn, lo, hi);
    int status = rootfinder.Solve();
    if ( status != 1 )
      {
        cout << "*** Post *** RootFinder failed"
             << endl;
        return false;
      }
    xmin = rootfinder.Root();

    CL = oldCL + alpha;
    lo = D - B;
    hi = XMAX;
    rootfinder.SetFunction(fn, lo, hi);
    status = rootfinder.Solve();
    if ( status != 1 )
      {
        cout << "*** Post *** RootFinder failed"
             << endl;
        return false;
      }
    xmax = rootfinder.Root();

    CL = oldCL;
    return true;
  }
  
  double chisq(double s)
  {
    double num = likelihood(s, bestf(s));
    double shat = D - B;
    double den = likelihood(shat, B);
    return -2*log(num / den);
  }

  double fchisq(double s)
  {
    return chisq(s) - 1;
  }

  double fm(double b)
  {
    double y = likelihood(S, b); 
    return y;
  }
  
  double fcdf(double x)
  {
    return cdf(x) - CL;
  }
  
  double lln(double x)
  {
    return -(1+K) + D/(S + x) + Q/x;
  }
  double cdf(double x) { return ifp.Integral(XMIN, x) / NORM; } 
  
  //----------------------------------------------------------------------
  // Compute numerically, best fit of b for a given s
  //----------------------------------------------------------------------
  double best(double s)
  {
    S = s;
    // function whose root is to be found
    ROOT::Math::WrappedMemFunction<Study,
      double (Study::*)(double)> fn(*this, &Study::lln);
    ROOT::Math::RootFinder rootfinder;
    
    rootfinder.SetFunction(fn, XMIN, XMAX);
    int status = rootfinder.Solve();
    if ( status != 1 )
      {
        cout << "*** Post *** RootFinder failed"
             << endl;
        return -1;
      }
    return rootfinder.Root();
  }
  
  //----------------------------------------------------------------------
  // Compute exactly best fit of b for a given s
  //----------------------------------------------------------------------
  double bestf(double s)
  {
      double q = D + Q - s*(1+K);
      double y = 0.5*(q + sqrt(q*q + 4*(1+K)*s*Q))/(1+K);
      return y;
  }
};

//----------------------------------------------------------------------------

Study* PTR=0;

double profile(double* x, double* p)
{
  return PTR->profile(x[0]);
}

double marginal(double* x, double* p)
{
  return PTR->marginalExact(x[0]);
}

double likelihood(double* x, double* p)
{
  return PTR->likelihood(x[0], x[1]);
}
//----------------------------------------------------------------------------
// MAIN PROGRAM
//----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // set 
  TApplication app("app", &argc, argv);
  setStyle();
  gStyle->SetPalette(1);

  char record[80];

  // D0 top quark discovery results

  int  D = 17;
  double B = 3.8;
  double dB = 0.6;
  Study study(D, B, dB);

  // need pointer to study object
  PTR = &study;

  // histogram bounds
  double xmin=0;
  double xmax=40;

  double ymin=0;
  double ymax=10;


  double x[200];
  double y[200];

  int xbins=80;
  double xstep=(xmax-xmin)/xbins;
  ofstream fout("post.txt");
  sprintf(record, "%10s %10s %10s %10s", "signal", "approx", "exact", "cdf");
  fout << record << endl;

  ofstream bout("best.txt");
  sprintf(record, "%10s %10s %10s", "signal", "approx", "exact");
  bout << record << endl;

  for(int i=0; i <= xbins; ++i)
    {
      double s  = xmin + i*xstep; 
      double y1 = study.marginal(s);         // approximate marginal 
      double y2 = study.marginalExact(s);    // exact marginal
      double y3 = study.cdf(s);              // cumulative dist. function
      sprintf(record, "%10.2f %10.3e %10.3e %10.3e", s, y1, y2, y3);
      fout << record << endl;

      double b1 = study.best(s);
      double b2 = study.bestf(s);
      sprintf(record, "%10.2f %10.3f %10.3f", s, b1, b2);
      bout << record << endl;
      x[i] = s;
      y[i] = b2;
    }
  fout.close();
  bout.close();

  // ---- check a few calculations

  double s  = study.D - study.B;             // best fit signal
  double b1 = study.best(s);                 // best fit background (approx)
  double b2 = study.bestf(s);                // best fit background (exact)
  cout << endl;

  sprintf(record, "%10s %10s %10s", "B", "approx", "exact");
  cout << record << endl;

  sprintf(record, "%10.3f %10.3f %10.3f", study.B, b1, b2);
  cout << record << endl;
  
  // Now compute 68.3% interval

  double lower;
  double upper;
  if ( study.limits(lower, upper) )
    {
      double width = upper - lower;
      sprintf(record, "\t%10.2f, %-10.2f - width = %-10.2f", 
              lower, upper, width);
      cout << endl
           << "interval based on profile likelihood" << endl;
      cout << record << endl;


      double x = sqrt(study.D + study.dB*study.dB);
      lower = s - x;
      upper = s + x;
      width = upper - lower;
      sprintf(record, "\t%10.2f, %-10.2f - width = %-10.2f", 
              lower, upper, width);
      cout << "interval using D - B +/- sqrt(D + db**2)" << endl;
      cout << record << endl;
    }

  if ( study.blimits(lower, upper) )
    {
      double width = upper - lower;
      sprintf(record, "\t%10.2f, %-10.2f - width = %-10.2f", 
              lower, upper, width);
      cout << "central interval using Bayes" << endl;
      cout << record << endl;
    }

  double B10 = study.marginalExact(s) / study.marginalExact(0);
  cout << endl;
  sprintf(record, "B10 = %10.3f, sqrt[2ln(B10)] = %10.3f", B10, 
          sqrt(2*log(B10)));
  cout << record << endl << endl;

  // Margins:
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.22);
  gStyle->SetPadRightMargin(0.08);
  gStyle->SetTitleYOffset(1.70);    //(1.25);

  TCanvas cprofile("fig_profile", "profile", 10, 10, 500, 500);
  TF1 fprofile("profile", profile, xmin, xmax);
  fprofile.SetLineColor(kGreen+1);
  fprofile.SetLineWidth(2);
  fprofile.SetMinimum(0);
  fprofile.SetMaximum(0.1);
  fprofile.GetHistogram()->GetXaxis()->SetTitle("expected signal (s)");
  fprofile.GetHistogram()->GetXaxis()->CenterTitle();
  cprofile.cd();
  fprofile.Draw();
  cprofile.Update();
  cprofile.SaveAs(".gif");

  TCanvas cmarginal("fig_marginal", "profile/marginal", 10, 10, 500, 500);
  TF1 fmarginal("marginal", marginal, xmin, xmax);
  fmarginal.SetLineColor(kRed+1);
  fmarginal.SetLineWidth(2);
  fmarginal.SetMinimum(0);
  fmarginal.SetMaximum(0.1);
  fmarginal.GetHistogram()->GetXaxis()->SetTitle("expected signal (s)");
  fmarginal.GetHistogram()->GetXaxis()->CenterTitle();
  fmarginal.GetHistogram()->GetYaxis()->SetTitle("p(17|s, H_{1})");
  fmarginal.GetHistogram()->GetYaxis()->CenterTitle();
  cmarginal.cd();
  fmarginal.Draw();
  cmarginal.Update();
  cmarginal.SaveAs(".gif");


  // Margins:
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.16);
  gStyle->SetTitleYOffset(1.35);    //(1.25);


  TGraph graph(xbins, x, y);
  graph.SetLineWidth(2);

  TCanvas clike("fig_likelihood", "likelihood", 520, 10, 500, 500);
  TF2 flike("likelihood", likelihood, 0, 30, 0, 6);
  flike.SetLineColor(kRed+1);
  flike.SetLineWidth(2);
  flike.SetMinimum(0);
  flike.GetHistogram()->GetXaxis()->SetTitle("expected signal (s)");
  flike.GetHistogram()->GetXaxis()->CenterTitle();
  flike.GetHistogram()->GetYaxis()->SetTitle("expected background (b)");
  flike.GetHistogram()->GetYaxis()->CenterTitle();
  clike.cd();
  flike.Draw("contz");
  graph.Draw("C");

  double xpos = 0.20;
  double ypos = 0.40;
  Scribe scribe(xpos, ypos);
  scribe.write("top quark discovery (1995)");
  scribe.write("D0 results", 0.05);
  scribe.write("D = 17 events",0.10);
  scribe.write("B = 3.8 #pm 0.6 events",0.10);
  clike.Update();
  clike.SaveAs(".gif");

  app.Run();
  return 0;  
}

