//
//  thalamicell.h
//  cortex
//
//  Created by Diogo Porfirio de Castro Vieira on 27/11/11.
//  Copyright 2011 Universidade de Sao Paulo. All rights reserved.
//

#ifndef cortex_thalamicell_h
#define cortex_thalamicell_h

#include "neuron.h"

namespace synapse{
	class alphasynapse;}
using namespace synapse;
using namespace std;
namespace Cells {
class thalamicell: public neuron {
public:
    thalamicell();
    thalamicell(double aArg, double bArg, double cArg, double dArg);
    void setpar(double aArg, double bArg, double cArg, double dArg);        
    void setvrest(double vrest);
    void setvtresh(double vtresh);
    void setk(double k_);
    void setcap(double cap_);     
    
protected:
    void checkPeak(double time);
    void fx(double inj,double time);
    double a, b, c, d;
    double cap, k, vr, vt;
};
}

#endif
