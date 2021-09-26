/* 
 * File:   alphasynapse.h
 * Author: diogopcv
 *
 * Created on 4 de Fevereiro de 2010, 13:24
 */

#ifndef _ALPHASYNAPSE_H
#define	_ALPHASYNAPSE_H

#include <vector>
#include <math.h>
#include "neuron.h"

namespace Cells{
	class neuron;}
using namespace std;
using namespace Cells;

namespace synapse {
	
    class alphasynapse {
    public:
	    alphasynapse();
	    alphasynapse(double tauArg, double gmaxArg, double delayArg, double EsynArg);
	    alphasynapse(double tauArg, double gmaxArg, double delayArg, double EsynArg, bool shortArg, bool longArg);
	    void setpar(double tauArg, double gmaxArg, double delayArg, double EsynArg);
	    void setShortTerm(int argShort, double ap, double ataux);
	    void evaluate(double time);
	    void addevent(double spk);
	    double getE();
	    double getGsyn();
	    double getGmax();
	    void seth(double hArg);   
	    void setLongTerm(int argLong, double Aminus, double Aplus, double mu, double tauMinxus, double tauPlus);
	    void refreshM();
	    double getXfac();
	    double getPa(); 
	    double getM(); 
	    void calcK1();
	    void calcK2();
	    void calcK3();    
	    void calcK4(); 
	    double getGaux();
	protected:    	  
	    vector <double> * spikes;  
	    void startKs();     
	    double h;
	    double * k1, * k2, * k3, * k4;
	    double tau, delay, gmax, gsyn, Esyn;
	    bool SHORTERM, LONGTERM;
	    double xfac, taux, p, Pa, M;
	    double Aminus, Aplus, tauMinxus, tauPlus, mu;   
	    double gsynaux, xfacaux, paux, Paaux, Maux;  
    };
}

#endif	/* _ALPHASYNAPSE_H */

