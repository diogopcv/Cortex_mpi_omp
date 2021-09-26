/*
 *  neuron.h
 *  
 *
 *  Created by Diogo Porfirio de Castro Vieira on 02/12/09.
 *  Copyright 2009 Universidade de Sao Paulo. All rights reserved.
 *
 */
#ifndef _NEURON_
#define _NEURON_
#include <iostream>
#include <sstream>
#include <math.h>
#include <string>
#include <complex>
#include <vector>
#include <map>
#include "alphasynapse.h"

namespace synapse{
	class alphasynapse;}
using namespace synapse;
using namespace std;
namespace Cells {
	
    class neuron {
    public:
        neuron();
        ~neuron();     
        void evaluate(double inj,double time);
        void setw0(double * w);
        double getw(int ind);
        void seth(double hArg);
        void setpar(double aArg, double bArg, double cArg, double dArg);
       	void setGroup(int group_);
       	int getGroupID();
        vector<double> * getevents();
        void setExcitatory();
        void setInhibitory();
        int getTypeSyn();   
        void setId(int myId);
        int getID();
        void setSpkBuffer(double ** buffer, int preSpk);
        void addsyndend(alphasynapse * syn, int idPre);
        void recEvent(int idPre, double spkTime);
        
    protected:
        virtual void checkPeak(double time) = 0;    	
        virtual void fx(double inj,double time) = 0;
    	void startKs();  
    	void backevent();      
        void sendevent(double time);
        double calcsyncurrent(double time);
        map<int,alphasynapse *> * sdend;
        vector <double> * events;
        double ** ks;
        double * fbuf;
        double * w, * waux;
        double h;
        int _ID;
        double ** spkBuffer;
        int preSpk;        
        int group;
        int numeq;
        int typesyn;
    };
}
#endif
