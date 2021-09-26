//
//  fscell.cpp
//  cortex
//
//  Created by Diogo Porfirio de Castro Vieira on 27/11/11.
//  Copyright 2011 Universidade de Sao Paulo. All rights reserved.
//

#include "fscell.h"  

fscell::fscell() : neuron() {
    numeq = 2;
    fbuf = new double[numeq];
    w = new double[numeq];
    waux = new double[numeq];
    ks = new double * [numeq]; 
    for (int i = 0; i < numeq; i++)
    	ks[i] = new double[4];         
    startKs();
    cap = 20.0, k = 1.0, vr = -55.0, vt = -40.0, vb = -55.0;
    setpar(0.2, 0.025, -45.0, 0.0);
    seth(0.1);
	double w_[2] = {vr,0.0};
    setw0(w_);        
    setInhibitory();
}    

fscell::fscell(double aArg, double bArg, double cArg, double dArg) : neuron() {
    numeq = 2;
    fbuf = new double[numeq];
    w = new double[numeq];
    waux = new double[numeq];
    ks = new double * [numeq]; 
    for (int i = 0; i < numeq; i++)
    	ks[i] = new double[4];         
    startKs();
    cap = 20.0, k = 1.0, vr = -55.0, vt = -40.0, vb = -55.0;
    setpar(aArg, bArg, cArg, dArg);
    seth(0.1);
	double w_[2] = {vr,0.0};
    setw0(w_);          
    setInhibitory();
}

void fscell::fx(double inj, double time) {
    inj -= calcsyncurrent(time);
    fbuf[0] = (k * (waux[0] - vr) * (waux[0] - vt) - waux[1] + inj)/cap;
    if (w[0] < vb)
        fbuf[1] = - a * waux[1];
    else
        fbuf[1] = a * (b * pow((waux[0] - vb),3) - waux[1]);
}

void fscell::checkPeak(double time) {
    if (w[0] >= 25.0) {
        w[0] = c;
        sendevent(time);
        backevent();
    }
}

void fscell::setpar(double aArg, double bArg, double cArg, double dArg) {
    a = aArg;
    b = bArg;
    c = cArg;
    d = dArg;
}

void fscell::setvrest(double vrest){
    vr = vrest;
	double w_[2] = {vrest,0.0};
    setw0(w_);
}

void fscell::setvtresh(double vtresh){
    vt = vtresh;
}

void fscell::setk(double k_){
    k = k_;
}

void fscell::setcap(double cap_){
    cap = cap_;
}
