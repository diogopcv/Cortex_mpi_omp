//
//  thalamicell.cpp
//  cortex
//
//  Created by Diogo Porfirio de Castro Vieira on 27/11/11.
//  Copyright 2011 Universidade de Sao Paulo. All rights reserved.
//

#include "thalamicell.h"  

thalamicell::thalamicell() : neuron() {
    numeq = 2;
    fbuf = new double[numeq];
    w = new double[numeq];
    waux = new double[numeq];
    ks = new double * [numeq]; 
    for (int i = 0; i < numeq; i++)
    	ks[i] = new double[4];         
    startKs(); 
    cap = 200.0, k = 1.6, vr = -60.0, vt = -50.0;
    setpar(0.01, 15.0, -60.0, 10.0);
    seth(0.1);
	double w_[2] = {vr,0.0};
    setw0(w_);        
}    

thalamicell::thalamicell(double aArg, double bArg, double cArg, double dArg) : neuron() {
    numeq = 2;
    fbuf = new double[numeq];
    w = new double[numeq];
    waux = new double[numeq];
    ks = new double * [numeq]; 
    for (int i = 0; i < numeq; i++)
    	ks[i] = new double[4];         
    startKs(); 
    cap = 200.0, k = 1.6, vr = -60.0, vt = -50.0;
    setpar(aArg, bArg, cArg, dArg);
    seth(0.1);
	double w_[2] = {vr,0.0};
    setw0(w_);           
}

void thalamicell::fx(double inj, double time) {
    inj -=calcsyncurrent(time);
    fbuf[0] = (k * (waux[0] - vr) * (waux[0] - vt) - waux[1] + inj)/cap;
    if (w[0] > -65.0)
        fbuf[1] = - a * waux[1];
    else
        fbuf[1] = a * (b * (waux[0] - vr) - waux[1]);
}

void thalamicell::checkPeak(double time) {
    if (w[0] >= 35.0 + 0.1*w[1]) {
        w[0] = c - 0.1*w[1];
        w[1] = w[1] + d;
        sendevent(time);
        backevent();
    }
}

void thalamicell::setpar(double aArg, double bArg, double cArg, double dArg) {
    a = aArg;
    b = bArg;
    c = cArg;
    d = dArg;
}

void thalamicell::setvrest(double vrest){
    vr = vrest;
	double w_[2] = {vrest,0.0};
    setw0(w_);
}

void thalamicell::setvtresh(double vtresh){
    vt = vtresh;
}

void thalamicell::setk(double k_){
    k = k_;
}

void thalamicell::setcap(double cap_){
    cap = cap_;
}
