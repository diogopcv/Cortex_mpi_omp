//
//  lscell.cpp
//  cortex
//
//  Created by Diogo Porfirio de Castro Vieira on 24/11/11.
//  Copyright 2011 Universidade de Sao Paulo. All rights reserved.
//

#include "lscell.h"  

lscell::lscell() : neuron() {
    numeq = 3;
    fbuf = new double[numeq];
    w = new double[numeq];
    waux = new double[numeq];
    ks = new double * [numeq]; 
    for (int i = 0; i < numeq; i++)
    	ks[i] = new double[4];         
    startKs();   
    cap = 20.0, k = 0.3, vr = -66.0, vt = -40.0, capd = 100.0; gcoup = 1.2;
    setpar(0.17, 5.0, -45.0, 100.0);
    seth(0.1);
	double w_[3] = {vr,0.0,vr};
    setw0(w_);        
    setInhibitory();
}    

lscell::lscell(double aArg, double bArg, double cArg, double dArg) : neuron() {
    numeq = 3;
    fbuf = new double[numeq];
    w = new double[numeq];
    waux = new double[numeq];
    ks = new double * [numeq]; 
    for (int i = 0; i < numeq; i++)
    	ks[i] = new double[4];         
    startKs();  
    cap = 20.0, k = 0.3, vr = -66.0, vt = -40.0, capd = 100.0; gcoup = 1.2;
    setpar(aArg, bArg, cArg, dArg);
    seth(0.1);
	double w_[3] = {vr,0.0,vr};
    setw0(w_);           
    setInhibitory();
}

void lscell::fx(double inj, double time) {
    inj -= calcsyncurrent(time);
    fbuf[0] = (k * (waux[0] - vr) * (waux[0] - vt)  + gcoup*(waux[2] - waux[0]) - waux[1] + inj)/cap;
    fbuf[1] = a * (b * (waux[0] - vr) - waux[1]);
    fbuf[2] = (waux[0]-waux[2])/capd;
    return;
}

void lscell::checkPeak(double time) {
    if (w[0] >= 30.0) {
        w[0] = c;
        w[1] = w[1] + d;
        sendevent(time);
        backevent();
    }
}

void lscell::setpar(double aArg, double bArg, double cArg, double dArg) {
    a = aArg;
    b = bArg;
    c = cArg;
    d = dArg;
}

void lscell::setvrest(double vrest){
    vr = vrest;
	double w_[3] = {vrest,0.0,vrest};
    setw0(w_);
}

void lscell::setvtresh(double vtresh){
    vt = vtresh;
}

void lscell::setk(double k_){
    k = k_;
}

void lscell::setcap(double cap_){
    cap = cap_;
}

void lscell::setcapd(double cap_){
    capd = cap_;
}

void lscell::setgcoup(double gcoup_) {
    gcoup = gcoup_;
}