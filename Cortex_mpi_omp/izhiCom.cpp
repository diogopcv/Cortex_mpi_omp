/* 
 * File:   izhiCom.cpp
 * Author: diogopcv
 * 
 * Created on 4 de Fevereiro de 2010, 13:22
 */

#include "izhiCom.h"
namespace Cells {
	
	izhiCom::izhiCom() : neuron() {
		numeq = 2;
	    fbuf = new double[numeq];
	    w = new double[numeq];
	    waux = new double[numeq];
        ks = new double * [numeq]; 
        for (int i = 0; i < numeq; i++)
        	ks[i] = new double[4];         
        startKs();
	    cap = 100.0, k = 0.7, vr = -60.0, vt = -40.0;
	    setpar(0.03, -2.0, -50.0, 100.0);
		double w_[2] = {vr,0.0};
	    setw0(w_);   
	    setExcitatory();        
    }
	
    izhiCom::izhiCom(double aArg, double bArg, double cArg, double dArg) : neuron(){
		numeq = 2;
	    fbuf = new double[numeq];
	    w = new double[numeq];
	    waux = new double[numeq];
        ks = new double * [numeq]; 
        for (int i = 0; i < numeq; i++)
        	ks[i] = new double[4];         
        startKs();
	    cap = 100.0, k = 0.7, vr = -60.0, vt = -40.0;
	    setpar(aArg, bArg, cArg, dArg);
		double w_[2] = {vr,0.0};
	    setw0(w_);   
	    setExcitatory();     
    }        
	
    void izhiCom::setpar(double aArg, double bArg, double cArg, double dArg) {
        a = aArg;
        b = bArg;
        c = cArg;
        d = dArg;
    }
	
    void izhiCom::fx(double inj, double time) {
	    inj -= calcsyncurrent(time);
	    fbuf[0] = (k * (waux[0] - vr) * (waux[0] - vt) - waux[1] + inj)/cap;
	    fbuf[1] = a * (b * (waux[0] - vr) - waux[1]);
    }
    
	void izhiCom::checkPeak(double time) {
	    if (w[0] >= 30.0) {
	        w[0] = c;
	        w[1] = w[1] + d;
	        sendevent(time);
	        backevent();
	    }
	}    
    
	void izhiCom::setvrest(double vrest){
	    vr = vrest;
		double w_[2] = {vrest,0.0};
	    setw0(w_);
	}
	
	void izhiCom::setvtresh(double vtresh){
	    vt = vtresh;
	}
	
	void izhiCom::setk(double k_){
	    k = k_;
	}
	
	void izhiCom::setcap(double cap_){
	    cap = cap_;
	}
	
	void izhiCom::turnRS(){
	    setpar(0.03, -2.0, -50.0, 100.0);
	    setk(0.7);
	    setvrest(-60.0);
	    setvtresh(-40.0);
	    setcap(100.0);
	    
	}
	
	void izhiCom::turnIBS(){
	    setpar(0.01, 5.0, -56.0, 130.0);
	    setk(1.2);
	    setvrest(-75.0);
	    setvtresh(-45.0);
	    setcap(150.0);
	}
	
	void izhiCom::turnCHS()   {
	    setpar(0.03, 1.0, -40.0, 150.0);
	    setk(1.5);
	    setvrest(-60.0);
	    setvtresh(-40.0);
	    setcap(50.0);
	    
	}    
}

