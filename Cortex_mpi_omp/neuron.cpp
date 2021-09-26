/* 
 * File:   neuron.cpp
 * Author: diogopcv
 * 
 * Created on 4 de Fevereiro de 2010, 13:22
 */

#include "neuron.h"
namespace Cells {
	
    neuron::neuron() {
        _ID = 0;
        seth(0.01);
        sdend = new map<int,alphasynapse *>;
        events = new vector <double>;
        setGroup(0);
    }    
    
    neuron::~neuron(){
    	delete[] fbuf;
    	delete[] w;
    	delete[] waux;
    	delete[] ks;
    	delete sdend;
    	delete events; 	    	    	
    }     
	
	void neuron::setw0(double * w_) {
	    for (int i = 0; i < numeq; i++) {
	        w[i] = w_[i];
	    }
	}
	
    void neuron::seth(double hArg) {
        h = hArg;
    }
	
    double neuron::getw(int ind) {
        return w[ind - 1];
    }
	
	double neuron::calcsyncurrent(double time) {
	    double isyn = 0.0;
        map<int,alphasynapse *>::iterator it;
        for(it = sdend->begin(); it != sdend->end(); ++it){    
			isyn += (*it).second->getGaux()*(waux[0] - (*it).second->getE());
        }
	    return isyn;
	}    
	
    void neuron::evaluate(double inj, double time) {
		
		// Calculando k1 
		for (int i = 0; i < numeq; i++)
			waux[i] = w[i];
		for(int i = 0; i < sdend->size(); i++)	
			(sdend->at(i))->calcK1();			            
		fx(inj, time);       
		for (int i = 0; i < numeq; i++)
			ks[i][0] = h * fbuf[i];
		
		// Calculando k2
		for (int i = 0; i < numeq; i++)
			waux[i] = w[i] + ks[i][0]/2;
		for(int i = 0; i < sdend->size(); i++)	
			(sdend->at(i))->calcK2();	              
		fx(inj, time);
		for (int i = 0; i < numeq; i++)
			ks[i][1] = h * fbuf[i];
		
		// Calculando k3
		for (int i = 0; i < numeq; i++)
			waux[i] = w[i] + ks[i][1]/2;
		for(int i = 0; i < sdend->size(); i++)	
			(sdend->at(i))->calcK3();  	             
		fx(inj, time);
		for (int i = 0; i < numeq; i++)
			ks[i][2] = h * fbuf[i];
		
		// Calculando k4
		for (int i = 0; i < numeq; i++)
			waux[i] = w[i] + ks[i][2];
		for(int i = 0; i < sdend->size(); i++)	
			(sdend->at(i))->calcK4(); 	        
		fx(inj, time);
		for (int i = 0; i < numeq; i++)
			ks[i][3] = h * fbuf[i];
		
		// Calculando w's novos
		for (int i = 0; i < numeq; i++)
			w[i] += (ks[i][0] + 2*ks[i][1] + 2*ks[i][2] + ks[i][3])/6;                                              
		
		for(int i = 0; i < sdend->size(); i++)	
			(sdend->at(i))->evaluate(time);             
		
		checkPeak(time);
    }
	
    void neuron::addsyndend(alphasynapse * syn, int idPre) {
        sdend->insert(pair<int, alphasynapse *>(idPre, syn));
    }
    
    void neuron::sendevent(double time) {
        events->push_back(time);
        if (spkBuffer[0][preSpk] < 0.0)
        	spkBuffer[0][preSpk] = time;
        else
      		spkBuffer[1][preSpk] = time;        
    } 
	
    vector<double> *  neuron::getevents(){
        return events;
    }
    
    void neuron::setGroup(int group_){
    	group = group_;
    }
 
    int neuron::getGroupID(){
    	return group;
    }    
    
	void neuron::setExcitatory(){
		typesyn = 1;
	}
	
	void neuron::setInhibitory(){
		typesyn = 0;
	}
	
	int neuron::getTypeSyn(){
		return typesyn;
	}  
	
    void neuron::setId(int myId){
    	_ID = myId;
    }

    int neuron::getID(){
    	return _ID;
    }	  
	
	void neuron::backevent() {
	    int i = 0;
	    while (i < sdend->size()) {
	        (sdend->at(i))->refreshM();
	        i++;
	    }
	}	
	
    void neuron::startKs(){
    	for(int i = 0; i < numeq; i++)
	    	for(int j = 0; j < 4; j++)
    			ks[i][j] = 0.0; 		
    }
    
	void neuron::setSpkBuffer(double ** buffer, int preSpk){
		this->spkBuffer	= buffer;
		this->preSpk = preSpk;
	}   
	
    void neuron::recEvent(int idPre, double spkTime){
    	((*sdend)[idPre])->addevent(spkTime);
    }	   	
}

