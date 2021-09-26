#include "alphasynapse.h"

alphasynapse::alphasynapse() {
    k1 = new double[4];
    k2 = new double[4];
    k3 = new double[4];
    k4 = new double[4];
    gmax = 1.0;
    tau = 1.0;
    delay = 5.0;
    Esyn = 50.0;
    spikes = new vector <double>;
    gsyn = 0.0;
    SHORTERM = 0;
    LONGTERM = 0;
    xfac = 1.0;
    Pa = 0.0;
    M = 0.0;
    taux = 150;
    p = 0.6;         
    seth(0.1);
    mu = 0.0;
    startKs();
}

alphasynapse::alphasynapse(double tauArg, double gmaxArg, double delayArg, double EsynArg) {
    k1 = new double[4];
    k2 = new double[4];
    k3 = new double[4];
    k4 = new double[4];
    gmax = gmaxArg;
    tau = tauArg;
    delay = delayArg;
    gsyn = 0.0;
    Esyn = EsynArg;
    spikes = new vector <double>;
    gsyn = 0.0;
    SHORTERM = 0;   
    LONGTERM = 0;
    xfac = 1.0;
    taux = 150;
    p = 0.6;  
    seth(0.1);  
    mu = 0.0; 
    startKs();
}

alphasynapse::alphasynapse(double tauArg, double gmaxArg, double delayArg, double EsynArg, bool shortArg, bool longArg) {
    k1 = new double[4];
    k2 = new double[4];
    k3 = new double[4];
    k4 = new double[4];
    gmax = gmaxArg;
    tau = tauArg;
    delay = delayArg;
    gsyn = 0.0;
    Esyn = EsynArg;
    spikes = new vector <double>;
    gsyn = 0.0;
    SHORTERM = shortArg;  
    LONGTERM = longArg;
    xfac = 1.0;
    taux = 150;
    p = 0.6;         
    seth(0.1);  
    mu = 0.0; 
    startKs();     
}    

void alphasynapse::startKs(){
	for(int i = 0; i < 4; i++){
		k1[i] = 0.0; 	
		k2[i] = 0.0;
		k3[i] = 0.0;
		k4[i] = 0.0;
	}
}

void alphasynapse::setpar(double tauArg, double gmaxArg, double delayArg, double EsynArg) {
    gmax = gmaxArg;
    tau = tauArg;
    delay = delayArg;
    Esyn = EsynArg;
}

void alphasynapse::addevent(double spk) {
    spikes->push_back(spk);
}

void alphasynapse::evaluate(double time) {
    
    if (!LONGTERM && !SHORTERM){
        gsyn += (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6; 
        if (!spikes->empty()){
            double s = time - spikes->at(0) - delay;          
            if (s>0){
                gsyn += gmax;
                spikes->erase(spikes->begin());
            }
        }        
        return;
    } 
    
    if (SHORTERM && !LONGTERM){ 
        gsyn += (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6; 
        xfac += (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6; 
        if (!spikes->empty()){
            double s = time - spikes->at(0) - delay;          
            if (s>0){	        		 
                gsyn += xfac*gmax;
                xfac = p*xfac;
                spikes->erase(spikes->begin());
            }
        }    
        return;
    }
    
    if (LONGTERM && !SHORTERM){ 
        gsyn += (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6; 
        Pa += (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6; 
        M += (k1[3] + 2*k2[3] + 2*k3[3] + k4[3])/6; 
        if (!spikes->empty()){
            double s = time - spikes->at(0) - delay;          
            if (s>0){	        		 
                gsyn += gmax;
//                printf("%lf\t%lf\n",gsyn,gmax);
                Pa += Aplus;
//                gmax += Aplus*pow((1-getGmax()),mu)*M*0.5;
                gmax += M*0.5;	        		
                if(gmax<0.0) gmax=0.0;
                spikes->erase(spikes->begin());
            }
        }            
        return;
    }	
    
    if (LONGTERM && SHORTERM){
        gsyn += (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6; 
        xfac += (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6; 
        Pa += (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6; 
        M += (k1[3] + 2*k2[3] + 2*k3[3] + k4[3])/6;         
        if (!spikes->empty()){
            double s = time - spikes->at(0) - delay;          
            if (s>0){	        		 
                gsyn += xfac*gmax;
                xfac = p*xfac;
                Pa += Aplus;
//                gmax += Aplus*pow((1-getGmax()),mu)*M*0.5;
                gmax += M*0.5;	        		
                if(gmax<0.0) gmax=0.0;
                spikes->erase(spikes->begin());
            }
        }            
        return;
    }    		
}

void alphasynapse::seth(double hArg) {
    h = hArg;
}

double alphasynapse::getGsyn() {
    return gsyn;
}

double alphasynapse::getE() {
    return Esyn;
}

void alphasynapse::setShortTerm(int argShort, double ap, double ataux){
    SHORTERM = argShort;
    p = ap;
    taux = ataux;
}

void alphasynapse::setLongTerm(int argLong_, double Aminus_, double Aplus_, double mu_, double tauMinxus_, double tauPlus_){
    LONGTERM = argLong_;
    Aminus = Aminus_;
    Aplus = Aplus_;
    tauMinxus = tauMinxus_;
    tauPlus = tauPlus_;
    mu = mu_;
}    

double alphasynapse::getGmax(){
    return gmax;	
}

void alphasynapse::refreshM(){
    if (LONGTERM){
        M -= Aminus;
//        gmax += Aminus*pow(getGmax(),mu)*Pa*0.5;	
        gmax += Pa*0.5;
        if(gmax>200.0) gmax=200.0;
    }    	
}   

double alphasynapse::getPa(){
    return Pa;
}

double alphasynapse::getXfac(){
    return xfac;
}

double alphasynapse::getM(){
    return M;
}  

void alphasynapse::calcK1(){
    gsynaux = gsyn;
    k1[0] = - h * (gsynaux/tau);	
    if (!LONGTERM && SHORTERM){
        xfacaux = xfac;
        k1[1] = - h * ((xfacaux-1)/taux);	
        return;
    }
    if (LONGTERM && !SHORTERM){   
        Paaux = Pa;
        Maux = M;
        k1[2] = - h * (Paaux/tauPlus);
        k1[3] = - h * (Maux/tauMinxus);
        return;
    }
    if (LONGTERM && SHORTERM){
        xfacaux = xfac;
        Paaux = Pa;
        Maux = M;
        k1[1] = - h * ((xfacaux-1)/taux);
        k1[2] = - h * (Paaux/tauPlus);
        k1[3] = - h * (Maux/tauMinxus);
        return;
    }
}

void alphasynapse::calcK2(){
    gsynaux = gsyn + k1[0]/2;
    k2[0] = - h * (gsynaux/tau);	
    if (!LONGTERM && SHORTERM){
        xfacaux = xfac + k1[1]/2;
        k2[1] = - h * ((xfacaux-1)/taux);	
        return;
    }
    if (LONGTERM && !SHORTERM){   
        Paaux = Pa + k1[2]/2;
        Maux = M + k1[3]/2;
        k2[2] = - h * (Paaux/tauPlus);
        k2[3] = - h * (Maux/tauMinxus);
        return;
    }
    if (LONGTERM && SHORTERM){
        xfacaux = xfac + k1[1]/2;
        Paaux = Pa + k1[2]/2;
        Maux = M + k1[3]/2;
        k2[1] = - h * ((xfacaux-1)/taux);
        k2[2] = - h * (Paaux/tauPlus);
        k2[3] = - h * (Maux/tauMinxus);
        return;
    }
    
}

void alphasynapse::calcK3(){
    gsynaux = gsyn + k2[0]/2;
    k3[0] = - h * (gsynaux/tau);	
    if (!LONGTERM && SHORTERM){
        xfacaux = xfac + k2[1]/2;
        k3[1] = - h * ((xfacaux-1)/taux);	
        return;
    }
    if (LONGTERM && !SHORTERM){   
        Paaux = Pa + k2[2]/2;
        Maux = M + k2[3]/2;
        k3[2] = - h * (Paaux/tauPlus);
        k3[3] = - h * (Maux/tauMinxus);
        return;
    }
    if (LONGTERM && SHORTERM){
        xfacaux = xfac + k2[1]/2;
        Paaux = Pa + k2[2]/2;
        Maux = M + k2[3]/2;
        k3[1] = - h * ((xfacaux-1)/taux);
        k3[2] = - h * (Paaux/tauPlus);
        k3[3] = - h * (Maux/tauMinxus);
        return;
    }													
}

void alphasynapse::calcK4(){
    gsynaux = gsyn + k3[0];
    k4[0] = - h * (gsynaux/tau);	
    if (!LONGTERM && SHORTERM){
        xfacaux = xfac + k3[1];
        k4[1] = - h * ((xfacaux-1)/taux);	
        return;
    }
    if (LONGTERM && !SHORTERM){   
        Paaux = Pa + k3[2];
        Maux = M + k3[3];
        k4[2] = - h * (Paaux/tauPlus);
        k4[3] = - h * (Maux/tauMinxus);
        return;
    }
    if (LONGTERM && SHORTERM){
        xfacaux = xfac + k3[1];
        Paaux = Pa + k3[2];
        Maux = M + k3[3];
        k4[1] = - h * ((xfacaux-1)/taux);
        k4[2] = - h * (Paaux/tauPlus);
        k4[3] = - h * (Maux/tauMinxus);
        return;
    }													
}

double alphasynapse::getGaux(){
    return gsynaux;
}

