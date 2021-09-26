#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "tools.cpp"
#include "neuron.h"
#include "izhiCom.h"
#include "fscell.h"
#include "lscell.h"
#include "ltscell.h"
#include "thalamicell.h"
#include "alphasynapse.h"
#include <mpi.h>
#include <omp.h>

using namespace std;

/* 

Area simulada 4mm2
Numero de neuronios (real) 312,9 * 10^3
Fator de escala 1:1

tipos de neurônios 
0	nb1
1	p23
2	b23
3	nb23
4	ss4L4
5	ss4L23
6	p4
7	b4
8	nb4
9	p5L23
10	p5L56
11	b5
12	nb5
13	p6L4
14	p6L56
15	b6
16	nb6	
*/

long rand0::idum;

typedef struct {
	int * posSyn;
	unsigned char type;
	int length;	
} preNeuron;

typedef struct {
	int * preSyn;
	unsigned char type;
	short int * delay;
	float * weight;	
	int length;
} posNeuron;

int createNet (int nneuron, float scale, float sizeNet, preNeuron * mtxPre, posNeuron * mtxPos);
void rasterdata(vector<neuron *> * listNeuron, const char * pchar);

int main(int argc, char * argv[]) {

	int nneuron = 312900, i, j, k, numtasks, rank;
	double t = 0.0f, tmax = 100.0f, h = 0.01f;
	int Nstep = (int) (tmax - t)/h + 1;
	double scale = 1.0, sizeNet = 2.0;
	preNeuron * mtxPre;
	posNeuron * mtxPos;
	
	MPI::Status status; 
	
	nneuron = createNet (nneuron, scale, sizeNet, mtxPre, mtxPos);

	MPI::Init(argc,argv);
	numtasks = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();
	omp_set_num_threads(2);

	neuron * cell;
    alphasynapse * syn;
    double ** spkTimeNo, ** spkTime;
    vector<neuron *> * listNeuron = new vector<neuron *>;
   	unsigned char type;
	int nneuronNo = nneuron/numtasks;
	int displs[numtasks], recv_counts[numtasks];
	
	for (i=0; i< numtasks - 1; i++) {
		displs[i] = i * nneuronNo;
		recv_counts[i] = nneuronNo;
	}
	displs[i] = i * nneuronNo;
	recv_counts[i] = nneuronNo + (nneuron % numtasks); 	
	
	if (rank != numtasks-1){
		spkTime = new double*[2];
		spkTimeNo = new double*[2];
		for (i = 0; i < 2; i++){
			spkTimeNo[i] = new double[nneuronNo];
			spkTime[i] = new double[nneuron];
		}
        for(j = 0; j < nneuronNo; j++){
			spkTimeNo[0][j] = -1.0f;
			spkTimeNo[1][j] = -1.0f;
		}
		int count = 0;
		for (i = rank*nneuronNo; i < (rank+1)*nneuronNo; i++){
			if (type == 0){
				cell = new lscell();
			}
			else if (type == 2 || type == 7 || type == 11 || type == 15){
				cell = new fscell();
			}
			else if (type == 3 || type == 8 || type == 13 || type == 16){
				cell = new ltscell();
			}
			else{
				cell = new izhiCom();
			}
			cell->setId(i);
			cell->seth(h);
			cell->setSpkBuffer(spkTimeNo, count);
			for(j = 0; j < mtxPos[i].length; j++){
				syn = new alphasynapse();
				syn->seth(h);
				type = mtxPos[mtxPos[i].preSyn[j]].type;
				if (type == 0 || type == 2 || type == 3 || type == 7 || type == 8 || type == 11 || type == 12 || type == 15 || type == 16)
					syn->setpar(6.0, mtxPos[i].weight[j], mtxPos[i].delay[j], -80.0);
				else
					syn->setpar(5.0, mtxPos[i].weight[j], mtxPos[i].delay[j], 0.0);
				cell->addsyndend(syn, mtxPos[i].preSyn[j]);
			}
			listNeuron->push_back(cell);
			count++;
		}	

		stringstream sID;
		sID << "spkTime_" << rank  << "_np_" << numtasks <<".dat";

		for(i = 0; i < Nstep; i++){
			t += h;
			#pragma omp parallel for 
			for(j = 0; j < nneuronNo; j++){
				(listNeuron->at(j))->evaluate(5.0, t);
			}
			if((i+1)%100 == 0){
		        
		        MPI_Allgatherv(&spkTimeNo[0][0], nneuronNo, MPI::DOUBLE, 
	                         &spkTime[0][0], recv_counts, displs, 
	                         MPI::DOUBLE, MPI_COMM_WORLD); 
	                         
		        MPI_Allgatherv(&spkTimeNo[1][0], nneuronNo, MPI::DOUBLE,
	                         &spkTime[1][0], recv_counts, displs, 
	                         MPI::DOUBLE, MPI_COMM_WORLD); 
				int idx;
				for(j = 0; j < nneuron; j++){
					if(spkTime[0][j] > 0.0){
						for(k = 0; k < mtxPre[j].length; k++){
							idx = mtxPre[j].posSyn[k] - rank*nneuronNo;
							if (idx >= 0 && idx < nneuronNo){
								(listNeuron->at(idx))->recEvent(j,spkTime[0][j]);
								if(spkTime[1][j] > 0.0)
									(listNeuron->at(idx))->recEvent(j,spkTime[1][j]);
							}
							else
								continue;
						}
					}
				}	
				
				for(j = 0; j < nneuronNo; j++){
					spkTimeNo[0][j] = -1.0f;
					spkTimeNo[1][j] = -1.0f;
				}
			}
		}
		sID.str("");
		sID << "rasterRank_" << rank  << "_np_" << numtasks <<".dat";
		rasterdata(listNeuron, (sID.str()).c_str());	
	}
	else{
		int addNeuron = nneuron % numtasks;
        spkTime = new double*[2];
        spkTimeNo = new double*[2];
        for (i = 0; i < 2; i++){
                spkTimeNo[i] = new double[nneuronNo + addNeuron];
                spkTime[i] = new double[nneuron];
        }	
        for(j = 0; j < nneuronNo + addNeuron; j++){
                spkTimeNo[0][j] = -1.0f;
                spkTimeNo[1][j] = -1.0f;
        }	
		int count = 0;
		for (i = rank*nneuronNo; i < (rank+1)*nneuronNo + addNeuron; i++){
			if (type == 0){
				cell = new lscell();
			}
			else if (type == 2 || type == 7 || type == 11 || type == 15){
				cell = new fscell();
			}
			else if (type == 3 || type == 8 || type == 13 || type == 16){
				cell = new ltscell();
			}
			else{
				cell = new izhiCom();
			}
			cell->setId(i);
			cell->seth(h);
			cell->setSpkBuffer(spkTimeNo, count);
			for(j = 0; j < mtxPos[i].length; j++){
				syn = new alphasynapse();
				syn->seth(h);
				type = mtxPos[mtxPos[i].preSyn[j]].type;
				if (type == 0 || type == 2 || type == 3 || type == 7 || type == 8 || type == 11 || type == 12 || type == 15 || type == 16)
					syn->setpar(6.0, mtxPos[i].weight[j], mtxPos[i].delay[j], -80.0);
				else
					syn->setpar(5.0, mtxPos[i].weight[j], mtxPos[i].delay[j], 0.0);					
				cell->addsyndend(syn, mtxPos[i].preSyn[j]);
			}
			listNeuron->push_back(cell);
			count++;
		}	
		
		for(i = 0; i < Nstep; i++){
			t += h;
			#pragma omp parallel for   
			for(j = 0; j < nneuronNo + addNeuron; j++){
				(listNeuron->at(j))->evaluate(5.0, t);
			}
			if((i+1)%100 == 0){
		        
		        MPI_Allgatherv(&spkTimeNo[0][0], nneuronNo, MPI::DOUBLE, 
	                         &spkTime[0][0], recv_counts, displs, 
	                         MPI::DOUBLE, MPI_COMM_WORLD); 
	                         
		        MPI_Allgatherv(&spkTimeNo[1][0], nneuronNo, MPI::DOUBLE, 
	                         &spkTime[1][0], recv_counts, displs, 
	                         MPI::DOUBLE, MPI_COMM_WORLD); 	                         
                         	
				int idx;
				for(j = 0; j < nneuron; j++){
					if(spkTime[0][j] > 0.0){
						for(k = 0; k < mtxPre[j].length; k++){
							idx = mtxPre[j].posSyn[k] - rank*nneuronNo;
							if (idx >= 0 && idx < nneuronNo + addNeuron){
								(listNeuron->at(idx))->recEvent(j,spkTime[0][j]);
								if(spkTime[1][j] > 0.0)
									(listNeuron->at(idx))->recEvent(j,spkTime[1][j]);
							}
							else
								continue;
						}
					}
				}
				
	    		for(j = 0; j < nneuronNo + addNeuron; j++){
	           		spkTimeNo[0][j] = -1.0f;
				spkTimeNo[1][j] = -1.0f;
	    		}
			}
		}
		stringstream sID;
		sID.str("");
		sID << "rasterRank_" << rank  << "_np_" << numtasks <<".dat";
                rasterdata(listNeuron, (sID.str()).c_str());		
	}
	
	MPI::Finalize();

	return 0;
}

int createNet (int nneuron, float scale, float sizeNet, preNeuron * mtxPre, posNeuron * mtxPos) {
    
    long idum = -1089432789; 
    int i, j, k, l, p, count;
    int nL1, nL23, nL4, nL5, nL6;
    double scaleL1, scaleL23, scaleL4, scaleL5, scaleL6, scalePre, scalePos, sort;	
    
    rand0::setSeed(idum);    
    
    //definindo numero de neuronios por camada e escala por camada
    nL1 = (int) round(sqrt(0.015*nneuron));
    scaleL1 = sizeNet/nL1;
    nL23 = (int) round(sqrt(0.337*nneuron));
    scaleL23 = sizeNet/nL23;    
    nL4 = (int) round(sqrt(0.349*nneuron));
    scaleL4 = sizeNet/nL4;    
    nL5 = (int) round(sqrt(0.076*nneuron));
    scaleL5 = sizeNet/nL5;    
    nL6 = (int) round(sqrt(0.223*nneuron));
    scaleL6 = sizeNet/nL6;    
    
    nneuron = nL1*nL1 + nL23*nL23 + nL4*nL4 + nL5*nL5 + nL6*nL6;
	mtxPre = new preNeuron[nneuron];
	mtxPos = new posNeuron[nneuron];  
	
	//carrega arquivos de conexao e arborizacao axonal
    FILE * connFile, * axonFile;
	double connDat[33][21], axonDat[17][6];    
    connFile = fopen("dataConn.dat","r");
    axonFile = fopen("dataAxon.dat","r");    
    for(i = 0; i < 33*21; i++){
    	fscanf(connFile,"lf", &connDat[i/21][i%21]);
    }   
    for(i = 0; i < 17*6; i++){
    	fscanf(axonFile,"lf", &axonDat[i/6][i%6]);
    }       
    fclose(connFile);
	fclose(axonFile);
    
    //Definindo por camada, tipo de cada neurônio
    
    count = 0;
    for (i = 0; i < nL1*nL1; i++) { 
    	mtxPos[count].type = 0;
	    mtxPre[count].length = 0;
	    mtxPos[count].length = 0;    	
    	count++;
    }
    
    for (i = 0; i < nL23*nL23; i++) { 
    	sort = rand0::next();
    	if (sort < 0.78)
	    	mtxPos[count].type = 1;
    	else if (sort >= 0.78 && sort < 0.87)
	    	mtxPos[count].type = 2;    		
    	else
	    	mtxPos[count].type = 3;    		
    	mtxPre[count].length = 0;
	    mtxPos[count].length = 0;	    	
    	count++;  	
    }
    
    for (i = 0; i < nL4*nL4; i++) { 
    	sort = rand0::next();
    	if (sort < 0.27)
	    	mtxPos[count].type = 4;
    	else if (sort >= 0.27 && sort < 0.54)
	    	mtxPos[count].type = 5;    		
    	else if (sort >= 0.54 && sort < 0.81)
	    	mtxPos[count].type = 6;    		
    	else if (sort >= 0.81 && sort < 0.85)
	    	mtxPos[count].type = 7;    		
    	else
	    	mtxPos[count].type = 8;    		
	    mtxPre[count].length = 0;
	    mtxPos[count].length = 0;    	 	
    	count++;     	
    }
    
    for (i = 0; i < nL5*nL5; i++) {
    	sort = rand0::next();
    	if (sort < 0.64)
	    	mtxPos[count].type = 9;
    	else if (sort >= 0.64 && sort < 0.81)
	    	mtxPos[count].type = 10;    		
    	else if (sort >= 0.81 && sort < 0.89)
	    	mtxPos[count].type = 11;    		
    	else
	    	mtxPos[count].type = 12;    		
	    mtxPre[count].length = 0;
	    mtxPos[count].length = 0;    	  	
    	count++;     	
    }
    
    for (i = 0; i < nL6*nL6; i++) {
    	sort = rand0::next();
    	if (sort < 0.62)
	    	mtxPos[count].type = 13;
    	else if (sort >= 0.62 && sort < 0.82)
	    	mtxPos[count].type = 14;	
    	else if (sort >= 0.82 && sort < 0.91)
	    	mtxPos[count].type = 15;	
    	else
	    	mtxPos[count].type = 16;
	    mtxPre[count].length = 0;
	    mtxPos[count].length = 0;
    	count++;     	
    }
    
	/* A info do arquivo de conexao é lida linha a linha p. Para cada linha
	é varre-se todos o neuronios, aqueles que pertencerem ao tipo especificado
	pela linha é analisado a info da linha em questao. Com a leitura é coletado 
	o indice dos neuronios potencialmente pre-sinapticos (analisando o tipo e
	distância (arquivo de arborizacao axonal) ). Após isso é sorteado dessa lista
	a quantidade de conexoes especificadas no arquivo de conexao */    
    
	int xPre, yPre, xPos, yPos, num;
	int initPre, initPos, nlayerPre, nlayerPos, layerPre, layerPos;		
	double re, dist;
	vector<int> listNeuron;
	vector<int> listDelay;
	
	int * _posSyn;
	int * _preSyn;
	short int * _delay;
	float * _weight;
	
	// varredura das linhas do arquivo de conexao
	for (p = 0; p < 33; p++){
		// varredura de todos o neuronios
		for (i = 0; i < nneuron; i++){
			//verifica se neuronio pos-sinaptico é o tipo especificado na linha
			if (mtxPos[i].type == connDat[p][0]){
				// define a camada, escala e offset na lista de neuronios respectivo ao neuronio pos em questao
				layerPos = (int) connDat[p][1];
				if (mtxPos[i].type == 0){
					initPos = 0;
					nlayerPos = nL1;
					scalePos = scaleL1;
				}
				else if (mtxPos[i].type == 1 || mtxPos[i].type == 2 || mtxPos[i].type == 3){
					initPos = nL1*nL1;
					nlayerPos = nL23;
					scalePos = scaleL23;				
				}
				else if (mtxPos[i].type == 4 || mtxPos[i].type == 5 || mtxPos[i].type == 6 || mtxPos[i].type == 7 || mtxPos[i].type == 8){
					initPos = nL1*nL1 + nL23*nL23;
					nlayerPos = nL4;	
					scalePos = scaleL4;				
				}
				else if (mtxPos[i].type == 9 || mtxPos[i].type == 10 || mtxPos[i].type == 11 || mtxPos[i].type == 12){
					initPos = nL1*nL1 + nL23*nL23 + nL4*nL4;
					nlayerPos = nL5;	
					scalePos = scaleL5;				
				}
				else{
					initPos = nL1*nL1 + nL23*nL23 + nL4*nL4 + nL5*nL5;
					nlayerPos = nL6;	
					scalePos = scaleL6;		
				}			
				//define posicao (sem escala) do neuronio pos	
				xPos = (i - initPos)%nlayerPos;
				yPos = (i - initPos)/nlayerPos;	
				//varre lina do arquivo de conexao, coletando os neuronios pre-sinapticos de um determinado tipo no qual axonio atinge pos sinaptico
				for (j = 4; j < 21; j++){
					if (connDat[p][j] > 0){
						//coleta informacao da arborizacao axonal do neuronio pre na camaa informada pelo arquivo de conexao
						re = axonDat[j - 4][layerPos];
						for (k = 0; k < nneuron; k++){		
							// define a camada, escala e offset na lista de neuronios respectivo ao neuronio pos em questao
							if(mtxPos[k].type == j - 4){
								if (mtxPos[k].type == 0){
									initPre = 0;
									nlayerPre = nL1;
									layerPre = 1;
									scalePre = scaleL1;
								}
								else if (mtxPos[k].type == 1 || mtxPos[k].type == 2 || mtxPos[k].type == 3){
									initPre = nL1*nL1;
									nlayerPre = nL23;
									layerPre = 2;	
									scalePre = scaleL23;				
								}
								else if (mtxPos[k].type == 4 || mtxPos[k].type == 5 || mtxPos[k].type == 6 || mtxPos[k].type == 7 || mtxPos[k].type == 8){
									initPre = nL1*nL1 + nL23*nL23;
									nlayerPre = nL4;	
									layerPre = 3;	
									scalePre = scaleL4;			
								}
								else if (mtxPos[k].type == 9 || mtxPos[k].type == 10 || mtxPos[k].type == 11 || mtxPos[k].type == 12){
									initPre = nL1*nL1 + nL23*nL23 + nL4*nL4;
									nlayerPre = nL5;	
									layerPre = 4;	
									scalePre = scaleL5;			
								}
								else{
									initPre = nL1*nL1 + nL23*nL23 + nL4*nL4 + nL5*nL5;
									nlayerPre = nL6;	
									layerPre = 5;	
									scalePre = scaleL6;		
								}	
								//define posicao (sem escala) do neuronio pre				
								xPre = (k - initPre)%nlayerPre;
								yPre = (k - initPre)/nlayerPre;	
								//verifica se neuronio pre atinge pos		
								dist = pow(xPre*scalePre-xPos*scalePos,2) + pow(yPre*scalePre-yPos*scalePos,2);				
								if (dist < re*re) {
									listNeuron.push_back(k);
									dist = sqrt(pow(dist,2) + pow(abs(layerPre - layerPos)*0.25,2));
									listDelay.push_back((int) round(dist/0.1));									
								}
							}
						}
						
						//quanto neuronios pre de um detrminado tipo realizam sinapse neste neuronio pos
						count = (int) connDat[p][3]*(connDat[p][j]/100.0);						
						//coletados neuronios potencialmente pre-sinapticos é realizado sorteio 
						for (k = 0; k < count; k++){
							sort = round(rand0::next()*listNeuron.size());
							num = listNeuron.at((int) sort);
							
							mtxPos[i].length++;		
							mtxPre[num].length++;
							
							_posSyn = new int[mtxPre[num].length];
							_preSyn = new int[mtxPos[i].length];
							_delay = new short int[mtxPos[i].length];
							_weight = new float[mtxPos[i].length];
							
							for (l = 0; l < (mtxPre[num].length - 1); l++){			
								_posSyn[l] = mtxPre[num].posSyn[l];
							}
							
							_posSyn[l] = i;
							
							for (l = 0; l < (mtxPos[i].length - 1); l++){			
								_preSyn[l] = mtxPos[i].preSyn[l];
								_weight[l] = mtxPos[i].weight[l];
								_delay[l] = mtxPos[i].delay[l];
							}
							
							_preSyn[l] = num;
							_weight[l] = 1.0;								
							_delay[l] = listDelay.at((int) sort);
							
							delete[] mtxPre[num].posSyn;
							delete[] mtxPos[i].preSyn;
							delete[] mtxPos[i].weight;
							delete[] mtxPos[i].delay;
							mtxPre[num].posSyn = _posSyn;
							mtxPos[i].preSyn = _preSyn;
							mtxPos[i].weight = _weight;
							mtxPos[i].delay = _delay;					
						}				
					}
					listNeuron.clear();
					listDelay.clear();
				}	
			}
		}	
	}	
	
	return nneuron;		
}  

void rasterdata(vector<neuron *> * listneuron, const char * pchar) {
    ofstream fout, fplot;
    fout.open(pchar);
    int i, j,k, max_size = 0;    
    
    for (i = 0; i < listneuron->size(); i++)
        if (!((listneuron->at(i))->getevents())->empty())
        	if (((listneuron->at(i))->getevents())->size() > max_size)
        		max_size = ((listneuron->at(i))->getevents())->size();    
	
    for (i = 0; i < listneuron->size(); i++) {
        if (!((listneuron->at(i))->getevents())->empty()) {
        	fout << (listneuron->at(i))->getID()  << "\t";
            for (j = 0; j < ((listneuron->at(i))->getevents())->size(); j++) {
                fout << ((listneuron->at(i))->getevents())->at(j) << "\t";
            }
            if (((listneuron->at(i))->getevents())->size() < max_size){
            	k = max_size - ((listneuron->at(i))->getevents())->size();
            	for (j = 0;j < k; j++)
            		fout << 0.0 << "\t";	
            }
            fout << "\n";
        }
    }
    fout.close();
}
