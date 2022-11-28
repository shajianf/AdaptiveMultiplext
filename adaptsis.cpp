//adaptive networks for sis model. we use Gross paper's values and we have good agreement between them.25.3.1401
//we try to use adaptive network in my multiplex networks.4/5/1401

#include<iostream>
#include<time.h>
#include<stdlib.h>
#include<fstream>
#include<math.h>
#include<cstdlib>
#include<vector>

#define N 100000

using namespace std;

ofstream Time;
ofstream ER;

//int adj_arr[1000][1000];          // Adjecancy matrix
vector <vector<int>> A(N);          // Adjecancy matrix
vector <vector<int>> B(N);

double p = 0.0001;                   // p= k/ n-1;
int N0 = 10000;                       // number of infected node at time 0.

int ensemble = 1;
int tmax = 4000;                     // time steps.

double beta = 0.0;                    // probability of a node get infected
double muo = 0.002;                    // probability of recovery

double gama = 0.98;
double kapa = 0.5;

double rew = 0.04;

double phi = 0.2;

int sis[N]={0};                     //state of each node ; 1= infected , 0= suseptible.at time t. 
int sis_up[N]={0};                  //state of each node at time t+1.
int thresh[N]={0};
int thresh_up[N]={0};


int main(){

    clock_t runtime = clock();

    Time.open("Inf-forw rew0.04.txt");
    ER.open("actf.04.txt");

    srand(time(NULL));

    /*for(int b =0; b < N; b++){          //reset

        A[b].clear();
                //sis [b]=0;
                //sis_up [b]=0;

    }
*/

    for (int i=0; i< N;i++){            //generate ER network
        for (int j=0; j<i; j++){

            float ran = rand()/(1.+RAND_MAX);
            if( ran< p){

                A[i].push_back(j);
                A[j].push_back(i); 
                        // A[i][j]++;
                       // A[j][i]++; 
            }
            float ran1 = rand()/(1.+RAND_MAX);
            if( ran1 < p){

                B[i].push_back(j);
                B[j].push_back(i);
            }
        }
    }
    
    double var =0;
    double beta_A=0;

    while ( beta < 0.02 ){
        //cout<<"++++++++++++++++++"<< " " << beta <<endl;
        beta_A = gama * beta;

        int ens[ensemble] ={0};
        int ens_A[ensemble] ={0};
     
        for(int e=0 ; e <ensemble ; e++){  
            
            if(var == 0){
                int y=0;
                while(y<N0 ){                       //initialization
                    int O=0;
                    O = rand()%N;
                    if (sis_up[O]==0){
                        sis[O]=1;
                        sis_up[O]=1;
                        y++;
                    }
                }
            }
            for (int c1=0; c1<N; c1++){
                thresh[c1]=0;
                thresh_up[c1]=0;
            }

            int y=0;
            while(y<N0 ){                       //initialization
                int O=0;
                O = rand()%N;
                if (thresh_up[O]==0){
                    thresh[O]=1;
                    thresh_up[O]=1;
                    y++;
                }
            }
            //cout << "aaa" << y <<endl;

           // double control1 = 10000;                //time window
            //double control2 = 0;

           // int width =40;

            int t=0;  
            while(t<tmax){

                for (int f=0; f<N; f++){
                    int act_nei =0;
                    for (int f0=0; f0 < B[f].size();f0++){
                        
                        int q0=0;
                        q0 = B[f][f0];

                        if (thresh[q0]){

                            act_nei++;
                           // cout << ">>>" <<act_nei <<endl;
                        }

                    }

                    int k =0;
                    k = B[f].size();
                    //cout << phi <<'\t' << ((double)act_nei/k) <<endl;
                    if(phi < ((double) act_nei/k)){

                        thresh_up[f]=1;
                    }else{

                        thresh_up[f]=0;
                    }


                    float ran1= rand()/(1.+RAND_MAX);

                    if ((sis[f]) && (ran1 < muo)){              //recovery with mu prob.
                        //cout <<"hello" << endl;
                        sis_up[f]=0;
                    }  

                    for(int h=0;h< A[f].size(); h++){

                        bool flag =false;
                        int q=0;
                        q = A[f][h];

                        float ran2 = rand()/(1.+RAND_MAX);

                        if ((thresh[f]==1)&&(sis[f]==0)&&(ran2<(beta_A *sis[q]))){

                            sis_up[f]=1;
                            flag=true;
                        }
                        if (flag) break;

                        if((thresh[f]==0)&&(sis[f]==0) && ran2 < (beta * sis[q])){          //infect with beta prob.

                            sis_up[f]=1;
                            flag = true;

                        }  
                        if(flag) break; 
                    }

                    double ranR= rand()/(1.+RAND_MAX);
                    if((sis[f]==1)&&(thresh[f]==0)&& kapa<ranR){

                        thresh_up[f]=1;
                    }


                }//for f. 


                /*for (int c=0;c < N; c++){
                    
                    Inf_nei[c].clear();
                } 

                Inf_nei.resize(N);*/      
                
                
                //var =0;
                vector <int> sus_to(N);
                sus_to.clear();

                for( int u=0; u<N ; u++){                               //change step.
                    
                    sis[u] = sis_up[u];
                    thresh[u] = thresh_up[u];
                    //var += sis_up[u];

                    if ((sis[u] ==0)&&(thresh[u]==1)){

                        sus_to.push_back(u);
                    }


                }


                int oo =0;
                for (int c1=0;c1 <N;c1++){
                    
                    
                    for (int c2=0; c2 < A[c1].size(); c2++){

                        int q0 = 0;
                        q0 = A[c1][c2];

                        double ran3 = 0;
                        ran3 = rand()/(1.+RAND_MAX);

                        if ((sis[q0]) && (sis[c1]==0)&&(thresh[c1]==1)&& (ran3< rew)){

                            //bool flag = true;
                            int counter =0;
                            for (int s=0 ; s< sus_to.size() ;s++){

                                double R = rand()%(sus_to.size());
                                counter ++;
                                int q3 =0;

                                q3 = sus_to[R];

                                bool multi_edge = true;

                                for ( int v1=0; v1 < A[c1].size();v1++){

                                    if (A[c1][v1] ==q3){

                                        multi_edge = false;

                                        break;
                                    }
                                }

                                if (q3 != c1 && multi_edge){
                                
                                    //flag = false;

                                    A[c1].erase(A[c1].begin()+c2);

                                    for (int v=0;v < A[q0].size();v++){
                                
                                        if (A[q0][v] == c1){

                                            A[q0].erase(A[q0].begin()+v);
                                        }

                                    }

                                    A[c1].push_back(q3);
                                    A[q3].push_back(c1);

                                   // ER << c1 << '\t' << q3 <<"\t" << A[c1].size()<<endl;

                                    //oo ++;
                                    break;
                                }

                               // if (counter == sus_to.size()) flag= false;
                            }//for sus_total
                            //oo += Inf_nei[c1][q0];

                        }
                    }//for c2

                    

                    
                }//for c1.

                //ER << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" <<endl;

               // cout << oo << "\t" << var <<endl;

               // cout << var  << endl;
                
               /* control2 += var;

                if(t% width ==0){

                    control2 /= (width);
                   // cout << control1 <<"+++++++++++" << control2  << ">>>>>>" << var<< endl;
                    if(abs(control1-control2)< .001){

                        control2 =0;
                        t = tmax;
                    }else{

                        control1 = control2;
                        control2 =0;
                    }
                }*/
                

                t++;
              //  cout << t <<endl;
            }//end of while t.
          //  cout << "????????????????????????????????????????????" << t <<endl;
            var =0;
            int densityofinf =0;
            int densityofact =0;
            for(int b=0; b < N; b++){

                densityofinf += sis_up[b];
                densityofact += thresh_up[b];
                var += sis_up[b];
            }

            ens[e] = densityofinf;
            ens_A[e] = densityofact;
        }//end of ensemble.

        double inf =0;
        double act =0;
        for(int a=0; a<ensemble ;a++){
            inf +=ens[a];
            act +=ens_A[a];

            ens[a]=0;
            ens_A[a]=0;
        }

        Time << beta << '\t' << inf/double(N * ensemble )<< endl;
        ER << beta << '\t' << act/double(N * ensemble )<< endl;

        cout << beta <<endl;

        beta += 0.0001;
    }

    cout << "run time = " << (double) (clock() - runtime) / (CLOCKS_PER_SEC * 60) << " min" << endl;

 return 0;

}