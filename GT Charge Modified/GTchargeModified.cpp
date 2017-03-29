#include <bits/stdc++.h>
using namespace std;
#define pdd pair<double,double>
#define pii pair<int,int>
#define vi vector<int>
#define vvi vector<vi >
#define vpii vector<pii>
#define vvpii vector<vpii>
#define vstrategyStruc vector<strategyStruc>
#define vvstrategyStruc vector<vstrategyStruc>
#define pdi pair<double,int>
#define vpdi vector<pdi>
#define vd vector<double>
#define vvd vector<vd >
#define vpack vector<packet>
#define vvpack vector<vpack>
#define pb push_back
#define SQ(x) ((x)*(x))
#define PI 3.14159


#define NETWORK_SIZE 500
#define NETWORK_RADIUS NETWORK_SIZE/2
#define BASE_STATION pdd(0,0)
#define TIME_SLOT 60
#define NO_OF_GAMES 1
#define DEAD_PERCENT 0.2


#define NO_OF_NODES 500
#define NODE_THRESHOLD 1
#define NODE_TRANSMISSION_RANGE 25
#define NODE_ENERGY_CAPACITY 5
#define MIN_NEIGHBOUR_DIS 15

//CHARGERS
#define NO_OF_CHARGERS 16
#define CHARGER_THRESHOLD 150
#define CHARGER_ENERGY_CAPACITY 1000
#define CHARGER_SPEED 1
#define CHARGER_RADIUS ((CHARGER_SPEED)*(TIME_SLOT))

//energy computation parameters
#define Eelec 5e-8  //in nJ/bit
#define efs 1e-11      // in pJ/bit/m2
#define emp 1.3e-15 // in pJ/bit/m4
#define d0 87.0  //in m
#define message_length 10000 //in bits
#define packet_size 10000   //in bits
#define Eda 0  //in nJ/bit

//GAME THEORY PARAMETERS
#define ALPHA 0.3
#define BETA 0.1
#define GAMMA 0.1
#define ETTA 0.001

#define TOTAL_TIME (NO_OF_GAMES*TIME_SLOT)







//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//WS NODE
typedef struct wsNode
{
	pdd location;
	double remainingEnergy;
	double energyDataRecvTransmit; //receiving + transmission
	double energyForBroadcast;

	wsNode(pdd loc,double remEnergy,double enRT,double enB) :
	    location(loc), remainingEnergy(remEnergy), energyDataRecvTransmit(enRT), energyForBroadcast(enB){}

}wsNode;


typedef struct  packet
{
	pdd location;
	double CN;
	pii requestNo;  //requestNo.first=node no requestNo.second=game

}packet;

typedef struct  strategyStruc
{
    int strategy;
    int nodeNo;
	pdd targetLocation;
	//double CN;
	//pii requestNo;

}strategyStruc;

double randomFloat(){
	double random = (NETWORK_SIZE) * ((double)rand() / (double)RAND_MAX );
    return random;
}

double calculateDistance(pdd a,pdd b){
	return sqrt(SQ(a.first-b.first)+SQ(a.second-b.second));
}

int inRange(double x,double y){
    if(x<0 || x>NETWORK_SIZE || y<0 || y>NETWORK_SIZE)return 0;
    return 1;
}


void generateRandomLocation(std::vector<wsNode> &nodes){
	int nodeCount=NO_OF_NODES;
	while(nodeCount){
		int notPresentFlag=1;
		double xCoordinate = randomFloat();
		double yCoordinate = randomFloat();
		for (int i = 0; i < nodes.size(); ++i)
		{
			if(calculateDistance(pdd(xCoordinate,yCoordinate) , nodes[i].location) < MIN_NEIGHBOUR_DIS
				&& calculateDistance(pdd(xCoordinate,yCoordinate) , BASE_STATION) < MIN_NEIGHBOUR_DIS
                &&inRange(xCoordinate,yCoordinate)){
				notPresentFlag=0;
				break;
			}
		}

		if(notPresentFlag){
			wsNode *nodeObject=new wsNode(pdd(xCoordinate,yCoordinate),NODE_ENERGY_CAPACITY,0.0,0.0);
			nodes.pb(*nodeObject);
			//cout<<xCoordinate<<" "<<yCoordinate<<endl;
			nodeCount--;
			//notPresentFlag=1;
		}
	}

	//n+1 th is the BASE STATION
	nodes.pb({BASE_STATION,INT_MAX,0.0,0.0});
}

void computecalculateDistanceBetweenNodes(vector<wsNode> &nodes,
	vector<vector<double> > &interNodecalculateDistanceMatrix){

	for (int i = 0; i < NO_OF_NODES+1; ++i)//+1 for the base station
	{
		for (int j = 0; j < NO_OF_NODES+1; ++j)
		{
			interNodecalculateDistanceMatrix[i][j]=calculateDistance(nodes[i].location,nodes[j].location);
			//cout<<interNodecalculateDistanceMatrix[i][j]<<" ";
		}
		//cout<<endl;
	}

}

void createGraph(vvi &graph,vvd &interNodecalculateDistanceMatrix ){

	for (int i = 0; i < NO_OF_NODES+1; ++i)
	{
		for (int j = 0; j < NO_OF_NODES+1; ++j)
		{
			if(i!=j && interNodecalculateDistanceMatrix[i][j] < NODE_TRANSMISSION_RANGE){
				graph[i].pb(j);
				graph[j].pb(i);
			}
		}
	}
}


struct comparator
 {
   bool operator()(const pdi& l, const pdi& r)
   {
       return l.first > r.first;
   }
 };



void dijkstraShortestPath(vvi &graph,vi &parent,vd &dist,vvd &interNodecalculateDistanceMatrix){
	int src=NO_OF_NODES;
	dist[src]=0;
	parent[src]=-1;


	priority_queue<pdi,vpdi,comparator > heap;
	heap.push(pdi(dist[src],src));
	while(!heap.empty()){
		pdi top=heap.top();
		heap.pop();

		int currNode=top.second;
		for (int i = 0; i < graph[currNode].size(); ++i)
		{
			int j=graph[currNode][i];
			if(dist[j]>dist[currNode]+interNodecalculateDistanceMatrix[currNode][j]){
				parent[j]=currNode;
				dist[j]=dist[currNode]+interNodecalculateDistanceMatrix[currNode][j];
				heap.push(pdi(dist[j],j));
			}
		}
	}

	//for(int i=0;i<NO_OF_NODES+1;i++)cout<<parent[i]<<" ";
	//cout<<endl;

}


bool double_equals(double a, double b, double epsilon = 0.001)
{
    return std::abs(a - b) < epsilon;
}

double transmissionEnergy(double d){
	double ans=0;
	ans=message_length*Eelec;
	if(d<d0){
		ans+=message_length*efs*SQ(d);
	}
	else{
		ans+=message_length*emp*SQ(SQ(d));
	}
	return ans;
}

double receivingEnergy(){
	return message_length*Eelec;
}

void computeReceivingTransmittingEnergy(vi &parent,vector<wsNode> &nodes){
	vd transEnergy(NO_OF_NODES+1,0);
	vd recvEnergy(NO_OF_NODES+1,0);
	for (int i = 0; i < NO_OF_NODES; ++i)
	{
		int j=i;

		while(parent[j]!=-1){
			transEnergy[j]+=transmissionEnergy(calculateDistance(nodes[j].location,nodes[parent[j]].location));
			recvEnergy[parent[j]]+=receivingEnergy();
			j=parent[j];
		}
	}

	for (int i = 0; i < NO_OF_NODES; ++i)
	{
		nodes[i].energyDataRecvTransmit=transEnergy[i]+recvEnergy[i];
	}
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//WS CHARGER
typedef struct charger
{
	pdd location;
	double remainingEnergy;

}charger;

void initialChargersDeployment(vector<charger> &chargers){
	double radius=NETWORK_SIZE/2;
	double theta=(2.0*PI)/NO_OF_CHARGERS;
	for (int i = 0; i < NO_OF_CHARGERS; ++i)
	{
		double xCoordinate=radius + (5*radius/6)*cos(theta*i);
		double yCoordinate=radius + (5*radius/6)*sin(theta*i);
		//cout<<xCoordinate<<" "<<yCoordinate<<endl;
		chargers.pb({pdd(xCoordinate,yCoordinate),CHARGER_ENERGY_CAPACITY});
	}
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//GAME THEORY
//calculate broadcasting energy
void calculateBroadcastingEnergy(vector<wsNode> &nodes,vvd &RE,vvd &interNodecalculateDistanceMatrix,int game){
    for(int i=0;i<NO_OF_NODES;i++){
        nodes[i].energyForBroadcast=0.0;
    }
    for(int i=0;i<NO_OF_NODES;i++){

            if(RE[i][game]< NODE_THRESHOLD){
                    for(int j=0;j<NO_OF_NODES;j++){
                    double d=interNodecalculateDistanceMatrix[i][j];
                    if(i!=j && (d<2*CHARGER_RADIUS || double_equals(d,2*CHARGER_RADIUS))){
                        nodes[j].energyForBroadcast+=0.001*(receivingEnergy()+transmissionEnergy(d)); //taking directly d for transmission?
                    }
                }
            }
    }
}



//calculating remaining energy for each game
void calculateRemainingEnergy(vector<wsNode> &nodes,vvd &RE,vvd &normalizedRE,vi &noOfRequests,
                                vi &deadNodes,vvi &reqGame,vi &flag,vi &deadTime,int game){
	double minEnergy=INT_MAX,maxEnergy=INT_MIN;
	for (int i = 0; i < NO_OF_NODES; ++i)
	{
		RE[i][game]=RE[i][game-1] -
		 			(nodes[i].energyDataRecvTransmit+nodes[i].energyForBroadcast);
        if(RE[i][game]<0 || double_equals(RE[i][game],0.0)){
            if(deadTime[i]==NO_OF_GAMES+1)
                deadTime[i]=game;
            deadNodes[game]++;
            RE[i][game]=0.0;
            continue;
        }

        if(RE[i][game]< NODE_THRESHOLD){
            noOfRequests[game]++;
            if(flag[i]==0){
                reqGame[i].pb(game);
                flag[i]=1;


            }
        }

		minEnergy=min(minEnergy,RE[i][game]);
		maxEnergy=max(maxEnergy,RE[i][game]);

	}
	for (int i = 0; i < NO_OF_NODES; ++i)
	{
	    if(RE[i][game]<0 || double_equals(RE[i][game],0.0)){
            normalizedRE[i][game]=0;
	    }
	    else
            normalizedRE[i][game]=(RE[i][game]-minEnergy)/(maxEnergy - minEnergy);
	}
}

void initializeCnEnImp(vvd &CN,vvd &EN,vvd &IMP){
	for (int i = 0; i < NO_OF_NODES; ++i)
	{
		CN[i][0]=EN[i][0]=IMP[i][0]=0;
	}
}

void computeCnEnImp(vvd &RE,vvd &normalizedRE,vvd &CN,vvd &EN,vvd &IMP,vi &freq,int game){
	double minmCN=INT_MAX,maxmCN=INT_MIN;
	for (int i = 0; i < NO_OF_NODES; ++i)
	{
	    if(RE[i][game]<0 || double_equals(RE[i][game],0.0))
            continue;
		EN[i][game]=1 - normalizedRE[i][game];
		IMP[i][game]=(freq[i]*1.0)/(game);

		CN[i][game] = (CN[i][game-1]*(1+ALPHA) +
						(EN[i][game] - EN[i][game-1])*(BETA) +
						IMP[i][game]*(GAMMA));
		maxmCN=max(maxmCN,CN[i][game]);
		minmCN=min(minmCN,CN[i][game]);
	}

	for (int i = 0; i < NO_OF_NODES; ++i)
	{

	    if(RE[i][game]<0 || double_equals(RE[i][game],0.0))
            continue;
		CN[i][game]=(CN[i][game]-minmCN)/(maxmCN - minmCN);
	}
}




void updateRequestsForChargers(vector<charger> &chargers,vector<wsNode> &nodes,vvd &RE,vvd &CN,
					vvpack &requestsLessThanR,vvpack &requestsGreaterThanR,int game){
	for (int i = 0; i < NO_OF_NODES; ++i)
	{
		 if(RE[i][game]>0 && RE[i][game]<NODE_THRESHOLD)
		 {
		 	packet packetObj={nodes[i].location,CN[i][game], pii(i,game)};
		 	for (int j = 0; j < NO_OF_CHARGERS; ++j)
		 	{
		 		double d=calculateDistance(nodes[i].location,chargers[j].location);

		 		if(d<CHARGER_RADIUS || double_equals(d,CHARGER_RADIUS)){
		 			requestsLessThanR[j].pb(packetObj);
		 			//cout<<"less than R// "<<d<<" "<<j<<"location: "<<chargers[j].location.first<<" "<<chargers[j].location.second<<endl;
		 		}
		 		else if(d>CHARGER_RADIUS && (d<2*CHARGER_RADIUS|| double_equals(d,2*CHARGER_RADIUS))){
		 			requestsGreaterThanR[j].pb(packetObj);
		 			//cout<<"greater than R//"<<d<<" "<<j<<" location: "<<chargers[j].location.first<<" "<<chargers[j].location.second<<endl;
		 		}

		 	}
		 //	cout<<endl;

		 }
	}
}

void calculatePayoffUpdateSelectedStrategy(vector<charger> &chargers,vector<wsNode> &nodes,vvpack &requestsLessThanR,
                                            vvpack &requestsGreaterThanR,vvstrategyStruc &selectedStrategy,vi &freq,
                                            vvd &CN,vvd &RE,vi &noOfServed,vvi &servedGame,vd &energyTransferred,vi &flag,int game){
	for (int i = 0; i < NO_OF_CHARGERS; ++i)
	{
		int strategy=0;
		int PcNODE=-1;
		int PmNODE=-1;
		double PcMAX=-1;
		double PmMAX=-1;

		for (int req = 0; req < requestsLessThanR[i].size(); ++req)
		{
			double d=calculateDistance(chargers[i].location,requestsLessThanR[i][req].location);

			double profit=(requestsLessThanR[i][req].CN*SQ(CHARGER_RADIUS))/SQ(d);
			if(PcMAX < profit)
			{
				PcMAX=profit;
				PcNODE=requestsLessThanR[i][req].requestNo.first;
			}
		}


		for (int req = 0; req < requestsGreaterThanR[i].size(); ++req)
		{
			double d=calculateDistance(chargers[i].location,requestsGreaterThanR[i][req].location);

			double profit=ETTA*(requestsGreaterThanR[i][req].CN*SQ(CHARGER_RADIUS))/SQ(CHARGER_RADIUS-d);
			if(PmMAX < profit)
			{
				PmMAX=profit;
				PmNODE=requestsGreaterThanR[i][req].requestNo.first;
			}
		}


        if(PcNODE==-1 && PmNODE==-1){
            selectedStrategy[i][game].strategy=0;
			selectedStrategy[i][game].targetLocation=selectedStrategy[i][game-1].targetLocation;
            //cout<<"returning to main"<<endl;
            //return;
        }
		else if(PcNODE!=-1 &&(PcMAX > PmMAX || double_equals(PcMAX,PmMAX)) && (RE[PcNODE][game]>0 && RE[PcNODE][game]< NODE_THRESHOLD))
		{
                selectedStrategy[i][game].strategy=1;
                selectedStrategy[i][game].nodeNo =PcNODE;
                selectedStrategy[i][game].targetLocation =nodes[PcNODE].location;
                freq[PcNODE]++;
                CN[PcNODE][game]=0.0;
                energyTransferred[game]+=NODE_ENERGY_CAPACITY-RE[PcNODE][game];
                RE[PcNODE][game]=NODE_ENERGY_CAPACITY;

                noOfServed[game]++;
                servedGame[PcNODE].pb(game);
                flag[PcNODE]=0;

                //cout<<"<R d: "<<calculateDistance(chargers[i].location,nodes[PcNODE].location)<<endl;
                chargers[i].location=nodes[PcNODE].location;

		}
		else if(PmNODE!=-1)
		{
			selectedStrategy[i][game].strategy=2;
			selectedStrategy[i][game].nodeNo =PmNODE;


			double x2=nodes[PmNODE].location.first;
            double y2=nodes[PmNODE].location.second;
            double temp_x1=chargers[i].location.first;
            double temp_y1=chargers[i].location.second;
            double slope;
            slope=(y2-temp_y1)/(x2-temp_x1);  //if slope is infinite? handle the case
            double theta=atan(abs(slope));

            double temp_x,temp_y;
            if(slope>=0){
                if(x2-temp_x1<0)theta=PI+ theta;
            }
            else{
                theta=PI-theta;
                if(x2-temp_x1>0)theta=PI+theta;
            }

            temp_x = temp_x1 + CHARGER_RADIUS*cos(theta);
            temp_y = temp_y1 + CHARGER_RADIUS*sin(theta);
//            if(temp_x<0 || temp_y<0 || temp_x>500 || temp_y>500){
//                    temp_x=temp_x1;
//                    temp_y=temp_y1;
//                    selectedStrategy[i][game].strategy=0;
//            }
            //if(temp_y<0)temp_y=temp_y1;
            //if(temp_x>500)temp_x=temp_x1;
            //if(temp_y>500)temp_y=temp_y1;
            selectedStrategy[i][game].targetLocation = pdd(temp_x,temp_y);
            //cout<<">R d: "<<calculateDistance(chargers[i].location,pdd(temp_x,temp_y))<<endl;
			chargers[i].location=pdd(temp_x,temp_y);
			//cout<<">R d: "<<d<<endl;
		}
		else{

            selectedStrategy[i][game].strategy=0;
			selectedStrategy[i][game].targetLocation=selectedStrategy[i][game-1].targetLocation;

		}
	}

}


void calculateAverageWaitingTime(vvi &reqGame,vvi &servedGame,vd &waitingTime){

    for(int i=0;i<NO_OF_NODES;++i){

        int j;
        for(j=0;j<servedGame[i].size();++j){
            waitingTime[i]=waitingTime[i]+servedGame[i][j]-reqGame[i][j];
        }

        if(reqGame[i].size()>servedGame[i].size())
            waitingTime[i]=waitingTime[i]+NO_OF_GAMES+1-reqGame[i][j];
        if(reqGame[i].size()>0)
            waitingTime[i]/=reqGame[i].size();

    }

}

void initializeSelectedStrategy(vvstrategyStruc &selectedStrategy,vector<charger> &chargers ){
    for(int i=0;i<NO_OF_CHARGERS;i++){
        selectedStrategy[i][0].targetLocation=chargers[i].location;
        selectedStrategy[i][0].strategy=0;
        selectedStrategy[i][0].nodeNo=-1;
    }
}

void calculateDistanceTravelledByCharger(vd &distanceTravelled,vvstrategyStruc &selectedStrategy){
    for(int i=0;i<NO_OF_GAMES;i++){
        double avg=0;
        for(int j=0;j<NO_OF_CHARGERS;j++){
            avg+=calculateDistance(selectedStrategy[j][i].targetLocation,selectedStrategy[j][i+1].targetLocation);
        }
        avg/=NO_OF_CHARGERS;
        distanceTravelled[i]=avg;
    }

}

void calculateTravellingEfficiency(vi &nodesCharged,vd &distanceTravelled,vd &travellingEfficiency){
    for(int i=0;i<NO_OF_GAMES;i++){
        if(double_equals(distanceTravelled[i],0.0))travellingEfficiency[i]=0;
        else
            travellingEfficiency[i]=nodesCharged[i]/distanceTravelled[i];
    }

}

void calculateTravellingEfficiencyCumulative(vi &noOfServed,vd &distanceTravelled,vd &travellingEfficiencyCumulative){
    int totalNodesCharged=0;
    double totalDistTravelled=0;
    for(int i=0;i<NO_OF_GAMES;i++){
        totalDistTravelled+=distanceTravelled[i];
        totalNodesCharged+=noOfServed[i];
        if(double_equals(totalDistTravelled,0.0))travellingEfficiencyCumulative[i]=0;
        else
            travellingEfficiencyCumulative[i]=totalNodesCharged/totalDistTravelled;
    }

}

void calculateAverageEnergyTransferred(vd &energyTransferred,vd &avgEnergyTransferred){
    for(int i=0;i<NO_OF_GAMES;i++){
        avgEnergyTransferred[i]=energyTransferred[i]/NO_OF_CHARGERS;
    }
}

//#####################################################################


//graphs
void calculateNetworkLifetime(vi &deadNodes,double &lifetime){
    for(int i=0;i<=NO_OF_GAMES;i++){
        if(deadNodes[i]>=(int)(DEAD_PERCENT*NO_OF_NODES)){
                lifetime=i*TIME_SLOT;
                break;
        }
    }
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void readNodeLocationFromFile(vector<wsNode> &nodes){

    FILE *fx_coord,*fy_coord;

    fx_coord=fopen("F:\\project8\\node location\\Put_NODE_X.txt","r");
    fy_coord=fopen("F:\\project8\\node location\\Put_NODE_Y.txt","r");

    float xCoordinate,yCoordinate;

    for(int i=0;i<NO_OF_NODES;++i)
    {

        fscanf(fx_coord,"%f ",&xCoordinate);
        fscanf(fy_coord,"%f ",&yCoordinate);

        wsNode *nodeObject=new wsNode(pdd(xCoordinate,yCoordinate),NODE_ENERGY_CAPACITY,0.0,0.0);
        nodes.pb(*nodeObject);
        //cout<<i<<" "<<nodes[i].location.first<<" "<<nodes[i].location.second<<endl;
    }
    nodes.pb({BASE_STATION,INT_MAX,0.0,0.0});

    fclose(fx_coord);
    fclose(fy_coord);

}


void writeNodeLocationToFile(vector<wsNode> &nodes)
{
    FILE *fn,*fx_coord,*fy_coord;
    fn=fopen("F:\\project8\\GT Charge Modified\\Put_NO_OF_NODES.txt","w");
    fprintf(fn,"%d ",NO_OF_NODES);

    fx_coord=fopen("F:\\project8\\GT Charge Modified\\Put_NODE_X.txt","w");
    fy_coord=fopen("F:\\project8\\GT Charge Modified\\Put_NODE_Y.txt","w");


    for(int i=0;i<NO_OF_NODES;++i)
    {

        fprintf(fx_coord,"%f ",nodes[i].location.first);
        fprintf(fy_coord,"%f ",nodes[i].location.second);
    }


    fclose(fn);
    fclose(fx_coord);
    fclose(fy_coord);
}

void writeChargerLocationToFile(vector<charger> &chargers)
{
    FILE *fn,*fx_coord,*fy_coord;
    fn=fopen("F:\\project8\\GT Charge Modified\\Put_NO_OF_CHARGERS.txt","w");
    fprintf(fn,"%d ",NO_OF_CHARGERS);

    fx_coord=fopen("F:\\project8\\GT Charge Modified\\Put_CHARGER_X.txt","w");
    fy_coord=fopen("F:\\project8\\GT Charge Modified\\Put_CHARGER_Y.txt","w");


    for(int i=0;i<NO_OF_CHARGERS;++i)
    {

        fprintf(fx_coord,"%f ",chargers[i].location.first);
        fprintf(fy_coord,"%f ",chargers[i].location.second);
    }


    fclose(fn);
    fclose(fx_coord);
    fclose(fy_coord);
}


void writeRequestZoneLessThanRToFile(vvpack &requestsLessThanR,int game){


    FILE *freq;

    freq=fopen("F:\\project8\\GT Charge Modified\\Put_REQ_LESS_THAN_R.txt","a");



    for(int i=0;i<NO_OF_CHARGERS;++i)
    {
        fprintf(freq,"%d %d :",game,i);
        for(int j=0;j<requestsLessThanR[i].size();++j){
                fprintf(freq,"%d ",requestsLessThanR[i][j].requestNo.first);
        }
        fprintf(freq,"\n");
    }


    fclose(freq);

}


void writeRequestZoneGreaterThanRToFile(vvpack &requestsGreaterThanR,int game){


    FILE *freq;

    freq=fopen("F:\\project8\\GT Charge Modified\\Put_REQ_GREATER_THAN_R.txt","a");



    for(int i=0;i<NO_OF_CHARGERS;++i)
    {
        fprintf(freq,"%d %d :",game,i);
        for(int j=0;j<requestsGreaterThanR[i].size();++j){
                fprintf(freq,"%d ",requestsGreaterThanR[i][j].requestNo.first);
        }
        fprintf(freq,"\n");
    }


    fclose(freq);

}

void writeSelectedStrategyToFile(vvstrategyStruc &selectedStrategy)
{
    FILE *fgames,*fstrategy;
    fgames=fopen("F:\\project8\\GT Charge Modified\\Put_NO_OF_GAMES.txt","w");
    fprintf(fgames,"%d ",NO_OF_GAMES);

    fstrategy=fopen("F:\\project8\\GT Charge Modified\\Put_STRATEGY.txt","w");



    for(int i=0;i<NO_OF_CHARGERS;++i)
    {

            double temp_x,temp_y;

            for(int j=0;j<=NO_OF_GAMES;++j)
            {
                int st=selectedStrategy[i][j].strategy;

                    temp_x=selectedStrategy[i][j].targetLocation.first;
                    temp_y=selectedStrategy[i][j].targetLocation.second;
                    fprintf(fstrategy,"%f %f %d %d\n",temp_x,temp_y,st,i);
                   // cout<<"selected strategy of charger "<<i<<" in game "<<j<<"is "<<st<<" and target location is ";
                    //cout<<"("<<temp_x<<") , ("<<temp_y<<")"<<endl;

            }
    }


    fclose(fgames);
    fclose(fstrategy);

}



void writeRequestsToFile(vi &noOfRequests, vi &noOfServed){


    FILE *freq,*freqServed,*ftemp;
    freq=fopen("F:\\project8\\GT Charge Modified\\Put_REQUESTS.txt","w");
    freqServed=fopen("F:\\project8\\GT Charge Modified\\Put_REQ_SERVED.txt","w");
    ftemp=fopen("F:\\project8\\GT Charge Modified\\Put_REQ_TEMP.txt","w");

    for(int i=1;i<=NO_OF_GAMES;++i)
    {
        fprintf(freq,"%d\n",noOfRequests[i]);
        fprintf(freqServed,"%d\n",noOfServed[i]);
        fprintf(ftemp,"%d %d\n",noOfRequests[i],noOfServed[i]);
    }


    fclose(freq);
    fclose(freqServed);
    fclose(ftemp);
}


void writeRemEnergyToFile(vvd &RE){


    FILE *fre;
    fre=fopen("F:\\project8\\GT Charge Modified\\Put_REMENERGY.txt","w");


    for(int i=0;i<=NO_OF_GAMES;++i)
    {

        for(int j=0;j<NO_OF_NODES;++j){

            fprintf(fre,"%f ",RE[j][i]);
//            if(i==NO_OF_GAMES)
//                cout<<j<<" "<<RE[j][i]<<endl;
        }
        fprintf(fre,"\n");
    }


    fclose(fre);
}



void writeDeadNodesToFile(vi &deadNodes){


    FILE *fdead;
    fdead=fopen("F:\\project8\\GT Charge Modified\\Put_NO_OF_DEAD_NODES.txt","w");

    for(int i=1;i<=NO_OF_GAMES;++i)
    {
        fprintf(fdead,"%d\n",deadNodes[i]);
    }


    fclose(fdead);
}



void writeResponseTimeToFile(vvi &reqGame,vvi &servedGame){

    FILE *freq,*fserved;
    freq=fopen("F:\\project8\\GT Charge Modified\\Put_REQ_GAME.txt","w");
    fserved=fopen("F:\\project8\\GT Charge Modified\\Put_SERVED_GAME.txt","w");

    for(int i=0;i<NO_OF_NODES;++i)
    {
//        cout<<reqGame[i].size()<<" "<<servedGame[i].size()<<endl;
//        fprintf(freq,"%d :",i);
        fprintf(freq,"%d ",reqGame[i].size());
//        cout<<"req: "<<reqGame[i].size()<<" served: ";
        for(int j=0;j<reqGame[i].size();++j){

            fprintf(freq,"%d ",reqGame[i][j]);

        }

        fprintf(freq,"\n");
//        fprintf(fserved,"%d :",i);


        fprintf(fserved,"%d ",servedGame[i].size());
//        cout<<servedGame[i].size()<<endl;

         for(int j=0;j<servedGame[i].size();++j){

            fprintf(fserved,"%d ",servedGame[i][j]);

        }

//        if(servedGame[i].size()<reqGame[i].size())
//            fprintf(fserved,"%d ",NO_OF_GAMES+1);
        fprintf(fserved,"\n");
    }


    fclose(freq);
    fclose(fserved);

}


void writeAverageWaitingTimeToFile(vd &waitingTime,vi &deadTime){


    FILE *fwait,*fdead;
    fwait=fopen("F:\\project8\\GT Charge Modified\\Put_AVERAGE_WAITING.txt","w");
    fdead=fopen("F:\\project8\\GT Charge Modified\\Put_DEAD_TIME.txt","w");

    for(int i=0;i<NO_OF_NODES;++i)
    {
        fprintf(fwait,"%f\n",waitingTime[i]);
        fprintf(fdead,"%d\n",deadTime[i]);
    }


    fclose(fwait);
    fclose(fdead);
}

void writeDistanceTravelledByChargerToFile(vd &distanceTravelled){

    FILE *fDistTravelled;

    fDistTravelled=fopen("F:\\project8\\GT Charge Modified\\Put_DISTANCE_TRAVELLED.txt","w");

    for(int i=0;i<NO_OF_GAMES;++i)
    {
        cout<<distanceTravelled[i]<<endl;
        fprintf(fDistTravelled,"%f ",distanceTravelled[i]);

    }

    fclose(fDistTravelled);


}

void writeTravellingEfficiencyOfChargerToFile(vd &travellingEfficiency){

    FILE *ftravellingEfficiency;

    ftravellingEfficiency=fopen("F:\\project8\\GT Charge Modified\\Put_TRAVELLING_EFFICIENCY.txt","w");

    for(int i=0;i<NO_OF_GAMES;++i)
    {

        fprintf(ftravellingEfficiency,"%f ",travellingEfficiency[i]);

    }

    fclose(ftravellingEfficiency);


}

void writeTravellingEfficiencyOfChargerCumulativeToFile(vd &travellingEfficiencyCumulative){

    FILE *ftravellingEfficiency;

    ftravellingEfficiency=fopen("F:\\project8\\GT Charge Modified\\Put_TRAVELLING_EFFICIENCY_CUMULATIVE.txt","w");

    for(int i=0;i<NO_OF_GAMES;++i)
    {

        fprintf(ftravellingEfficiency,"%f ",travellingEfficiencyCumulative[i]);

    }

    fclose(ftravellingEfficiency);


}

void writeAverageEnergyTransferredToFile(vd &avgEnergyTransferred){
    FILE *favgEnergyTransferred;

    favgEnergyTransferred=fopen("F:\\project8\\GT Charge Modified\\Put_AVG_ENERGY_TRANSFERRED.txt","w");

    for(int i=0;i<NO_OF_GAMES;++i)
    {

        fprintf(favgEnergyTransferred,"%f ",avgEnergyTransferred[i]);

    }

    fclose(favgEnergyTransferred);
}


//################################################################################
//GRAPHS files
void appendNetworkLifetimeVaryingNodesToFile(double &lifetime){
    FILE *fnetLife;

    fnetLife=fopen("F:\\project8\\comparison files\\GT Charge Modified\\Put_NETWORK_LIFETIME_VARYING_NODES.txt","a");

    fprintf(fnetLife,"%d %f\n",NO_OF_NODES,lifetime);

    fclose(fnetLife);
}

void appendNetworkLifetimeVaryingChargersToFile(double &lifetime){

    FILE *fnetLife;

    fnetLife=fopen("F:\\project8\\comparison files\\GT Charge Modified\\Put_NETWORK_LIFETIME_VARYING_CHARGERS.txt","a");

    fprintf(fnetLife,"%d %f\n",NO_OF_CHARGERS,lifetime);

    fclose(fnetLife);
}
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//printing all information computed

void printNodesLocation(vector<wsNode> &nodes){
    for(int i=0;i<NO_OF_NODES;i++){
        cout<<nodes[i].location.first<<" "<<nodes[i].location.second<<endl ;

    }

}

void printGraph(vvi &graph){

    for(int i=0;i<=NO_OF_NODES;++i){

        for(int j=0;j<graph[i].size();++j){
            cout<<graph[i][j]<<" ";
        }
        cout<<endl;
    }
}


void printDijkstra(vi &parent){

    for(int i=0;i<=NO_OF_NODES;++i){
        cout<<"parent["<<i<<"] : "<<parent[i]<<endl;
    }
}

void printRecvTransEnergy(vector<wsNode> &nodes){

    for(int i=0;i<NO_OF_NODES;++i){
        if(nodes[i].energyDataRecvTransmit>(1.0/60))
        cout<<i<<" "<<nodes[i].energyDataRecvTransmit<<endl;
    }
}

void printLessThanRZone(vvpack &requestsLessThanR){

    for(int i=0;i<NO_OF_CHARGERS;++i){

        for(int j=0;j<requestsLessThanR[i].size();++j){
            cout<<i<<" "<<requestsLessThanR[i][j].requestNo.first<<endl;
        }
    }
}


void printGreaterThanRZone(vvpack &requestsGreaterThanR){

    for(int i=0;i<NO_OF_CHARGERS;++i){

        for(int j=0;j<requestsGreaterThanR[i].size();++j){
            cout<<i<<" "<<requestsGreaterThanR[i][j].requestNo.first<<endl;
        }
    }
}

void printChargersinR(vector<wsNode> &nodes,vvd &RE,vector<charger> &chargers,int game){

    for(int i=0;i<NO_OF_NODES;++i){
        cout<<"node "<<i<<": RE "<<RE[i][game]<<" -> ";
        for(int j=0;j<NO_OF_CHARGERS;++j){

            double d=calculateDistance(nodes[i].location,chargers[j].location);
            if(d< CHARGER_RADIUS || double_equals(d,CHARGER_RADIUS))
                cout<<j<<" ";
        }
        cout<<endl;
    }
}

void printall(vvd &RE,vvd &normalizedRE,vvd &CN,vvpack &requestsLessThanR,vvpack &requestsGreaterThanR,int game){

    cout<<"Game no: "<<game<<endl;

    cout<<"printing RE:\n";
    for(int i=0;i<NO_OF_NODES;++i)
    {
        cout<<"node no: "<<i<<" "<<RE[i][game]<<" "<<normalizedRE[i][game]<<" CN "<<CN[i][game]<<endl;
    }


//   cout<<"printing normalized RE:\n";
//   for(int i=0;i<NO_OF_NODES;++i)
//   {
//       cout<<"node no: "<<i<<" "<<normalizedRE[i][game]<<endl;
//   }
//
//
//
//    cout<<"printing CN:\n";
//    for(int i=0;i<NO_OF_NODES;++i)
//    {
//        cout<<"node no: "<<i<<" "<<CN[i][game]<<endl;
//    }
//
    cout<<"printing requestsLessThanR:\n";
    for(int i=0;i<NO_OF_CHARGERS;++i)
    {
        cout<<"charger no: "<<i<<" size: "<<requestsLessThanR[i].size()<<endl;
        for(int j=0;j<requestsLessThanR[i].size();++j){

        cout<<"loc: "<<requestsLessThanR[i][j].location.first<<" "<<requestsLessThanR[i][j].location.second<<" ";
        cout<<"Cn: "<<requestsLessThanR[i][j].CN<<" ";
        cout<<"node no: "<<requestsLessThanR[i][j].requestNo.first<<endl;
        }

    }
//
//
//    cout<<"printing requestsGreaterThanR:\n";
//
//    for(int i=0;i<NO_OF_CHARGERS;++i)
//    {
//
//        cout<<"charger no: "<<i<<"size: "<<requestsGreaterThanR[i].size()<<endl;
//         for(int j=0;j<requestsGreaterThanR[i].size();++j){
//        cout<<"loc: "<<requestsGreaterThanR[i][j].location.first<<" "<<requestsGreaterThanR[i][j].location.second<<" ";
//        cout<<"Cn: "<<requestsGreaterThanR[i][j].CN<<" ";
//        cout<<"node no: "<<requestsGreaterThanR[i][j].requestNo.second<<endl;
//         }
//    }
//
//
//



}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
int main(int argc, char const *argv[])
{
	//##############################################
	//GENERATE RANDOM LOCATION FOR WS NODES
	std::vector<wsNode> nodes;
	srand(time(NULL));
	//generateRandomLocation(nodes);
    readNodeLocationFromFile(nodes);

    //printNodesLocation(nodes);
	//calculateDistance between nodes
	vvd interNodecalculateDistanceMatrix(NO_OF_NODES+1,vd(NO_OF_NODES+1));
	computecalculateDistanceBetweenNodes(nodes,interNodecalculateDistanceMatrix);

	//##############################################
	//DEploy chargers initially
	std::vector<charger> chargers;
	initialChargersDeployment(chargers);


    //writing nodes and chargers location to files
    writeNodeLocationToFile(nodes);
    writeChargerLocationToFile(chargers);

	//##############################################
	//Compute path from each node to base station
	//(using Dijkstra's algo)

	//create graph for packet forwarding
	vvi graph(NO_OF_NODES+1);
	createGraph(graph,interNodecalculateDistanceMatrix);
    //printGraph(graph);

	//single source shortest path
	vector<int> parent(NO_OF_NODES+1);
	vector<double> dist(NO_OF_NODES+1,INT_MAX);
	dijkstraShortestPath(graph,parent,dist,interNodecalculateDistanceMatrix);
    //printDijkstra(parent);

	//calculate data packets forwarded by each node
	computeReceivingTransmittingEnergy(parent,nodes);
	printRecvTransEnergy(nodes);



	//#########################################################
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	//GAME THEORY

	//calculate remaining energy
	vvd RE(NO_OF_NODES,vd(NO_OF_GAMES+1));
	vvd normalizedRE(NO_OF_NODES,vd(NO_OF_GAMES+1));

	for (int i = 0; i < NO_OF_NODES; ++i)
	{
		RE[i][0]=NODE_ENERGY_CAPACITY;
		normalizedRE[i][0]=1;
	}



	//compute Charging necessity(CN[]),
	//Energy Needed(EN[])
	//Importance of a node(IMP[])
	vvd CN(NO_OF_NODES,vd(NO_OF_GAMES+1));
	vvd EN(NO_OF_NODES,vd(NO_OF_GAMES+1));
	vvd IMP(NO_OF_NODES,vd(NO_OF_GAMES+1));
	vi freq(NO_OF_NODES,0);
	initializeCnEnImp(CN,EN,IMP);


////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//game has been started

	//data for final trajectory graph

	vvstrategyStruc selectedStrategy(NO_OF_CHARGERS,vstrategyStruc(NO_OF_GAMES+1));
	initializeSelectedStrategy(selectedStrategy,chargers );

	vi noOfRequests(NO_OF_GAMES+1,0);
    vi noOfServed(NO_OF_GAMES+1,0);
    vi deadNodes(NO_OF_GAMES+1,0);
    vvi reqGame(NO_OF_NODES);
    vvi servedGame(NO_OF_NODES);
    vi flag(NO_OF_NODES,0);
    vd waitingTime(NO_OF_NODES,0);
    vi deadTime(NO_OF_NODES,NO_OF_GAMES+1);
    vd energyTransferred(NO_OF_GAMES+1,0);

	for (int game = 1; game <= NO_OF_GAMES; ++game)
	{

		//calculate remaining energy
		calculateRemainingEnergy(nodes,RE,normalizedRE,noOfRequests,deadNodes,reqGame,flag,deadTime,game);
//        printChargersinR(nodes,RE,chargers,game);

        //calculate broadcasting energy
        calculateBroadcastingEnergy(nodes,RE,interNodecalculateDistanceMatrix,game);

		//compute Charging necessity(CN[]),
		//Energy Needed(EN[])
		//Importance of a node(IMP[])
		computeCnEnImp(RE,normalizedRE,CN,EN,IMP,freq,game);



		//mobile charger point of view
		vvpack requestsLessThanR(NO_OF_CHARGERS);
		vvpack requestsGreaterThanR(NO_OF_CHARGERS);
		//findRequestsForChargers(chargers,nodes,path,requestBuffer,requestsLessThanR,requestsGreaterThanR,game);
		updateRequestsForChargers(chargers,nodes,RE,CN,requestsLessThanR,requestsGreaterThanR,game);
//		if(game<50)
//		{
//		    printLessThanRZone(requestsLessThanR);
//            printGreaterThanRZone(requestsGreaterThanR);
//		}

        //printall(RE,normalizedRE,CN,requestsLessThanR,requestsGreaterThanR,game);

		//###################
		//calculate payoff
		calculatePayoffUpdateSelectedStrategy(chargers,nodes,requestsLessThanR,requestsGreaterThanR,selectedStrategy,
                                            freq,CN,RE,noOfServed,servedGame,energyTransferred,flag,game);

//        writeRequestZoneLessThanRToFile(requestsLessThanR,game);
//        writeRequestZoneGreaterThanRToFile(requestsGreaterThanR,game);

        requestsLessThanR.clear();
        requestsGreaterThanR.clear();

        //cout<<"Game "<<game<<" over"<<endl;


	}

calculateAverageWaitingTime(reqGame,servedGame,waitingTime);

//writeSelectedStrategyToFile(selectedStrategy);
writeRequestsToFile(noOfRequests,noOfServed);
//writeRemEnergyToFile(RE);
writeDeadNodesToFile(deadNodes);
writeResponseTimeToFile(reqGame,servedGame);
writeAverageWaitingTimeToFile(waitingTime,deadTime);

//calculate distance travelled
vd distanceTravelled(NO_OF_GAMES,0);
calculateDistanceTravelledByCharger(distanceTravelled,selectedStrategy);

//write distance travelled to file
writeDistanceTravelledByChargerToFile(distanceTravelled);

//travelling efficiency
vd travellingEfficiency(NO_OF_GAMES,0);
calculateTravellingEfficiency(noOfServed,distanceTravelled,travellingEfficiency);
writeTravellingEfficiencyOfChargerToFile(travellingEfficiency);

//cumulative travelling efficiency
vd travellingEfficiencyCumulative(NO_OF_GAMES,0);
calculateTravellingEfficiencyCumulative(noOfServed,distanceTravelled,travellingEfficiencyCumulative);
writeTravellingEfficiencyOfChargerCumulativeToFile(travellingEfficiencyCumulative);

vd avgEnergyTransferred(NO_OF_GAMES+1);
calculateAverageEnergyTransferred(energyTransferred,avgEnergyTransferred);
writeAverageEnergyTransferredToFile(avgEnergyTransferred);

//graphs
//NETWORK LIFETIME
double lifetime=TOTAL_TIME+1;
calculateNetworkLifetime(deadNodes,lifetime);
appendNetworkLifetimeVaryingNodesToFile(lifetime);
//appendNetworkLifetimeVaryingChargersToFile(lifetime);


return 0;
}
