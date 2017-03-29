#include <bits/stdc++.h>
using namespace std;
#define pdd pair<double,double>
#define pii pair<int,int>
#define vi vector<int>
#define vvi vector<vi>
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
#define NO_OF_GAMES 1000
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

#define ALPHA 0.3
#define BETA 0.667
#define DEL_B -0.005
#define DEL_L -0.005
#define DEL_U_NODE 4
#define DEL_U_CHARGER 100
#define DEL_R_NODE -1
#define DEL_R_CHARGER -4

#define TOTAL_TIME (NO_OF_GAMES*TIME_SLOT)

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//WS NODE
typedef struct wsNode
{
	pdd location;
	double remainingEnergy;
	//
	double energyDataRecvTransmit; //receiving + transmission
	double energyForBroadcast;
    double contribution;
	double priority;

	wsNode(pdd loc,double remEnergy,double enRT,double enB,double contri,double priori) :
	    location(loc), remainingEnergy(remEnergy), energyDataRecvTransmit(enRT), energyForBroadcast(enB),
	    contribution(contri), priority(priori){}

}wsNode;

typedef struct  strategyStruc
{
    int strategy;
	pdd targetLocation;
	//double CN;
	//pii requestNo;

}strategyStruc;


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//WS CHARGER
typedef struct charger
{
	pdd location;
	double remainingEnergy;
	double contribution;
	double priority;

}charger;

void initialChargersDeployment(vector<charger> &chargers){
//	for (int i = 0; i < NO_OF_CHARGERS; ++i)
//	{
//		chargers.pb({BASE_STATION,CHARGER_ENERGY_CAPACITY,0.0,0.0});
//	}

    double radius=NETWORK_SIZE/2;
	double theta=(2.0*PI)/NO_OF_CHARGERS;
	for (int i = 0; i < NO_OF_CHARGERS; ++i)
	{
		double xCoordinate=radius + (5*radius/6)*cos(theta*i);
		double yCoordinate=radius + (5*radius/6)*sin(theta*i);
		//coutcout<<xCoordinate<<" "<<yCoordinate<<endl;
		chargers.pb({pdd(xCoordinate,yCoordinate),CHARGER_ENERGY_CAPACITY});
	}
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//

bool double_equals(double a, double b, double epsilon = 0.001)
{
    return std::abs(a - b) < epsilon;
}

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
			wsNode *nodeObject=new wsNode(pdd(xCoordinate,yCoordinate),NODE_ENERGY_CAPACITY,0.0,0.0,0.0,0.0);
			nodes.pb(*nodeObject);
			//cout<<xCoordinate<<" "<<yCoordinate<<endl;
			nodeCount--;
			//notPresentFlag=1;
		}
	}

	//n+1 th is the BASE STATION
	nodes.pb({BASE_STATION,INT_MAX,0.0,0.0,0.0,0.0});
}

void computecalculateDistanceBetweenNodes(vector<wsNode> &nodes,
	vector<vector<double> > &interNodecalculateDistanceMatrix){

	for (int i = 0; i < NO_OF_NODES+1; ++i)//+1 for the base station
	{
		for (int j = 0; j < NO_OF_NODES+1; ++j)
		{
			interNodecalculateDistanceMatrix[i][j]=calculateDistance(nodes[i].location,nodes[j].location);
			//cout<<interNodecalculateDistanceMatrix[i][j]<<endl;
		}
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

void initializeSelectedStrategy(vvstrategyStruc &selectedStrategy,vector<charger> &chargers ){
    for(int i=0;i<NO_OF_CHARGERS;i++){
        selectedStrategy[i][0].targetLocation=chargers[i].location;
        selectedStrategy[i][0].strategy=0;
    }
}

//calculating remaining energy for each game
void calculateRemainingEnergyOfNode(vector<wsNode> &nodes,vvd &RENode,vvd &normalizedRENode,vi &noOfRequests,
                                    vi &deadNodes,vvi &reqGame,vi &flag,vi &deadTime,int game){
	double minEnergy=INT_MAX,maxEnergy=INT_MIN;
	for (int i = 0; i < NO_OF_NODES; ++i)
	{
		RENode[i][game]=RENode[i][game-1] -
		 			(nodes[i].energyDataRecvTransmit+nodes[i].energyForBroadcast);

        if(RENode[i][game]<0 || double_equals(RENode[i][game],0.0)){
            if(deadTime[i]==NO_OF_GAMES+1)
                deadTime[i]=game;
            RENode[i][game]=0.0;
            deadNodes[game]++;
            continue;
        }
        if(RENode[i][game] < NODE_THRESHOLD){
            noOfRequests[game]++;
            if(flag[i]==0){
                reqGame[i].pb(game);
                //cout<<"node "<<i<<" req game: "<<game<<endl;
                flag[i]=1;
            }
        }

		minEnergy=min(minEnergy,RENode[i][game]);
		maxEnergy=max(maxEnergy,RENode[i][game]);
//		if(RENode[i][game]<0 && game<10)
//		cout<<"game "<<game<<" RE["<<i<<"] : "<<RENode[i][game]<<endl;
	}
//    if(game<60)
//	cout<<"game "<<game<<" min: "<<minEnergy<<" max: "<<maxEnergy<<endl;
	for (int i = 0; i < NO_OF_NODES; ++i)
	{
	    if(RENode[i][game]<0 || double_equals(RENode[i][game],0.0)){
                normalizedRENode[i][game]=0.0;
	    }
	    else
            normalizedRENode[i][game]=(RENode[i][game]-minEnergy)/(maxEnergy - minEnergy);
	}
}


void calculateNormalizedREOfCharger(vector<charger> &chargers,vvd &RECharger,vvd &normalizedRECharger,int game){
	double minEnergy=INT_MAX,maxEnergy=INT_MIN;
	for (int i = 0; i < NO_OF_CHARGERS; ++i)
	{
	    RECharger[i][game]=RECharger[i][game-1];

	    if(RECharger[i][game]<0 || double_equals(RECharger[i][game],0.0)){
                RECharger[i][game]=0.0;
                continue;
	    }
		minEnergy=min(minEnergy,RECharger[i][game]);
		maxEnergy=max(maxEnergy,RECharger[i][game]);
	}

	if(minEnergy!=maxEnergy){
        for (int i = 0; i < NO_OF_CHARGERS; ++i)
        {
             if(RECharger[i][game]<0 || double_equals(RECharger[i][game],0.0)){
                        normalizedRECharger[i][game]=0.0;
             }
             else
                    normalizedRECharger[i][game]=(RECharger[i][game]-minEnergy)/(maxEnergy - minEnergy);
        }
	}
	else{
        for (int i = 0; i < NO_OF_CHARGERS; ++i)
        {
            normalizedRECharger[i][game]=RECharger[i][game]/CHARGER_ENERGY_CAPACITY;
        }
	}
}


void computeContributionOfNode(vector<wsNode> &nodes,vvd &contributionNode,
                                vvd &RENode,vvd &normalizedRENode,vvd &profitNode,int game){
    int x=1;//for charged x=1
    for(int i=0;i<NO_OF_NODES;i++){
        if(RENode[i][game] < NODE_THRESHOLD)
            x=0;

        contributionNode[i][game]=(0.5)*(2-x)*(1-ALPHA)*(1-normalizedRENode[i][game]) + ALPHA*(profitNode[i][game-1]);
//        if(isnan(contributionNode[i][game]) && game<350)
//            cout<<"game "<<game<<" node "<<i<<" normalizedRENode "<<normalizedRENode[i][game]<<" profit "<<profitNode[i][game-1]<<endl;
    }
}

void computePriorityOfNode(vector<wsNode> &nodes,vvd &contributionNode,
                           vvd &priorityNode,vvd &RENode,vvd &normalizedRENode,int game){
    for(int i=0;i<NO_OF_NODES;i++){
        if(RENode[i][game]<NODE_THRESHOLD){
            double temp=pow(contributionNode[i][game],BETA);
            priorityNode[i][game]=temp/(1+temp);
        }
        else priorityNode[i][game]=0;
//        if(isnan(priorityNode[i][game]) && game<350)
//            cout<<"game "<<game<<" node "<<i<<"CN "<<contributionNode[i][game]<<endl;
    }
}


void computeContributionOfCharger(vector<charger> &chargers,vvd &contributionCharger,
                                vvd &RECharger,vvd &normalizedRECharger,vvd &profitCharger,int game){

    for(int i=0;i<NO_OF_CHARGERS;i++){
        contributionCharger[i][game]=(1-ALPHA)*contributionCharger[i][game-1]-(ALPHA*(profitCharger[i][game-1]));
    }
}

void computePriorityOfCharger(vector<charger> &chargers,vvd &contributionCharger,
                              vvd &priorityCharger,vvd &RECharger,vvd &normalizedRECharger,int game){
    for(int i=0;i<NO_OF_CHARGERS;i++){
        if(RECharger[i][game] < CHARGER_THRESHOLD){
            double temp=pow(contributionCharger[i][game],BETA);
            priorityCharger[i][game]=temp/(1+temp);
        }
        else priorityCharger[i][game]=0;
    }
}


void createRequestZone(vector<wsNode> &nodes,vector<charger> &chargers,vvd &RENode,
                       vvd &RECharger,vvi &requestZoneNode,vvi &requestZoneCharger,int game){
    for(int i=0;i<NO_OF_CHARGERS;i++){
        //cout<<game<<" "<<i<<" "<<requestZoneNode[i].size()<<" "<<requestZoneCharger[i].size()<<endl;
        if(RECharger[i][game]<CHARGER_THRESHOLD)continue;
        //cout<<"charger :"<<i<<", position :"<<chargers[i].location.first<<" "<<chargers[i].location.second<<endl;
        for(int j=0;j<NO_OF_NODES;j++){
                double d=calculateDistance(chargers[i].location,nodes[j].location);
            if( (d < CHARGER_RADIUS || double_equals(d,CHARGER_RADIUS) ) && RENode[j][game]>0 && RENode[j][game]<NODE_THRESHOLD){
                requestZoneNode[i].pb(j);

                //cout<<"node :"<<j<<", position :"<<nodes[j].location.first<<" "<<nodes[j].location.second<<endl;
                //cout<<"anu"<<endl;
            }
        }

        for(int j=0;j<NO_OF_CHARGERS;j++){
            if(i!=j){
                double d=calculateDistance(chargers[i].location,chargers[j].location);
                if( (d < CHARGER_RADIUS || double_equals(d,CHARGER_RADIUS) ) && RECharger[j][game]>0 && RECharger[j][game]<CHARGER_THRESHOLD)
                    requestZoneCharger[i].pb(j);
            }
        }
        //cout<<game<<" "<<i<<" "<<requestZoneNode[i].size()<<" "<<requestZoneCharger[i].size()<<endl;
    }
}

void calculateProfitAndStrategy(vector<wsNode> &nodes,vector<charger> &chargers,vvd &RENode,vvd &RECharger,
                                vvd &profitNode,vvd &profitCharger,vvi &requestZoneNode,vvi &requestZoneCharger,
                                vvd &priorityNode,vvd &priorityCharger,vvstrategyStruc &selectedStrategy,
                                vi &noOfServed,vvi &servedGame,vd &energyTransferred,vi &flag,int game){


        for(int i=0;i<NO_OF_CHARGERS;i++){
            //calculate total priority
            double totalPriority=0;
            double maxmProfit=-1;
            int index=-1;
            pdd target=chargers[i].location;
            int strategy=0;
            selectedStrategy[i][game]={strategy,target};
            if(RECharger[i][game]<CHARGER_THRESHOLD)continue;

            for(int j=0;j<requestZoneNode[i].size();j++){
//                    if(isnan(priorityNode[requestZoneNode[i][j]][game]))
//                        cout<<"game "<<game<<" node no: "<<j<<endl;
//                    if(game>343 && game<347){
//                        cout<<"game "<<game<<" node "<<requestZoneNode[i][j]<<" priority "<<priorityNode[requestZoneNode[i][j]][game];
//                        cout<<" RE "<<RENode[requestZoneNode[i][j]][game]<<endl;
//                    }
                    totalPriority+=priorityNode[requestZoneNode[i][j]][game];
            }

//            cout<<i<<" "<<totalPriority<<" ";
//            if(requestZoneNode[i].size()!=0)
//                cout<<endl;
            for(int j=0;j<requestZoneCharger[i].size();j++){

//                    if(isnan(priorityNode[requestZoneCharger[i][j]][game]))
//                        cout<<"game "<<game<<" charger no: "<<j<<endl;
//                    if(game>343 && game<347){
//                        cout<<"game "<<game<<" charger "<<requestZoneCharger[i][j]<<" priority "<<priorityCharger[requestZoneCharger[i][j]][game];
//                        cout<<" RE "<<RECharger[requestZoneCharger[i][j]][game]<<endl;
//                    }
                    totalPriority+=priorityCharger[requestZoneCharger[i][j]][game];
            }
//            if(requestZoneCharger[i].size()!=0)
//                cout<<endl;
//            cout<<totalPriority<<endl;
            if(double_equals(totalPriority,0.00,0.0000001))continue;

            //calculate profits
            //nodes
            for(int j=0;j<requestZoneNode[i].size();j++){
                double temp=priorityNode[requestZoneNode[i][j]][game]/totalPriority;
                profitNode[requestZoneNode[i][j]][game]=(DEL_B + DEL_U_NODE*temp);
                if(profitNode[requestZoneNode[i][j]][game]>maxmProfit){
                    strategy=1;
                    index=requestZoneNode[i][j];
                    target=nodes[index].location;
                    maxmProfit=profitNode[requestZoneNode[i][j]][game];
                }
            }

            //other chargers in the region
            for(int j=0;j<requestZoneCharger[i].size();j++){
                double temp=priorityCharger[requestZoneCharger[i][j]][game]/totalPriority;
                profitCharger[requestZoneCharger[i][j]][game]=(DEL_B + DEL_U_CHARGER*temp);
                if(profitCharger[requestZoneCharger[i][j]][game]>maxmProfit){
                    strategy=2;
                    index=requestZoneCharger[i][j];
                    target=chargers[index].location;
                    maxmProfit=profitCharger[requestZoneCharger[i][j]][game];
                }
            }

            //charger own's profit


            //save strategy

//            cout<<game<<" "<<i<<" "<<strategy<<" "<<target.first<<" "<<target.second<<endl;
            if(strategy==1 && RENode[index][game]>0.0 && RENode[index][game]<NODE_THRESHOLD){

                    energyTransferred[game]+=DEL_U_NODE-RENode[index][game];
                    RENode[index][game]=DEL_U_NODE;
                    noOfServed[game]++;
                    RECharger[i][game]-=DEL_U_NODE;
                    profitCharger[i][game]=(DEL_L + DEL_R_NODE*(priorityCharger[i][game]/totalPriority));
                    //cout<<"node "<<i<<" served game: "<<game<<endl;
                    servedGame[index].pb(game);
                    flag[index]=0;

                    selectedStrategy[i][game]={strategy,target};
                    chargers[i].location=target;

            }

            if(strategy==2 && RECharger[index][game]>0.0 && RECharger[index][game]<CHARGER_THRESHOLD){
                RECharger[index][game]=DEL_U_CHARGER;
                RECharger[i][game]-=DEL_U_CHARGER;
                profitCharger[i][game]=(DEL_L + DEL_R_CHARGER*(priorityCharger[i][game]/totalPriority));

                selectedStrategy[i][game]={strategy,target};
                chargers[i].location=target;
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


void storeRequests(vvi &requestZoneNode,vvstrategyStruc &selectedStrategy,int &noOfRequests,int &noOfServed,int game){

        vector<bool> flag(NO_OF_NODES,false);

        for(int i=0;i<NO_OF_CHARGERS;++i){

            for(int j=0;j<requestZoneNode[i].size();++j){
                if(!flag[requestZoneNode[i][j]]){
                    flag[requestZoneNode[i][j]]=true;
                    noOfRequests++;
                }
            }
        }

        for(int i=0;i<NO_OF_CHARGERS;++i){

            if(selectedStrategy[i][game].strategy==1)
                noOfServed++;
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

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void readNodeLocationFromFile(vector<wsNode> &nodes){

    FILE *fx_coord,*fy_coord;

    fx_coord=fopen("F:\\project8\\node location\\Put_NODE_X.txt","r");
    fy_coord=fopen("F:\\project8\\node location\\Put_NODE_Y.txt","r");

    float xCoordinate,yCoordinate;

    for(int i=0;i<NO_OF_NODES;++i)
    {

        fscanf(fx_coord,"%f ",&xCoordinate);
        fscanf(fy_coord,"%f ",&yCoordinate);

        wsNode *nodeObject=new wsNode(pdd(xCoordinate,yCoordinate),NODE_ENERGY_CAPACITY,0.0,0.0,0.0,0.0);
        nodes.pb(*nodeObject);
       // cout<<i<<" "<<nodes[i].location.first<<" "<<nodes[i].location.second<<endl;
    }
    nodes.pb({BASE_STATION,INT_MAX,0.0,0.0,0.0,0.0});

    fclose(fx_coord);
    fclose(fy_coord);

}


void writeNodeLocationToFile(vector<wsNode> &nodes)
{
    FILE *fn,*fx_coord,*fy_coord;
    fn=fopen("F:\\project8\\GT Charge\\Put_NO_OF_NODES.txt","w");
    fprintf(fn,"%d ",NO_OF_NODES);

    fx_coord=fopen("F:\\project8\\GT Charge\\Put_NODE_X.txt","w");
    fy_coord=fopen("F:\\project8\\GT Charge\\Put_NODE_Y.txt","w");


    for(int i=0;i<NO_OF_NODES;++i)
    {
//        if(i==21 || i==225 || i==295 || i==464 || i==841 || i==876)
//            cout<<i<<" "<<nodes[i].location.first<<" "<<nodes[i].location.second<<endl;
        fprintf(fx_coord,"%f ",nodes[i].location.first);
        fprintf(fy_coord,"%f ",nodes[i].location.second);
    }


    fclose(fn);
    fclose(fx_coord);
    fclose(fy_coord);
}


void writeNodeLocationToFileDead(vector<wsNode> &nodes,vvd &interNodecalculateDistanceMatrix)
{
    FILE *fn,*fx_coord,*fy_coord;
    fn=fopen("F:\\project8\\GT Charge\\Put_NO_OF_NODES_DEAD.txt","w");
    fprintf(fn,"%d ",10);

    fx_coord=fopen("F:\\project8\\GT Charge\\Put_NODE_X_DEAD.txt","w");
    fy_coord=fopen("F:\\project8\\GT Charge\\Put_NODE_Y_DEAD.txt","w");


    for(int i=0;i<NO_OF_NODES;++i)
    {
        if(i==21 || i==76 || i==225 || i==295 || i==464 || i==807 || i==849 || i==896 || i==962 || i==987)
//            cout<<i<<" "<<nodes[i].location.first<<" "<<nodes[i].location.second<<endl;
       {
           fprintf(fx_coord,"%f ",nodes[i].location.first);//           cout<<"node "<<i<<" distance from BS "<<interNodecalculateDistanceMatrix[i][NO_OF_NODES]<<endl;

            fprintf(fy_coord,"%f ",nodes[i].location.second);
       }
    }


    fclose(fn);
    fclose(fx_coord);
    fclose(fy_coord);
}


void writeChargerLocationToFile(vector<charger> &chargers)
{
    FILE *fn,*fx_coord,*fy_coord;
    fn=fopen("F:\\project8\\GT Charge\\Put_NO_OF_CHARGERS.txt","w");
    fprintf(fn,"%d ",NO_OF_CHARGERS);

    fx_coord=fopen("F:\\project8\\GT Charge\\Put_CHARGER_X.txt","w");
    fy_coord=fopen("F:\\project8\\GT Charge\\Put_CHARGER_Y.txt","w");


    for(int i=0;i<NO_OF_CHARGERS;++i)
    {

        fprintf(fx_coord,"%f ",chargers[i].location.first);
        fprintf(fy_coord,"%f ",chargers[i].location.second);
    }


    fclose(fn);
    fclose(fx_coord);
    fclose(fy_coord);
}




//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




void writeSelectedStrategyToFile(vvstrategyStruc &selectedStrategy)
{
    FILE *fgames,*fstrategy;
    fgames=fopen("F:\\project8\\GT Charge\\Put_NO_OF_GAMES.txt","w");
    fprintf(fgames,"%d ",NO_OF_GAMES);

    fstrategy=fopen("F:\\project8\\GT Charge\\Put_STRATEGY.txt","w");



    for(int i=0;i<NO_OF_CHARGERS;++i)
    {

            double temp_x,temp_y;

            for(int j=0;j<=NO_OF_GAMES;++j)
            {
                int st=selectedStrategy[i][j].strategy;

                    temp_x=selectedStrategy[i][j].targetLocation.first;
                    temp_y=selectedStrategy[i][j].targetLocation.second;
                    fprintf(fstrategy,"%f %f %d %d\n",temp_x,temp_y,st,i);
                    //cout<<"selected strategy of charger "<<i<<" in game "<<j<<"is "<<st<<" and target location is ";
                    //cout<<"("<<temp_x<<") , ("<<temp_y<<")"<<endl;

            }
    }


    fclose(fgames);
    fclose(fstrategy);

}


void writeRequestsToFile(vi &noOfRequests, vi &noOfServed){


    FILE *freq,*freqServed,*ftemp;
    freq=fopen("F:\\project8\\GT Charge\\Put_REQUESTS.txt","w");
    freqServed=fopen("F:\\project8\\GT Charge\\Put_REQ_SERVED.txt","w");
    ftemp=fopen("F:\\project8\\GT Charge\\Put_REQ_TEMP.txt","w");

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


void writeDeadNodesToFile(vi &deadNodes){


    FILE *fdead;
    fdead=fopen("F:\\project8\\GT Charge\\Put_NO_OF_DEAD_NODES.txt","w");

    for(int i=1;i<=NO_OF_GAMES;++i)
    {
        fprintf(fdead,"%d\n",deadNodes[i]);
    }


    fclose(fdead);
}

void writeRemEnergyNodeToFile(vvd &RE){


    FILE *fre;
    fre=fopen("F:\\project8\\GT Charge\\Put_REMENERGY.txt","w");


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


void writeRemEnergyChargerToFile(vvd &RE){


    FILE *fre;
    fre=fopen("F:\\project8\\GT Charge\\Put_REMENERGY_CHARGER.txt","w");


    for(int i=0;i<=NO_OF_GAMES;++i)
    {
        for(int j=0;j<NO_OF_CHARGERS;++j){

            fprintf(fre,"%f ",RE[j][i]);

        }
        fprintf(fre,"\n");
    }


    fclose(fre);
}


void writeResponseTimeToFile(vvi &reqGame,vvi &servedGame){

    FILE *freq,*fserved;
    freq=fopen("F:\\project8\\GT Charge\\Put_REQ_GAME.txt","w");
    fserved=fopen("F:\\project8\\GT Charge\\Put_SERVED_GAME.txt","w");

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

        if(servedGame[i].size()==0)
              fprintf(fserved,"0 ");

//        if(servedGame[i].size()<reqGame[i].size())
//            fprintf(fserved,"%d ",NO_OF_GAMES+1);

        fprintf(fserved,"\n");
    }


    fclose(freq);
    fclose(fserved);

}


void writeAverageWaitingTimeToFile(vd &waitingTime,vi &deadTime){


    FILE *fwait,*fdead;
    fwait=fopen("F:\\project8\\GT Charge\\Put_AVERAGE_WAITING.txt","w");
    fdead=fopen("F:\\project8\\GT Charge\\Put_DEAD_TIME.txt","w");

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

    fDistTravelled=fopen("F:\\project8\\GT Charge\\Put_DISTANCE_TRAVELLED.txt","w");

    for(int i=0;i<NO_OF_GAMES;++i)
    {
        cout<<distanceTravelled[i]<<endl;
        fprintf(fDistTravelled,"%f ",distanceTravelled[i]);

    }

    fclose(fDistTravelled);


}

void writeTravellingEfficiencyOfChargerToFile(vd &travellingEfficiency){

    FILE *ftravellingEfficiency;

    ftravellingEfficiency=fopen("F:\\project8\\GT Charge\\Put_TRAVELLING_EFFICIENCY.txt","w");

    for(int i=0;i<NO_OF_GAMES;++i)
    {

        fprintf(ftravellingEfficiency,"%f ",travellingEfficiency[i]);

    }

    fclose(ftravellingEfficiency);


}

void writeTravellingEfficiencyOfChargerCumulativeToFile(vd &travellingEfficiencyCumulative){

    FILE *ftravellingEfficiency;

    ftravellingEfficiency=fopen("F:\\project8\\GT Charge\\Put_TRAVELLING_EFFICIENCY_CUMULATIVE.txt","w");

    for(int i=0;i<NO_OF_GAMES;++i)
    {

        fprintf(ftravellingEfficiency,"%f ",travellingEfficiencyCumulative[i]);

    }

    fclose(ftravellingEfficiency);


}

void writeAverageEnergyTransferredToFile(vd &avgEnergyTransferred){
    FILE *favgEnergyTransferred;

    favgEnergyTransferred=fopen("F:\\project8\\GT Charge\\Put_AVG_ENERGY_TRANSFERRED.txt","w");

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

    fnetLife=fopen("F:\\project8\\comparison files\\GT Charge\\Put_NETWORK_LIFETIME_VARYING_NODES.txt","a");


    fprintf(fnetLife,"%d %f\n",NO_OF_NODES,lifetime);

    fclose(fnetLife);
}

void appendNetworkLifetimeVaryingChargersToFile(double &lifetime){

    FILE *fnetLife;

    fnetLife=fopen("F:\\project8\\comparison files\\GT Charge\\Put_NETWORK_LIFETIME_VARYING_CHARGERS.txt","a");


    fprintf(fnetLife,"%d %f\n",NO_OF_CHARGERS,lifetime);

    fclose(fnetLife);
}




//################################################################################

//temporary print

void printNodesLocation(vector<wsNode> &nodes){
    for(int i=0;i<NO_OF_NODES;i++){
        cout<<nodes[i].location.first<<" "<<nodes[i].location.second<<endl ;
    }
}

void printDijkstra(vi &parent){

    for(int i=0;i<=NO_OF_NODES;++i){
        cout<<"parent["<<i<<"] : "<<parent[i]<<endl;
    }
}

void printPathToBS(vector<wsNode> &nodes,vi &parent){

    int x=0;
    for(int i=0;i<NO_OF_NODES;++i){

        int j=i;
        cout<<"node "<<i<<": ";
        while(parent[j]!=-1){
            cout<<j<<"-->";
            j=parent[j];
            if(j==21 || j==76 || j==225 || j==295 || j==464 || j==807 || j==849 || j==896 || j==962 || j==987)x++;
        }
        if(x>0)
            cout<<" x:"<<x<<" node "<<i<<endl;
    }
}

void printRecvTransEnergy(vector<wsNode> &nodes,vvd &interNodecalculateDistanceMatrix){

    for(int i=0;i<NO_OF_NODES;++i){
//        if(i==21 || i==225 || i==295 || i==807 || i==962)
//        if(interNodecalculateDistanceMatrix[i][NO_OF_NODES]<CHARGER_RADIUS)
//        if(i==21 || i==76 || i==225 || i==295 || i==464 || i==807 || i==849 || i==896 || i==962 || i==987)
            cout<<i<<" "<<nodes[i].energyDataRecvTransmit<<" dis from BS "<<interNodecalculateDistanceMatrix[i][NO_OF_NODES]<<endl;
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
void printRENodeCharger(vvd &RENode,vvd &RECharger,int game){
    int c=0;
    for(int i=0;i<NO_OF_NODES;i++){
//        if(RENode[i][game]<0){
//            cout<<"game "<<game<<" node "<<i<<endl;
//            c++;
//        }
        cout<<"REnode "<<i<<" : "<<RENode[i][game]<<endl;
    }
//    if(game<100 && c>0)
//    cout<<"nodes dead "<<c<<endl;

//    for(int i=0;i<NO_OF_CHARGERS;i++){
//        cout<<"RECharger "<<i<<" : "<<RECharger[i][game]<<endl;
//    }
}

void printRENodeChargerNormalized(vvd &normalizedRENode,vvd &normalizedRECharger,int game){
    for(int i=0;i<NO_OF_NODES;i++){
        //cout<<"REnode "<<i<<" : "<<normalizedRENode[i][game]<<endl;
    }

    for(int i=0;i<NO_OF_CHARGERS;i++){
        cout<<"RECharger "<<i<<" : "<<normalizedRECharger[i][game]<<endl;
    }
}


void printRequestZoneNode(vvi &requestZoneNode){
    for(int i=0;i<NO_OF_CHARGERS;++i){
            cout<<i<<" ";
        for(int j=0;j<requestZoneNode[i].size();++j){
            cout<<requestZoneNode[i][j]<<" ";
        }
        cout<<endl;
    }
}


void printRequestZoneCharger(vvi &requestZoneCharger){
    for(int i=0;i<NO_OF_CHARGERS;++i){
            cout<<i<<" ";
        for(int j=0;j<requestZoneCharger[i].size();++j){
            cout<<requestZoneCharger[i][j]<<" ";
        }
        cout<<endl;
    }
}

void printContributionPriorityNode(vector<wsNode> &nodes,vvd &contributionNode,vvd &priorityNode,vvd &RENode,int game){


//    for(int i=0;i<NO_OF_NODES;++i){
//
//        cout<<"game "<<game<<" node "<<i <<": Priority "<<priorityNode[i][game]<<" Contri "<<contributionNode[i][game]<<" RE "<<RENode[i][game]<<endl;
//    }
//
    if(game>400)
        return;
    for(int i=0;i<20;++i){

        cout<<"game "<<game<<" node "<<i <<": Priority "<<priorityNode[i][game]<<" Contri "<<contributionNode[i][game]<<" RE "<<RENode[i][game]<<endl;
    }
}



void printContributionPriorityCharger(vector<charger> &chargers,vvd &contributionCharger,vvd &priorityCharger,int game){

    if(game>5)
        return;
    for(int i=0;i<NO_OF_CHARGERS;++i){

        cout<<"charger "<<i <<": Priority Contri "<<priorityCharger[i][game]<<" "<<contributionCharger[i][game]<<endl;
    }
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
int main(){
    //NETWORK DEPLOYMENT
    //Deploy Nodes
    //##############################################
	//GENERATE RANDOM LOCATION FOR WS NODES
	vector<wsNode> nodes;
	srand(time(NULL));
	//generateRandomLocation(nodes);
    readNodeLocationFromFile(nodes);


//	printNodesLocation(nodes);
	//calculateDistance between nodes
	vvd interNodecalculateDistanceMatrix(NO_OF_NODES+1,vd(NO_OF_NODES+1));
	computecalculateDistanceBetweenNodes(nodes,interNodecalculateDistanceMatrix);

	//##############################################
	//DEploy chargers initially
	std::vector<charger> chargers;
	initialChargersDeployment(chargers);

    //writing nodes and chargers location to files
    writeNodeLocationToFile(nodes);
    writeNodeLocationToFileDead(nodes,interNodecalculateDistanceMatrix);
    writeChargerLocationToFile(chargers);



	//##############################################
	//Compute path from each node to base station
	//(using Dijkstra's algo)

	//create graph for packet forwarding
	vvi graph(NO_OF_NODES+1);
	createGraph(graph,interNodecalculateDistanceMatrix);
	//printGraph(graph);

	//single source shortest path
	vi parent(NO_OF_NODES+1);
	vd dist(NO_OF_NODES+1,INT_MAX);
	dijkstraShortestPath(graph,parent,dist,interNodecalculateDistanceMatrix);

//    printDijkstra(parent);
//    cout<<interNodecalculateDistanceMatrix[0][NO_OF_NODES]<<endl;
//    printPathToBS(nodes,parent);

	//calculate data packets forwarded by each node
	computeReceivingTransmittingEnergy(parent,nodes);
	//printRecvTransEnergy(nodes,interNodecalculateDistanceMatrix);


	//#########################################################
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	//GAME THEORY



	//calculate remaining energy
	vvd RENode(NO_OF_NODES,vd(NO_OF_GAMES+1));
	vvd normalizedRENode(NO_OF_NODES,vd(NO_OF_GAMES+1));
	vvd contributionNode(NO_OF_NODES,vd(NO_OF_GAMES+1));
	vvd priorityNode(NO_OF_NODES,vd(NO_OF_GAMES+1,0.0));
    vvd profitNode(NO_OF_NODES,vd(NO_OF_GAMES+1,0.0));

	vvd RECharger(NO_OF_CHARGERS,vd(NO_OF_GAMES+1));
	vvd normalizedRECharger(NO_OF_CHARGERS,vd(NO_OF_GAMES+1));
	vvd contributionCharger(NO_OF_CHARGERS,vd(NO_OF_GAMES+1));
    vvd priorityCharger(NO_OF_CHARGERS,vd(NO_OF_GAMES+1,0.0));
    vvd profitCharger(NO_OF_CHARGERS,vd(NO_OF_GAMES+1,0.0));


	for (int i = 0; i < NO_OF_NODES; ++i)
	{
		RENode[i][0]=NODE_ENERGY_CAPACITY;
		normalizedRENode[i][0]=1;
	}



	for (int i = 0; i < NO_OF_CHARGERS; ++i)
	{
		RECharger[i][0]=CHARGER_ENERGY_CAPACITY;
		normalizedRECharger[i][0]=1;


		RECharger[i][1]=CHARGER_ENERGY_CAPACITY;
		normalizedRECharger[i][1]=1;
	}

    ////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//game has been started

	//data for final trajectory graph

	vvstrategyStruc selectedStrategy(NO_OF_CHARGERS,vstrategyStruc(NO_OF_GAMES+1));
	initializeSelectedStrategy(selectedStrategy,chargers);


    vi noOfRequests(NO_OF_GAMES+1,0);
    vi noOfServed(NO_OF_GAMES+1,0);
    vi deadNodes(NO_OF_GAMES+1,0);
    vvi reqGame(NO_OF_NODES);
    vvi servedGame(NO_OF_NODES);
    vi flag(NO_OF_NODES,0);
    vd waitingTime(NO_OF_NODES,0);
    vi deadTime(NO_OF_NODES,NO_OF_GAMES+1);
    vd energyTransferred(NO_OF_GAMES+1,0);

    for(int game=1;game<=NO_OF_GAMES;game++){

        //calculate remaining energy
		calculateRemainingEnergyOfNode(nodes,RENode,normalizedRENode,noOfRequests,deadNodes,reqGame,flag,deadTime,game);
		//first update RE of charger whenever getting charged or it charges some other charger
		calculateNormalizedREOfCharger(chargers,RECharger,normalizedRECharger,game);

		//printRENodeCharger(RENode,RECharger,game);
		//printRENodeChargerNormalized(normalizedRENode,normalizedRECharger,game);

		//find nodes in the expected zone of charger
        vvi requestZoneNode(NO_OF_CHARGERS);
        vvi requestZoneCharger(NO_OF_CHARGERS);


        createRequestZone(nodes,chargers,RENode,RECharger,requestZoneNode,requestZoneCharger,game);
//        if(game>5 &&  game<100){
//         cout<<"game "<<game<<endl;
//         printRequestZoneNode(requestZoneNode);

//        }
        //printRequestZoneCharger(requestZoneCharger);

        //compute contribution and priorities of nodes
        computeContributionOfNode(nodes,contributionNode,RENode,normalizedRENode,profitNode,game);
        computePriorityOfNode(nodes,contributionNode,priorityNode,RENode,normalizedRENode,game);

        //printContributionPriorityNode(nodes,contributionNode,priorityNode,RENode,game);


        computeContributionOfCharger(chargers,contributionCharger,RECharger,normalizedRECharger,profitCharger,game);
        computePriorityOfCharger(chargers,contributionCharger,priorityCharger,RECharger,normalizedRECharger,game);

        //printContributionPriorityCharger(chargers,contributionCharger,priorityCharger,game);


        //compute profit of nodes and chargers

        //calculate profit
        calculateProfitAndStrategy(nodes,chargers,RENode,RECharger,profitNode,profitCharger,
                            requestZoneNode,requestZoneCharger,priorityNode,priorityCharger,
                            selectedStrategy,noOfServed,servedGame,energyTransferred,flag,game);



    }

calculateAverageWaitingTime(reqGame,servedGame,waitingTime);
//
//writeSelectedStrategyToFile(selectedStrategy);
writeRequestsToFile(noOfRequests,noOfServed);
writeDeadNodesToFile(deadNodes);
//writeRemEnergyNodeToFile(RENode);
//writeRemEnergyChargerToFile(RECharger);
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

//average energy transferred
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
