#include <bits/stdc++.h>
using namespace std;
#define pdd pair<double,double>
#define pii pair<int,int>
#define vi vector<int>
#define vvi vector<vi>
#define vpii vector<pii>
#define vvpii vector<vpii>
#define pdi pair<double,int>
#define vpdi vector<pdi>
#define vd vector<double>
#define vvd vector<vd >
#define pb push_back
#define SQ(x) ((x)*(x))
#define PI 3.14159

#define NETWORK_SIZE 500
#define NETWORK_RADIUS (NETWORK_SIZE/2)
#define BASE_STATION pdd(0,0)
#define TIME_SLOT 60
#define NO_OF_GAMES 1000
#define DEAD_PERCENT 0.2


#define NO_OF_NODES 1000
#define NODE_THRESHOLD 1
#define NODE_TRANSMISSION_RANGE 25
#define NODE_ENERGY_CAPACITY 5
#define MIN_NEIGHBOUR_DIS 15

//CHARGERS
#define NO_OF_CHARGERS 16
#define CHARGER_THRESHOLD 150
#define CHARGER_ENERGY_CAPACITY 1500
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



double randomFloat(double high){
	double random = (high) * ((double)rand() / (double)RAND_MAX );
    return random;
}

double calculateDistance(pdd a,pdd b){
	return sqrt(SQ(a.first-b.first)+SQ(a.second-b.second));
}




void calculateDistanceBetweenNodes(vector<wsNode> &nodes,
	vector<vector<double> > &interNodeDistanceMatrix){

	for (int i = 0; i < NO_OF_NODES+1; ++i)//+1 for the base station
	{
		for (int j = 0; j < NO_OF_NODES+1; ++j)
		{
			interNodeDistanceMatrix[i][j]=calculateDistance(nodes[i].location,nodes[j].location);
//			cout<<interNodeDistanceMatrix[i][j]<<" ";
		}
//		cout<<endl;
	}

}

void createGraph(vvi &graph,vvd &interNodeDistanceMatrix ){

	for (int i = 0; i < NO_OF_NODES+1; ++i)
	{
		for (int j = 0; j < NO_OF_NODES+1; ++j)
		{
			if(i!=j && interNodeDistanceMatrix[i][j] < NODE_TRANSMISSION_RANGE){
				graph[i].pb(j);
				//graph[j].pb(i);
				//cout<<i<<" "<<j<<endl;
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



void dijkstraShortestPath(vvi &graph,vi &parent,vd &dist,vvd &interNodeDistanceMatrix){
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
			if(dist[j]>dist[currNode]+interNodeDistanceMatrix[currNode][j]){
				parent[j]=currNode;
				dist[j]=dist[currNode]+interNodeDistanceMatrix[currNode][j];
				heap.push(pdi(dist[j],j));
			}
		}
	}

	//for(int i=0;i<NO_OF_NODES+1;i++)cout<<parent[i]<<" ";
	//cout<<endl;

}


bool double_equals(double a, double b, double epsilon = 0.0000001)
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

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//game

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
		double xCoordinate=(radius/2)*cos(theta*i);
		double yCoordinate=(radius/2)*sin(theta*i);
		//cout<<xCoordinate<<" "<<yCoordinate<<endl;
		chargers.pb({pdd(xCoordinate,yCoordinate),CHARGER_ENERGY_CAPACITY});
	}
}

void computeNodeChargerDistance(vector<wsNode> &nodes,vector<charger> &chargers,vvd &nodeChargerDistance){
    for(int i=0;i<NO_OF_NODES+1;i++){
        for(int j=0;j<NO_OF_CHARGERS;j++){
            nodeChargerDistance[i][j]=calculateDistance(nodes[i].location, chargers[j].location);
        }
    }
}


void formClusters(vector<charger> &chargers,vvd &nodeChargerDistance,vvi &clusters){
    for(int i=0;i<NO_OF_NODES;i++){
        double minm=INT_MAX;
        int index=-1;
        for(int j=0;j<NO_OF_CHARGERS;j++){
            double prod=(1+(nodeChargerDistance[i][j]/(2*NETWORK_RADIUS)));
            prod=prod*(2-(chargers[j].remainingEnergy/CHARGER_ENERGY_CAPACITY));
           // cout<<i<<" "<<j<<" "<<prod<<endl;
            if(prod<minm){
                minm=prod;
                index=j;
            }
        }

        clusters[index].pb(i);
    }
}

void updateRENodeDeadNodes(vector<wsNode> &nodes,vvd &RENode,vi &deadNodes,int game){
    for(int i=0;i<NO_OF_NODES;i++){
        if(game==1 || (game>1 && !double_equals(RENode[i][game],NODE_ENERGY_CAPACITY))){
            RENode[i][game]=RENode[i][game-1]-nodes[i].energyDataRecvTransmit;
            //if(game ==999)
           // cout<<RENode[i][game]<<" ";
            if(double_equals(RENode[i][game],0) || (RENode[i][game]<0)){
                    deadNodes[game]++;cout<<RENode[i][game]<<" ";
            }
        }
    }
    //if(game==999)
    //cout<<endl;
}






void computeChargingTrajectory(vector<charger> &chargers,vector<wsNode> &nodes,vvd &nodeChargerDistance,vvd &interNodeDistanceMatrix,
                               vvi &clusters,vvi &trajectory,vvd &energyTransferredChargers,vvd &RENode,vi &nodesCharged,int game){

    for(int i=0;i<NO_OF_CHARGERS;i++){
        double minm=INT_MAX,timeInterval=0;
        int index=-1;

        double remEnergyNode;

        for(int j=0;j<clusters[i].size();++j){
            int nodeIndex=clusters[i][j];

            //if the node is already dead
            if(RENode[nodeIndex][game]<0 || double_equals(RENode[nodeIndex][game],0))
                continue;

            //if the trajectory for this game has already been decided in some previous game
            if(trajectory[i][game]!=-1)continue;

            //this is the node charged just now so check to move to other nodes
            if(trajectory[i][game-1]==nodeIndex)continue;

            //double RENode=nodes[nodeIndex].remainingEnergy;



            double t;
            if(game==1)
                t=nodeChargerDistance[nodeIndex][i]/CHARGER_SPEED;
            else
                t=interNodeDistanceMatrix[trajectory[i][game-1]][nodeIndex]/CHARGER_SPEED;

            double energy=RENode[nodeIndex][game]-nodes[nodeIndex].energyDataRecvTransmit*ceil(t/TIME_SLOT);
//            double energy=0;
            if(game+ceil(t/TIME_SLOT) > NO_OF_GAMES)continue;

            //if charger selects this node,energy of node will become<=0 before charger reaches it
            if(energy<0 || double_equals(energy,0))
                continue;



            double prod=(1+(nodeChargerDistance[nodeIndex][i]/(2*NETWORK_RADIUS)));
                prod*=(1+(RENode[nodeIndex][game]/NODE_ENERGY_CAPACITY));
//             cout<<nodeIndex<<" "<<prod<<" "<<t<<" "<<energy<<endl;
            if(prod<minm){
                    minm=prod;
                    index=nodeIndex;
                    timeInterval=t;
                    remEnergyNode=energy;
                   // cout<<index<<" "<<minm<<" "<<timeInterval<<endl;
            }

        }

        //update trajectory
        if(index!=-1){
            int k;
            for(k=game;k<=game+floor(timeInterval/TIME_SLOT) && k<=NO_OF_GAMES;k++){
                trajectory[i][k]=index;
            }
            RENode[index][k-1]=NODE_ENERGY_CAPACITY;
            energyTransferredChargers[i][k-1]=NODE_ENERGY_CAPACITY-remEnergyNode;
            chargers[i].remainingEnergy-=(NODE_ENERGY_CAPACITY-remEnergyNode);
            //cout<<index<<" "<<minm<<" "<<timeInterval<<endl;
            nodesCharged[k]++;
        }
    }
}



void computeChargedNodesDeadNodes(vvd &timeOfCharging,vd &deadTime,vi &nodesCharged,vi &nodesDead){
    for(int i=0;i<NO_OF_CHARGERS;i++){
        for(int j=0;j<timeOfCharging[i].size();j++){
            int counter=timeOfCharging[i][j]/TIME_SLOT;
            nodesCharged[counter]++;
        }
    }

    for(int i=0;i<NO_OF_NODES;i++){
        int counter=deadTime[i]/TIME_SLOT;
        nodesDead[counter]++;
    }
    for(int i=1;i<NO_OF_GAMES;i++){
        nodesDead[i]+=nodesDead[i-1];
    }
}

void calculateDistanceTravelledByCharger(vvd &timeOfCharging,vd &distanceTravelled){

    for(int i=0;i<NO_OF_CHARGERS;++i){


        int lastChargedGame=timeOfCharging[i].back()/TIME_SLOT;
        double distLeft=(timeOfCharging[i].back()-lastChargedGame*TIME_SLOT);
        int j;

        for(j=0;j<min(lastChargedGame,NO_OF_GAMES);j++){
            distanceTravelled[j]+=(TIME_SLOT*CHARGER_SPEED);
        }
        if(j<NO_OF_GAMES-1)
            distanceTravelled[j+1]+=distLeft;
    }

    for(int i=0;i<NO_OF_GAMES;++i){
        distanceTravelled[i]/=NO_OF_CHARGERS;
    }

}

void calculateTravellingEfficiency(vi &nodesCharged,vd &distanceTravelled,vd &travellingEfficiency){
    for(int i=0;i<NO_OF_GAMES;i++){
        if(double_equals(distanceTravelled[i],0.0))travellingEfficiency[i]=0;
        else
            travellingEfficiency[i]=nodesCharged[i]/distanceTravelled[i];
    }

}

void calculateTravellingEfficiencyCumulative(vi &nodesCharged,vd &distanceTravelled,vd &travellingEfficiencyCumulative){
    int totalNodesCharged=0;
    double totalDistTravelled=0;
    for(int i=0;i<NO_OF_GAMES;i++){
        totalDistTravelled+=distanceTravelled[i];
        totalNodesCharged+=nodesCharged[i];
        if(double_equals(totalDistTravelled,0.0))travellingEfficiencyCumulative[i]=0;
        else
            travellingEfficiencyCumulative[i]=totalNodesCharged/totalDistTravelled;
    }

}

void calculateAverageEnergyTransferred(vvd &energyTransferredChargers,vd &avgEnergyTransferred){



    for(int i=1;i<=NO_OF_GAMES;++i){
        double sum=0;
        for(int j=0;j<NO_OF_CHARGERS;j++){
            sum+=energyTransferredChargers[j][i];
        }
        avgEnergyTransferred[i]=sum/NO_OF_CHARGERS;
    }
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
    fn=fopen("F:\\project8\\CCGK\\Put_NO_OF_NODES.txt","w");
    fprintf(fn,"%d ",NO_OF_NODES);

    fx_coord=fopen("F:\\project8\\CCGK\\Put_NODE_X.txt","w");
    fy_coord=fopen("F:\\project8\\CCGK\\Put_NODE_Y.txt","w");


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
    fn=fopen("F:\\project8\\CCGK\\Put_NO_OF_CHARGERS.txt","w");
    fprintf(fn,"%d ",NO_OF_CHARGERS);

    fx_coord=fopen("F:\\project8\\CCGK\\Put_CHARGER_X.txt","w");
    fy_coord=fopen("F:\\project8\\CCGK\\Put_CHARGER_Y.txt","w");


    for(int i=0;i<NO_OF_CHARGERS;++i)
    {

        fprintf(fx_coord,"%f ",chargers[i].location.first);
        fprintf(fy_coord,"%f ",chargers[i].location.second);
    }


    fclose(fn);
    fclose(fx_coord);
    fclose(fy_coord);
}

void writeChargedNodesToFile(vi &nodesCharged){
    FILE *fNodesCharged;

    fNodesCharged=fopen("F:\\project8\\CCGK\\Put_NODES_CHARGED.txt","w");



    for(int i=0;i<NO_OF_GAMES;++i)
    {

        fprintf(fNodesCharged,"%d ",nodesCharged[i]);

    }


    fclose(fNodesCharged);

}

void writeDeadNodesToFile(vi &nodesDead){
    FILE *fNodesDead;

    fNodesDead=fopen("F:\\project8\\CCGK\\Put_NODES_DEAD.txt","w");



    for(int i=0;i<NO_OF_GAMES;++i)
    {

        fprintf(fNodesDead,"%d ",nodesDead[i]);

    }


    fclose(fNodesDead);

}


void writeDistanceTravelledByChargerToFile(vd &distanceTravelled){

    FILE *fDistTravelled;

    fDistTravelled=fopen("F:\\project8\\CCGK\\Put_DISTANCE_TRAVELLED.txt","w");

    for(int i=0;i<NO_OF_GAMES;++i)
    {

        fprintf(fDistTravelled,"%f ",distanceTravelled[i]);

    }

    fclose(fDistTravelled);


}

void writeTravellingEfficiencyOfChargerToFile(vd &travellingEfficiency){

    FILE *ftravellingEfficiency;

    ftravellingEfficiency=fopen("F:\\project8\\CCGK\\Put_TRAVELLING_EFFICIENCY.txt","w");

    for(int i=0;i<NO_OF_GAMES;++i)
    {

        fprintf(ftravellingEfficiency,"%f ",travellingEfficiency[i]);

    }

    fclose(ftravellingEfficiency);


}

void writeTravellingEfficiencyOfChargerCumulativeToFile(vd &travellingEfficiencyCumulative){

    FILE *ftravellingEfficiency;

    ftravellingEfficiency=fopen("F:\\project8\\CCGK\\Put_TRAVELLING_EFFICIENCY_CUMULATIVE.txt","w");

    for(int i=0;i<NO_OF_GAMES;++i)
    {

        fprintf(ftravellingEfficiency,"%f ",travellingEfficiencyCumulative[i]);

    }

    fclose(ftravellingEfficiency);


}


void writeAverageEnergyTransferredToFile(vd &avgEnergyTransferred){

    FILE *favgEnergyTransferred;

    favgEnergyTransferred=fopen("F:\\project8\\CCGK\\Put_AVG_ENERGY_TRANSFERRED.txt","w");

    for(int i=1;i<=NO_OF_GAMES;++i)
    {

        fprintf(favgEnergyTransferred,"%f ",avgEnergyTransferred[i]);

    }

    fclose(favgEnergyTransferred);


}

//################################################################################
//GRAPHS files
void appendNetworkLifetimeVaryingNodesToFile(double &lifetime){
    FILE *fnetLife;

    fnetLife=fopen("F:\\project8\\comparison files\\CCGK\\Put_NETWORK_LIFETIME_VARYING_NODES.txt","a");

    fprintf(fnetLife,"%d %f\n",NO_OF_NODES,lifetime);

    fclose(fnetLife);
}

void appendNetworkLifetimeVaryingChargersToFile(double &lifetime){

    FILE *fnetLife;

    fnetLife=fopen("F:\\project8\\comparison files\\CCGK\\Put_NETWORK_LIFETIME_VARYING_CHARGERS.txt","a");

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

void printNodeChargerDistance(vvd &nodeChargerDistance){
    for(int i=0;i<NO_OF_NODES+1;i++){
//        double dis=INT_MAX;
//        int index=-1;
        for(int j=0;j<NO_OF_CHARGERS;j++){
            cout<<nodeChargerDistance[i][j]<<" ";
//            if(nodeChargerDistance[i][j]<dis){
//                dis=nodeChargerDistance[i][j];
//                index=j;
//            }
        }
//        cout<<i<<" "<<index<<" "<<dis<<endl;
        cout<<endl;
    }
}

void printClusters(vvi &clusters){

    for(int i=0;i<NO_OF_CHARGERS;++i){

        cout<<i<<" : ";
        for(int j=0;j<clusters[i].size();++j){
                cout<<clusters[i][j]<<" ";
        }
        cout<<endl;
    }
}

void printChargingTrajectory(vvi &trajectory){

    for(int i=0;i<NO_OF_CHARGERS;++i){

        for(int j=1;j<=NO_OF_GAMES;++j){
            cout<<trajectory[i][j]<<" -> ";
        }
        cout<<endl;
    }
}

void printNodeDeadTime(vd &deadTime){

    for(int i=0;i<NO_OF_NODES;++i){
        cout<<i<<" "<<deadTime[i]<<endl;
    }
}


void printRENodes(vector<wsNode> &nodes){

    for(int i=0;i<NO_OF_NODES;++i){
        cout<<i<<" "<<nodes[i].remainingEnergy<<endl;
    }
}


void printChargedNodesDeadNodes(vi &nodesCharged,vi &nodesDead){
    for(int i=0;i<30;i++){
        cout<<i<<" nodesCharged "<<nodesCharged[i]<<"    nodesDead "<<nodesDead[i]<<endl;
    }
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




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
	vvd interNodeDistanceMatrix(NO_OF_NODES+1,vd(NO_OF_NODES+1,0));
	calculateDistanceBetweenNodes(nodes,interNodeDistanceMatrix);

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
	createGraph(graph,interNodeDistanceMatrix);
    //printGraph(graph);

	//single source shortest path
	vector<int> parent(NO_OF_NODES+1);
	vector<double> dist(NO_OF_NODES+1,INT_MAX);
	dijkstraShortestPath(graph,parent,dist,interNodeDistanceMatrix);
    //printDijkstra(parent);

	//calculate data packets forwarded by each node
	computeReceivingTransmittingEnergy(parent,nodes);
	//printRecvTransEnergy(nodes);




	//################################################################################
	vvi trajectory(NO_OF_CHARGERS,vi(NO_OF_GAMES+1,-1));
    vvd timeOfCharging(NO_OF_CHARGERS,vd(NO_OF_GAMES+1));
    vvd energyTransferredChargers(NO_OF_CHARGERS,vd(NO_OF_GAMES+1));
    vvd RENode(NO_OF_NODES,vd(NO_OF_GAMES+1));
    vi nodesCharged(NO_OF_GAMES+1,0);
    vi nodesDead(NO_OF_GAMES+1,0);

    for(int i=0;i<NO_OF_NODES;++i){
        RENode[i][0]=NODE_ENERGY_CAPACITY;
    }



	for(int game=1;game<=NO_OF_GAMES;game++){
        //calculate node charger distance
        vvd nodeChargerDistance(NO_OF_NODES+1,vd(NO_OF_CHARGERS));
        computeNodeChargerDistance(nodes,chargers,nodeChargerDistance);
        //printNodeChargerDistance(nodeChargerDistance);

        //coordination algorithm
        //form clusters
        vvi clusters(NO_OF_CHARGERS);
        formClusters(chargers,nodeChargerDistance,clusters);
        //printClusters(clusters);

        //update RE of nodes
        updateRENodeDeadNodes(nodes,RENode,nodesDead,game);
        computeChargingTrajectory(chargers,nodes,nodeChargerDistance,interNodeDistanceMatrix,
                              clusters,trajectory,energyTransferredChargers,RENode,nodesCharged,game);


	}
	//printChargingTrajectory(trajectory);

	//energy level of the network
	vd avgEnergyTransferred(NO_OF_GAMES+1,0);
    calculateAverageEnergyTransferred(energyTransferredChargers,avgEnergyTransferred);
    writeAverageEnergyTransferredToFile(avgEnergyTransferred);

    //chargedNodes and deadNodes
    writeChargedNodesToFile(nodesCharged);
    writeDeadNodesToFile(nodesDead);

    //calculateDistanceTravelledByCharger(distanceTravelled);
    //writeDistanceTravelledByChargerToFile(distanceTravelled);

	return 0;
}
