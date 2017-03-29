#include<bits/stdc++.h>
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
#define NO_OF_GAMES 1000

#define NO_OF_NODES 500
#define NODE_THRESHOLD 1
#define NODE_TRANSMISSION_RANGE 25
#define NODE_ENERGY_CAPACITY 5
#define MIN_NEIGHBOUR_DIS 15



using namespace std;


typedef struct wsNode
{
	pdd location;
	double remainingEnergy;
	double energyDataRecvTransmit; //receiving + transmission
	double energyForBroadcast;

	wsNode(pdd loc,double remEnergy,double enRT,double enB) :
	    location(loc), remainingEnergy(remEnergy), energyDataRecvTransmit(enRT), energyForBroadcast(enB){}

}wsNode;



double randomFloat(double low,double high){
	double random = low+((high-low) * ((double)rand() / (double)RAND_MAX ));
    return random;
}

double calculateDistance(pdd a,pdd b){
	return sqrt(SQ(a.first-b.first)+SQ(a.second-b.second));
}

int inRange(double x,double y){
    if(x<(-NETWORK_RADIUS) || x>NETWORK_RADIUS || y<(-NETWORK_RADIUS) || y>NETWORK_RADIUS)return 0;
    return 1;
}

//void generateRandomLocation(std::vector<wsNode> &nodes){
//	int nodeCount=NO_OF_NODES;
//	while(nodeCount){
//		int notPresentFlag=1;
//		double radius = randomFloat(0,NETWORK_RADIUS);
//		double angle = randomFloat(0,2*PI);
//
//		double xCoordinate=radius*cos(angle);
//		double yCoordinate=radius*sin(angle);
//
//
//		for (int i = 0; i < nodes.size(); ++i)
//		{
//			if(calculateDistance(pdd(xCoordinate,yCoordinate) , nodes[i].location) < MIN_NEIGHBOUR_DIS
//				&& calculateDistance(pdd(xCoordinate,yCoordinate) , BASE_STATION) < MIN_NEIGHBOUR_DIS
//                &&inRange(xCoordinate,yCoordinate)){
//				notPresentFlag=0;
//
//				break;
//			}
//		}
//
//		if(notPresentFlag){
//			wsNode *nodeObject=new wsNode(pdd(xCoordinate,yCoordinate),NODE_ENERGY_CAPACITY,0.0,0.0);
//			nodes.pb(*nodeObject);
//			//cout<<xCoordinate<<" "<<yCoordinate<<endl;
//			nodeCount--;
//			//notPresentFlag=1;
//		}
//	}
//
//	//n+1 th is the BASE STATION
//	nodes.pb({BASE_STATION,INT_MAX,0.0,0.0});
//}


void generateRandomLocation(std::vector<wsNode> &nodes){
	int nodeCount=NO_OF_NODES;
	while(nodeCount){
		int notPresentFlag=1;
		//double radius = randomFloat(NETWORK_RADIUS);
		//double angle = randomFloat(2*PI);

		double xCoordinate=randomFloat(-NETWORK_RADIUS,NETWORK_RADIUS);
		double yCoordinate=randomFloat(-NETWORK_RADIUS,NETWORK_RADIUS);
		double d=calculateDistance(pdd(0,0),pdd(xCoordinate,yCoordinate));
		while(d > NETWORK_RADIUS){
            xCoordinate=randomFloat(-NETWORK_RADIUS,NETWORK_RADIUS);
            yCoordinate=randomFloat(-NETWORK_RADIUS,NETWORK_RADIUS);
            d=calculateDistance(pdd(0,0),pdd(xCoordinate,yCoordinate));
		}


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
			nodeCount--;
			//notPresentFlag=1;
		}
	}

	//n+1 th is the BASE STATION
	nodes.pb({BASE_STATION,INT_MAX,0.0,0.0});
}

void writeNodeLocationToFile(vector<wsNode> &nodes)
{
    FILE *fn,*fx_coord,*fy_coord;
    fn=fopen("F:\\project8\\node location\\Put_NO_OF_NODES.txt","w");
    fprintf(fn,"%d ",NO_OF_NODES);

    fx_coord=fopen("F:\\project8\\node location\\Put_NODE_X.txt","w");
    fy_coord=fopen("F:\\project8\\node location\\Put_NODE_Y.txt","w");


    for(int i=0;i<NO_OF_NODES;++i)
    {

        fprintf(fx_coord,"%f ",nodes[i].location.first);
        fprintf(fy_coord,"%f ",nodes[i].location.second);
    }


    fclose(fn);
    fclose(fx_coord);
    fclose(fy_coord);
}


int main(int argc, char const *argv[]){

    vector<wsNode> nodes;
	srand(time(NULL));
	generateRandomLocation(nodes);

    writeNodeLocationToFile(nodes);

    return 0;
}
