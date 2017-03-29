clear;


noOfGames=textread('F:\\project8\\GT Charge\\Put_NO_OF_GAMES.txt','%d');
noOfNodes=textread('F:\\project8\\node location\\Put_NO_OF_NODES.txt','%d');


axis([0 noOfGames 0 noOfNodes/25]);
hold on;


distanceGTCharge=textread('F:\\project8\\GT Charge\\Put_REQ_SERVED.txt','%d');
distanceGTChargeModified=textread('F:\\project8\\GT Charge Modified\\Put_REQ_SERVED.txt','%d');
distanceCCGK=textread('F:\\project8\\CCGK\\Put_NODES_CHARGED.txt','%d');

plot(distanceGTCharge,'r-','LineWidth',1.5);
plot(distanceGTChargeModified,'g-','LineWidth',0.5);
plot(distanceCCGK,'b-','LineWidth',0.5);


   