clear;


noOfGames=textread('F:\\project8\\GT Charge\\Put_NO_OF_GAMES.txt','%d');
noOfNodes=textread('F:\\project8\\node location\\Put_NO_OF_NODES.txt','%d');


axis([0 noOfGames 0 noOfNodes]);
hold on;

chargedNodes=textread('F:\\project8\\CCGK\\Put_NODES_CHARGED.txt','%d');
deadNodes=textread('F:\\project8\\CCGK\\Put_NODES_DEAD.txt','%d');

plot(chargedNodes,'g-','LineWidth',1.5);
plot(deadNodes,'r-','LineWidth',1.5);


   