clear;
axis([-250 250 -250 250]);
hold on;
plot(0,0,'bs','LineWidth',3);



noOfNodes=textread('F:\\project8\\node location\\Put_NO_OF_NODES.txt','%d');
disp(noOfNodes);


node_X=textread('F:\\project8\\node location\\Put_NODE_X.txt','%f',noOfNodes);
node_Y=textread('F:\\project8\\node location\\Put_NODE_Y.txt','%f',noOfNodes);

plot(node_X,node_Y,'r.','LineWidth',3); 

%ANOMALY
% plot(node_X(624),node_Y(624),'b.','LineWidth',3); 
% noOfNodesDead=textread('Put_NO_OF_NODES_DEAD.txt','%d');
% disp(noOfNodesDead);
% 
% 
% node_X=textread('Put_NODE_X_DEAD.txt','%f',noOfNodesDead);
% node_Y=textread('Put_NODE_Y_DEAD.txt','%f',noOfNodesDead);
% 
% plot(node_X,node_Y,'b*','LineWidth',3); 


%GTcharge
%    OR
%GTcharge Modified
no_of_chargers=textread('F:\\project8\\GT Charge\\Put_NO_OF_CHARGERS.txt','%d');
disp(no_of_chargers);
charger_X=textread('F:\\project8\\GT Charge\\Put_CHARGER_X.txt','%f',no_of_chargers);
charger_Y=textread('F:\\project8\\GT Charge\\Put_CHARGER_Y.txt','%f',no_of_chargers);

plot(charger_X,charger_Y,'gV','LineWidth',3); 



%CCGK
no_of_chargers=textread('F:\\project8\\CCGK\\Put_NO_OF_CHARGERS.txt','%d');
disp(no_of_chargers);
charger_X=textread('F:\\project8\\CCGK\\Put_CHARGER_X.txt','%f',no_of_chargers);
charger_Y=textread('F:\\project8\\CCGK\\Put_CHARGER_Y.txt','%f',no_of_chargers);

plot(charger_X,charger_Y,'gV','LineWidth',3); 



   