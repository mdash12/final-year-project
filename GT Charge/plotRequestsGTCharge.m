clear;

no_of_games=textread('Put_NO_OF_GAMES.txt','%d');
disp(no_of_games);



no_of_nodes=textread('Put_NO_OF_NODES.txt','%d');
disp(no_of_nodes);

axis([0 200 0 120]);
hold on;


requests=textread('Put_REQUESTS.txt','%d');
reqServed=textread('Put_REQ_SERVED.txt','%d');

% disp(requests);
names={'requestsGTCharge','requestsServedGTCharge','deadNodesGTCharge','requestsGTChargeModified','requestsServedGTChargeModified','deadNodesGTChargeModified'};

h(1)=plot(requests,'b-','LineWidth',0.5);

hold on;
h(2)=plot(reqServed,'g-','LineWidth',0.5);

no_of_dead_nodes=textread('Put_NO_OF_DEAD_NODES.txt','%d');
h(3)=plot(no_of_dead_nodes,'c-','LineWidth',0.5);

requestsModified=textread('F:\\project8\\GT Charge Modified\\Put_REQUESTS.txt','%d');
reqServedModified=textread('F:\\project8\\GT Charge Modified\\Put_REQ_SERVED.txt','%d');

h(4)=plot(requestsModified,'m-','LineWidth',0.5);
hold on;

h(5)=plot(reqServedModified,'r-','LineWidth',0.5);

no_of_dead_nodes=textread('F:\\project8\\GT Charge Modified\\Put_NO_OF_DEAD_NODES.txt','%d');
h(6)=plot(no_of_dead_nodes,'y-','LineWidth',0.5);
legend(h,names,'Location','NorthWest');