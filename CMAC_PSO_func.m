%  Cerebellar model articulation controller based Particle Swarm Optimization Algorithm(CMAC-PSO)... 
%% This code was developed by Nazmiye Ebru Bulut, Alpaslan Duysak, Emre Dandil & Ugur Yuzgec
%% -----------------------------------------------------------------------------------------------------------------------------

function [gbest,gbestval,fitcount,t_c]= CMAC_PSO_func(fhd,D,ps,Max_iteration,VRmin,VRmax,X_suru,A,w1,w2,varargin)

rand('state',sum(100*clock));

gbest_position=zeros(D*10000/ps,D);

cc=[2 2];   %acceleration constants

iwt=0.9-(1:Max_iteration).*(0.5./Max_iteration);
% iwt=0.5.*ones(1,Max_iteration);
if length(VRmin)==1
    VRmin=repmat(VRmin,1,D); 
    VRmax=repmat(VRmax,1,D);
end
mv=0.5*(VRmax-VRmin);
lu = [VRmin;VRmax];
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);
Vmin=repmat(-mv,ps,1);
Vmax=-Vmin;

    
    % Generate the points
    pos = generate_points(ps, D, lu);
% pos=X_suru;
 
    fit =feval(fhd,pos',varargin{:});
	fitcount=ps;
    vel=Vmin+2.*Vmax.*rand(ps,D);%initialize the velocity of the particles
    pbest=pos;
    pbestval=fit; %initialize the pbest and the pbest's fitness value
    [gbestval,gbestid]=min(pbestval);
    gworstval=max(pbestval);
    gbest=pbest(gbestid,:);%initialize the gbest and the gbest's fitness value
    gbestrep=repmat(gbest,ps,1);
    gbest_position(1,:)=gbest;
    Distance = abs(pbest(1,1)-pbest(:,1)); % this is for analysis
    DistanceSum(1) = sum(Distance)/(ps-1);
    
    pTemp=pos;
    pTemp2=pos;
	for j=1:D
		for i=1:ps
			pTemp2(i,j)=pTemp(i,j)+abs(min(pTemp(:,j)));
		end      
		if max(pTemp2(:,j))>0
			pTemp2(:,j)=pTemp2(:,j)/max(pTemp2(:,j));          
		else
			pTemp2(:,j)=0;
		end
	end
        
         % Population diversity as a whole
        temp=[];
        for j=1:D
            temp(j)=mean(abs(pTemp2(:,j)-mean(pTemp2(:,j))));           
        end
        Div(1)= sum(temp)/ps;
        
        
%     [fbest,gbestid]=min(fit);
%     gbest=pos(gbestid,:);%initialize the gbest and the gbest's fitness value
% 	fworst = max(fit);
%     gbest_position(1,:) = gbest;
    Distance = abs(pos(1,1)-pos(:,1)); % for analysis
    DistanceSum(1) = sum(Distance)/(ps-1);
   
    % 2. Search History Analysis
if D==2 % to be done for only 2-dimension...
	refresh = 100;
	figure(11)

% initialize
% % %         x1 = linspace(lowerbound(1), upperbound(1), 101);
% % %         x2 = linspace(lowerbound(2), upperbound(2), 101);
% % %         x3 = zeros(length(x1), length(x2));
% % % 	for i = 1:length(x1)
% % % 		for j = 1:length(x2)
% % % 			x3(i,j) = feval(fhd,[x1(i);x2(j)],varargin{:});
% % % 		end
% % %     end
        
% tiles, labels, legend
% 	str = sprintf('Search History of FN%d',varargin{:});
% 	xlabel('x_1'); ylabel('x_2');  title(str);
% 	contour(x1', x2', x3'); 
    hold on;
	plot(pos(:,1),pos(:,2),'ro','MarkerSize',8,'MarkerFaceColor','k');
	drawnow
	% plot(gbest(1),gbest(2),'k*','MarkerSize',8);
	figure(13)
    %%% contour(x1', x2', x3'); 
    hold on; 
    plot(gbest(1),gbest(2),'ro','MarkerSize',8,'MarkerFaceColor','k');
end
  
    %Main loop
    t_c=[];t=1; m=100; %number of unsuccessful attempts to improve the fitness
    counter=zeros(ps,1); % the counter is for unsuccessful attempts
     % start from the second iteration since the first iteration was dedicated to calculating the fitness
while fitcount<D*10000
    %------------------------------------------------------------------------
    for k=1:ps
         r1(k,:)=mean(w1(1+(k-1)*A:A+(k-1)*A,:));  %calculate mean of the weights
         r2(k,:)=mean(w2(1+(k-1)*A:A+(k-1)*A,:));  %calculate mean of the weights
    end
    aa=cc(1).*r1.*(pbest-pos)+cc(2).*r2.*(gbestrep-pos);
    %------------------------------------------------------------------------
    vel=iwt(i).*vel+aa;
    vel=(vel>Vmax).*Vmax+(vel<=Vmax).*vel;
    vel=(vel<Vmin).*Vmin+(vel>=Vmin).*vel;
    pos=pos+vel;
 
    %% Check if a dimension of the candidate solution goes out of boundaries
	%------------------------------------------------------------------------
	pos = boundConstraint(pos,pbest,lu); % new boundary mechanism...
	%------------------------------------------------------------------------ 
    fit=feval(fhd,pos',varargin{:});
    fitcount=fitcount+ps;
    tmp=(pbestval<fit);
    temp=repmat(tmp',1,D);
    pbest=temp.*pbest+(1-temp).*pos;
    pbestval=tmp.*pbestval+(1-tmp).*fit;%update the pbest
    [gbestval,tmp]=min(pbestval);
    gworstval=max(pbestval);
    gbest=pbest(tmp,:);
    gbestrep=repmat(gbest,ps,1);%update the gbest
	
%     fworst = max(fit);
%     [fbest,index] = min(fit);
%     gbest = pos(index,:); 
	n=0.001;
	%% training weights...
        for i=1:20
            for k=1:ps
%                 for r=1:D
                w1(1+(k-1)*A:A+(k-1)*A,:)=w1(1+(k-1)*A:A+(k-1)*A,:)+...
                               n*((gbestval-pbestval(k))/(gbestval-gworstval))/A;
                w2(1+(k-1)*A:A+(k-1)*A,:)=w2(1+(k-1)*A:A+(k-1)*A,:)+...
                               n*((gbestval-pbestval(k))/(gbestval-gworstval))/A;
%                 end
            end
        end
        
	if fitcount==100*D||fitcount==200*D||fitcount==300*D||fitcount==500*D...
		||fitcount==1000*D||fitcount==2000*D||fitcount==3000*D||fitcount==4000*D||fitcount==5000*D...
		||fitcount==6000*D||fitcount==7000*D||fitcount==8000*D||fitcount==9000*D||fitcount==10000*D
			t_c=[t_c;abs(gbestval-varargin{:}*100)];              
    end
    % Increase the iteration counter
	t=t+1;

%------------------------------------------------------------------------------
if D==2    
    pTemp=pbest;
    pTemp2=pbest;
	for j=1:D
		for i=1:ps
			pTemp2(i,j)=pTemp(i,j)+abs(min(pTemp(:,j)));
		end      
		if max(pTemp2(:,j))>0
			pTemp2(:,j)=pTemp2(:,j)/max(pTemp2(:,j));          
		else
			pTemp2(:,j)=0;
		end
	end
        
    % Population diversity as a whole
	temp=[];
	for j=1:D
		temp(j)=mean(abs(pTemp2(:,j)-mean(pTemp2(:,j))));           
	end
	Div(t)= sum(temp)/ps;
end
%------------------------------------------------------------------------------

    gbest_position(t,:)=gbest;
% % %     disp(['Iteration ' num2str(t) ': Best Cost = ' num2str(gbestval)]);

% % %     gbestval=fbest;
% % %     gbest=gbest;
    
		Distance = abs(pbest(1,1)-pbest(:,1)); % for analysis
		DistanceSum(t) = sum(Distance)/(ps-1);
        % 2. Search History Analysis
		if D==2 && (rem(fitcount/ps,refresh) == 0)% to be done for only 2-dimension...
			figure(11)
% 			plot(pos(:,1),pos(:,2),'ro','MarkerSize',8,'MarkerFaceColor','k');
		%	plot(gbest(1),gbest(2),'k*','MarkerSize',8);
% % % 		xlabel('x_1'); ylabel('x_2'); title(str);
			drawnow
		end
end		
 
if D==2
  % 5. Diversity Analysis
    figure(15)
    hold on
    plot(ps:ps:fitcount,Div,'b-','LineWidth',2)
    xlabel('Function Evaluations','FontSize', 10)
    ylabel('Diversity Measurement', 'FontSize', 10)
	str = sprintf('Diversity Analysis of FN%d',varargin{:});
    title(str);      
    set(gca,'FontSize',10);
    box on
    drawnow

    % 6. Balance Analysis
    figure(16)
    plot(ps:ps:fitcount,((Div/max(Div))*100),...
         ps:ps:fitcount,((abs((Div-max(Div)))/max(Div))*100),'LineWidth',2)
    xlabel('Function Evaluations','FontSize', 10)
    ylabel('Percentage', 'FontSize', 10)
	str = sprintf('Balance Analysis of FN%d',varargin{:});
    title(str);
    legend('Exploration %','Exploitation %','Location','best') 
    set(gca,'FontSize',10);
    box on
    drawnow			
end

 if D==2
			figure(11)
            plot(gbest(1),gbest(2),'ro','MarkerSize',10,'MarkerFaceColor','r');
			% 3. Trajectory Analysis
            figure(12) 
            str = sprintf('Trajectory of Elite FN%d',varargin{:});
            subplot(211)
            hold on
            plot(ps:ps:fitcount,gbest_position(:,1),'r','LineWidth',2);
% % % 			xlabel('Iteration'); ylabel('x_1 position');title(str);
            legend('PSO','CMAC-PSO')
            box on
            subplot(212)
            hold on
            plot(ps:ps:fitcount,gbest_position(:,2),'r','LineWidth',2);
% % % 			xlabel('Iteration'); ylabel('x_2 position');
            legend('PSO','CMAC-PSO')
            box on
            figure(13)
% % %             contour(x1', x2', x3'); 
            hold on; 
            plot(gbest_position(:,1)',gbest_position(:,2)','LineWidth',2,'Color','r');
            plot(gbest_position(end,1),gbest_position(end,2),'rd','MarkerSize',8,'MarkerFaceColor','k');
			% 4. Average Distance Analysis
			figure(14) 
			hold on
			plot(ps:ps:fitcount,DistanceSum,'r-','LineWidth',2);
			legend('PSO','CMAC-PSO')
            drawnow
 end
end