function [gbest,gbestval,fitcount,t]=PSO_func(fhd,Dimension,Particle_Number,Max_Gen,VRmin,VRmax,X_suru,varargin)

rand('state',sum(100*clock));

ps=Particle_Number;
D=Dimension;
me=D*10000/ps;
gbest_position=zeros(D*10000/ps,D);

cc=[2 2];   %acceleration constants

iwt=0.9-(1:me).*(0.5./me);
% iwt=0.5.*ones(1,me);
if length(VRmin)==1
    VRmin=repmat(VRmin,1,D); 
    VRmax=repmat(VRmax,1,D);
end
mv=0.5*(VRmax-VRmin);
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);
Vmin=repmat(-mv,ps,1);
Vmax=-Vmin;

% % % pos=VRmin+(VRmax-VRmin).*rand(ps,D);
pos=X_suru;

e=feval(fhd,pos',varargin{:});

fitcount=ps;
vel=Vmin+2.*Vmax.*rand(ps,D);%initialize the velocity of the particles
pbest=pos;
pbestval=e; %initialize the pbest and the pbest's fitness value
[gbestval,gbestid]=min(pbestval);
gbest=pbest(gbestid,:);%initialize the gbest and the gbest's fitness value
gbestrep=repmat(gbest,ps,1);
gbest_position(1,:)=gbest;
Distance = abs(pbest(1,1)-pbest(:,1)); % Analiz icin eklendi...
DistanceSum(1) = sum(Distance)/(ps-1);
% 2. Search History Analysis
if D==2 % to be done for only 2-dimension...
	refresh = 100;
	figure(11)

% initialize
        x1 = linspace(VRmin(1), VRmax(1), 101);
        x2 = linspace(VRmin(2), VRmax(2), 101);
        x3 = zeros(length(x1), length(x2));
	for i = 1:length(x1)
		for j = 1:length(x2)
			x3(i, j) = feval(fhd,[x1(i);x2(j)],varargin{:});
		end
    end
    
    
% titles, labels, legend
	str = sprintf('Search History of FN%d',varargin{:});
    contour(x1', x2', x3'); 
    hold on;
	plot(pbest(:,1),pbest(:,2),'go','MarkerSize',8,'MarkerFaceColor','b');
	xlabel('x_1'); ylabel('x_2');  title(str);
    drawnow
	% plot(Xbest(1),Xbest(2),'k*','MarkerSize',8);
	figure(13)
    contour(x1', x2', x3'); 
    hold on; 
    plot(gbest(1),gbest(2),'go','MarkerSize',8,'MarkerFaceColor','k');
	str = sprintf('Trajectory of FN%d',varargin{:});
	xlabel('x_1'); ylabel('x_2');  title(str);
      
% % % 	xlabel('x_1'); ylabel('x_2');  title(str);
% % % 	contour(x1', x2', x3'); hold on;
% % % 	plot(pbest(:,1),pbest(:,2),'bs','MarkerSize',8);
% % % 	drawnow
	% plot(gbest(1),gbest(2),'k*','MarkerSize',8);
end
t=[];i=1;
    while fitcount<D*10000 
        aa=cc(1).*rand(ps,D).*(pbest-pos)+cc(2).*rand(ps,D).*(gbestrep-pos);
        vel=iwt(i).*vel+aa;
        vel=(vel>Vmax).*Vmax+(vel<=Vmax).*vel;
        vel=(vel<Vmin).*Vmin+(vel>=Vmin).*vel;
        pos=pos+vel;
        pos=((pos>=VRmin)&(pos<=VRmax)).*pos...
            +(pos<VRmin).*VRmin+(pos>VRmax).*VRmax;
        e=feval(fhd,pos',varargin{:});
        fitcount=fitcount+ps;
        tmp=(pbestval<e);
        temp=repmat(tmp',1,D);
        pbest=temp.*pbest+(1-temp).*pos;
        pbestval=tmp.*pbestval+(1-tmp).*e;%update the pbest
        [gbestval,tmp]=min(pbestval);
        gbest=pbest(tmp,:);
        gbestrep=repmat(gbest,ps,1);%update the gbest
		if fitcount==100*D||fitcount==200*D||fitcount==300*D||fitcount==500*D...
			||fitcount==1000*D||fitcount==2000*D||fitcount==3000*D||fitcount==4000*D||fitcount==5000*D...
			||fitcount==6000*D||fitcount==7000*D||fitcount==8000*D||fitcount==9000*D||fitcount==10000*D
				t=[t;abs(gbestval-varargin{:}*100)];              
        end
		i=i+1;
        gbest_position(i,:)=gbest;
		Distance = abs(pbest(1,1)-pbest(:,1)); % this is for analysis
		DistanceSum(i) = sum(Distance)/(ps-1);
		% 2. Search History Analysis
		if D==2 && (rem(fitcount/ps,refresh) == 0)% to be done for only 2-dimension...
			figure(11)
% 			plot(pbest(:,1),pbest(:,2),'go','MarkerSize',8,'MarkerFaceColor','b')
		%	plot(gbest(1),gbest(2),'k*','MarkerSize',8);
    % % % 	xlabel('x_1'); ylabel('x_2'); title(str);
			drawnow
		end
    end
		if D==2
            figure(11)
%             plot(gbest(1),gbest(2),'bs','MarkerSize',12,'MarkerFaceColor','b')
			% 3. Trajectory Analysis
			figure(12) 
            str = sprintf('Trajectory of Elite FN%d',varargin{:});
            subplot(211)
            hold on
            plot(ps:ps:fitcount,gbest_position(:,1),'g','LineWidth',2);
			xlabel('Iteration'); ylabel('x_1 position');title(str);
            box on
            subplot(212)
            hold on
            plot(ps:ps:fitcount,gbest_position(:,2),'g','LineWidth',2);
			xlabel('Iteration'); ylabel('x_2 position');
            box on
            figure(13)
            hold on; 
            plot(gbest_position(:,1)',gbest_position(:,2)','LineWidth',2,'Color','g');
            plot(gbest_position(end,1),gbest_position(end,2),'g^','MarkerSize',8,'MarkerFaceColor','k');
			% 4. Average Distance Analysis
			figure(14)
			plot(ps:ps:fitcount,DistanceSum,'g-','LineWidth',2);
			xlabel('Function Evaluations');
			ylabel('Average distance');
			str = sprintf('Average Distance of FN%d',varargin{:});
            title(str);
            drawnow
        end
end