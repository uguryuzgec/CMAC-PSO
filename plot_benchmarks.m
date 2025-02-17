clear 
close all
clc
% mex cec14_func.cpp -DWINDOWS
% 17-22 hibrid cizilmiyor ve 29-30...
% 1-3 unimodal func.
% 4-16 multimodal func.
% 23-28 composition func.
Xmin=-100;
Xmax=100;
D=2; % boyut sayisi
func_num=[2 8 13 25]; % fonk sayisi
fhd=str2func('cec14_func');

VRmin=repmat(Xmin,1,D);
VRmax=repmat(Xmax,1,D);
% initialize
        x11 = linspace(VRmin(1), VRmax(1), 101);
        x21 = linspace(VRmin(2), VRmax(2), 101);
        x3 = zeros(length(x11), length(x21));
iter = 1;
for m = func_num
    for i = 1:length(x11)
         for j = 1:length(x21)
            x3(i, j) = feval(fhd,[x11(i);x21(j)],m);
         end
    end
    % build simple grid for axes
        [x1, x2] = meshgrid(x11, x21);
        
        % draw surface
        subplot(2,4,iter);
        meshc(x1', x2', x3);
        hold on
        surf(x1', x2', x3, 'Linestyle', 'none'); 
       
        axis tight
        grid on
                    
        % tiles, labels, legend
        xlabel('x_1'), ylabel('x_2'), zlabel('F(x_1, x_2)')   
        view(37, 20)
        subplot(2,4,iter+4)
        contour(x1', x2', x3); 
        iter=iter+1;
end