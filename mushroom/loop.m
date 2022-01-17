%% Loops over a range of pole heights to extract force vs. displacement curve

%%

% finalSol = solution for final pole position
%   x = Sol(1,:), y = Sol(2,:), psi = Sol(3,:), h = Sol(4,:), l = Sol(5,:)
% FvsYp = matrix of values for force vs. displacement plot
%   yp =FvsYp(1,:), f = FvsYp(2,:)
% alpha = dimensionless patch area
% meshPt = number of mesh points, 1000 typically sufficient
% yprng = vector of pole postion values to solve for, i.e. yprng = 0:0.1:20
function [FvsZp, coatPullSol, compZpRng] = loop(alpha, lambda, acoat, rF, rIn, zpRng, f0, k0, gamma, fn, R0, initSol,C0,acoat2)
mesh=(0:0.00005:1).^4;

t=alpha*mesh;   % area mesh points

if isempty(initSol)
    initSol = Init(alpha, mesh, lambda, k0, R0);        % initial guess
end

FvsZp = zeros(2, length(zpRng));    % initialize FvsYp matrix
PullSol = zeros(6, length(mesh), length(zpRng));   % initialize solution matrix
compZpRng = zpRng;

% display a status bar for the calculation
h = waitbar(0,sprintf('Calculating... z_p = %0.1f nm/%0.1f nm', zpRng(1), zpRng(end)));

figure; % open a figure for the intermediate solutions

% loop over the yprng vector
for ii = 1:length(zpRng)
    
    lastwarn('');
   
    % update the status bar
    waitbar(ii/length(zpRng), h, sprintf('Calculating... z_p = %0.2f nm/%0.2f nm', zpRng(ii), zpRng(end)))
    

    try
    
    % solve for the iith value of yprng
    [~,Sol,f] = PullMembrane(alpha, mesh, lambda, acoat(ii), rF(ii), rIn, zpRng(ii), f0, k0, gamma, fn(ii), R0, initSol,C0(ii),acoat2);
    
    [warnStr, warnID] = lastwarn;
    
    if strcmp(warnID, 'MATLAB:bvp4c:RelTolNotMet')
        
        error(warnStr)
        
    end
    
    catch ME
        
        display(ME.message);
        
        PullSol = PullSol(:,:,1:ii-1);
        
        FvsZp = FvsZp(:,1:ii-1);
        
        compZpRng = compZpRng(1:ii-1);
        
        break;
        
    end
    
    % assign iith solution
    PullSol(:,:,ii) = Sol;
    
    % assign iith value of FvsYp
    FvsZp(:,ii) = [Sol(2,1)*R0, f];
    
    % set solution as initial guess for next iteration
    initSol = Sol;
    
end

close(h)    % close status bar

%display(sprintf('Final solution: z_p = %0.2f nm', compZpRng(end)));

