%--------------------------------------------------------------------------
% reduce_interp_1d_linear.m
% Top-level function that applies a recursive algorithm to find a reduced 
% number of X points that still linearly interpolates a 1D dataset within 
% an absolute error tolerance
%--------------------------------------------------------------------------
% [xi, yi] = reduce_interp_1d_linear(X,Y,tol,varargin)
% X : original independent sample points
% Y : original dependent sample points
% tol : absolute error tolerance
% varargin{1}: options structure (see below)
% xi : optimized list of X dimension sample points
% yi : optimized list of Y dimension sample points
%--------------------------------------------------------------------------
% Author: Daniel R. Herber, Graduate Student, University of Illinois at
% Urbana-Champaign
% Date: 08/17/2015
%--------------------------------------------------------------------------
function [xi, yi] = reduce_interp_1d_linear(X,Y,tol,varargin)

    % varargin contains the options structure
    if isempty(varargin)
        opts.interior_optflag = 1; % optimize interior point during spliting algorithm
        opts.post_optflag = 1; % QP for reducign MSE post spliting algorithm
        opts.display_flag = 1; % display (command window and plots)
    else
        opts = varargin{1};
        if ~isfield(opts,'interior_optflag') % is this field present?
            opts.interior_optflag = 1; % on by default
        end
        if ~isfield(opts,'post_optflag') % is this field present?
            opts.post_optflag = 1; % on by default
        end
        if ~isfield(opts,'display_flag') % is this field present?
            opts.display_flag = 1; % on by default
        end
    end

    %%% display original size and start reducing the data set
    if opts.display_flag
        disp(['Original Data Length: ',int2str(length(X))])
        disp('   Finding reduced number of points maintaning error tolerance...')
    end
    
    % call the recursive algorithm to find the reduced data set
    [xi, yi] = splitcheck_1d_linear(X,Y,tol,[],[],opts.interior_optflag);
    
    %%% display finished
    if opts.display_flag
        disp('   ...Finished!')
    end
    
    % the recursive algorithm contains duplicate xi points, remove them
    [xi,IA,~] = unique(xi); % extract only the unique xi points
    yi = yi(IA); % remove the corresponding yi values
    
    %%% display some plots
    if opts.display_flag
        plotflag = 1; % used in the script
        PLOTS_reduce_interp_1d_linear % script
    end
    
    % optionally find the optimal yi with respect to mean squared error at
    % fixed xi locations
    if opts.post_optflag
        
        % find the Hessian and gradient of the mean squared error objective 
        % function through complex step differentiation
        % note that the Hessian is tridiagonal with linear interpolation
        % also this could be found analytically (maybe in the future...)
        % uses a modified version of the following submission
        % http://www.mathworks.com/matlabcentral/fileexchange/18177
        msefun = @(xf,y,X,Y) mean( (Y - interp1(xf, y, X, 'linear')).^2 );
        [H,f] = hessiancsd_mod(@(y) msefun(xi,y,X,Y),zeros(size(yi')));
        
        % find the Jacobian of the linear maximum allowable error
        % constraint
        % note that A is mostly sparse
        % also this could be found analytically (maybe in the future...)
        % uses a modified version of the following submission
        % http://www.mathworks.com/matlabcentral/fileexchange/18176
        errorfun = @(xf,y,X) interp1(xf, y, X, 'linear');
        [A,~] = jacobiancsd_mod(@(y) errorfun(xi,y,X),yi');
        
        % create the absolute maximum error constraint
        A = [A;-A];
        b = [Y(:)+tol;-Y(:)+tol];
        
        % simple bounds
        lb = yi-2*tol; % allowable lower bound
        ub = yi+2*tol; % allowable upper bound
        
        % quadprog options (lowest tolerances for best yi points)
        options = optimoptions('quadprog','Display','none','TolFun',eps,...
            'TolX',eps,'MaxIter',150);
        
        %%% display start of QP problem solving
        if opts.display_flag
        disp('   Solving approximate QP to minimize MSE...')
        end
        
        % solve the quadratic program (should be convex)
        yi = quadprog(H,f,A,b,[],[],lb,ub,yi,options);
        
        %%% display some plots
        if opts.display_flag
            disp('   ...Finished!')
            plotflag = 2; % used in the script
            PLOTS_reduce_interp_1d_linear % script
        end
    end

end