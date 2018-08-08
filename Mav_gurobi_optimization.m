function [x, fval, exitflag] = Mav_gurobi_optimization(objective_function, x0, A, b, Aeq, beq, lb, ub, intcon, options)
% build gurobi model    
model.obj = objective_function;
model.A = [sparse(A); sparse(Aeq)];
n = size(model.A, 2);
model.vtype = repmat('C', n, 1);
model.vtype(intcon) = 'I';
model.sense = [repmat('<', size(A,1),1); repmat('=', size(Aeq,1),1)];
model.rhs = full([b(:); beq(:)]); % rhs must be dense 
if ~isempty(x0)
    model.start = x0;
end
if ~isempty(lb)
    model.lb =  lb;
else
    model.lb = -inf(n,1); % default lb for Matlab is -inf
end
if ~isempty(ub)
    model.ub =  ub;
end

% Extract relevant Gurobi parameters from (subset of) options
params = struct();
params.OutputFlag = 0;
% 
% if isfield(options,'Display') || isa(options,'optim.options.SolverOptions')
%     if any(strcmp(options.Display,{'off','none'}))
%         params.OutputFlag = 0;
%     end
% end
% 
% if isfield(options,'MaxTime') || isa(options,'optim.options.SolverOptions')
%     params.TimeLimit = options.MaxTime;
% end
% 
% if isfield(options,'MaxFeasiblePoints') ...
%         || isa(options,'optim.options.SolverOptions')
%     params.SolutionLimit = options.MaxFeasiblePoints;
% end
% 
% if isfield(options,'RelativeGapTolerance') ...
%         || isa(options,'optim.options.SolverOptions')
%     params.MIPGap = options.RelativeGapTolerance;
% end
% 
% if isfield(options,'AbsoluteGapTolerance') ...
%         || isa(options,'optim.options.SolverOptions')
%     params.MIPGapAbs = options.AbsoluteGapTolerance;
% end

% Solve model with Gurobi
result = gurobi(model, params);

% Resolve model if status is INF_OR_UNBD
if strcmp(result.status,'INF_OR_UNBD')
    params.DualReductions = 0;
    warning('Infeasible or unbounded, resolve without dual reductions to determine...');
    result = gurobi(model,params);
end

% Collect results
x = [];
output.message = result.status;
output.relativegap = [];
output.absolutegap = [];
output.numnodes = result.nodecount;
output.constrviolation = [];

if isfield(result,'x')
    x = result.x;
    if nargout > 3
        slack = model.A*x-model.rhs;
        violA = slack(1:size(A,1));
        violAeq = norm(slack((size(A,1)+1):end),inf);
        viollb = model.lb(:)-x;
        violub = 0;
        if isfield(model,'ub')
            violub = x-model.ub(:);
        end
        output.constrviolation = max([0; violA; violAeq; viollb; violub]);
    end
end

fval = [];

if isfield(result,'objval')
    fval = result.objval;
    if nargout > 3 && numel(intcon) > 0
        U = fval;
        L = result.objbound;
        output.relativegap = 100*(U-L)/(abs(U)+1);
        output.absolutegap = U-L;
    end
end

if strcmp(result.status, 'OPTIMAL')
    exitflag = 1;
elseif strcmp(result.status, 'INFEASIBLE') ...
        || strcmp(result.status, 'CUTOFF')
    exitflag = -2;
elseif strcmp(result.status, 'UNBOUNDED')
    exitflag = -3;
elseif isfield(result, 'x')
    exitflag = 2;
else
    exitflag = 0;
end

end

% Local Functions =========================================================

function [f,intcon,A,b,Aeq,beq,lb,ub,x0,options] = probstruct2args(s)
%PROBSTRUCT2ARGS Get problem structure fields ([] is returned when missing)

f = getstructfield(s,'f');
intcon = getstructfield(s,'intcon');
A = getstructfield(s,'Aineq');
b = getstructfield(s,'bineq');
Aeq = getstructfield(s,'Aeq');
beq = getstructfield(s,'beq');
lb = getstructfield(s,'lb');
ub = getstructfield(s,'ub');
x0 = getstructfield(s,'x0');
options = getstructfield(s,'options');
end
function f = getstructfield(s,field)
%GETSTRUCTFIELD Get structure field ([] is returned when missing)

if isfield(s,field)
    f = getfield(s,field);
else
    f = [];
end
end