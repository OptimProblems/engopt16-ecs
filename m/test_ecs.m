% Load ECS problem
problem = ECS_ENGOPT2016 ();
dim = size (problem.domain.bounds,2);
domain = problem.domain.bounds;
% domain = problem.domain;

% Random sampling
n = 1e6;
X = double (stk_sampling_randunif (n, dim, domain));
ind = problem.domain.indicator(X);
sum(ind)/n
% X = X(ind,:);

% Optimization 
% [x,fval] = fmincon(problem.F{1}, X, [], [], [], [], problem.domain.bounds(1,:), problem.domain.bounds(2,:), problem.nonlcon);

% options = gaoptimset(@gamultiobj);
% options.PopulationSize = 500;
% options.Generations = 200;
% options.ParetoFraction = 1;
% [x,fval] = gamultiobj(problem.objfunc, dim, [], [], [], [], problem.domain.bounds(1,:), problem.domain.bounds(2,:), problem.nonlcon, options);


% Simulation
[F, C] = problem.simu (X);
is_nan = any(isnan(F),2);
X_obs = X(~is_nan,:);
F_obs = F(~is_nan,:);
C_obs = C(~is_nan,:);

% Check feasibility
feasible = all(C_obs <= 0, 2);
sum(feasible)/size(X_obs,1)
% 
% % Modeling
% options = bmoo_get_default_options();
% models = buildModels (problem, options, X_obs, F_obs, C_obs);
% 
% % Validation set
% n = 1e3;
% Xcv = double (stk_sampling_randunif (n, dim, domain));
% ind = problem.domain.indicator(Xcv);
% Xcv = Xcv(ind,:);
% [Fcv, Ccv] = problem.simu (Xcv);
% is_nan = any(isnan(Fcv),2);
% Xcv = Xcv(~is_nan,:);
% Fcv = Fcv(~is_nan,:);
% Ccv = Ccv(~is_nan,:);
% 
% % Model prediction
% [objPred, constPred] = evaluate (Xcv, models);
% 
% % Box plot
% for i = 1:problem.nObj
%     figure, hold on
%     objMean = objPred{i}.mean;
%     objStd = sqrt(objPred{i}.var);
%     for t = 1:size(Fcv,1)
%         if (Fcv(t,i) <= objMean(t)+1.96*objStd(t)) && (Fcv(t,i) >= objMean(t)-1.96*objStd(t))
%             plot(Fcv(t,i),objMean(t),'ob', 'markerfacecolor', 'b', 'markersize', 3)
%             plot([Fcv(t,i);Fcv(t,i)], [objMean(t)-1.96*objStd(t), objMean(t)+1.96*objStd(t)],'-b')
%         else
%             plot(Fcv(t,i),objMean(t),'or', 'markerfacecolor', 'r', 'markersize', 3)
%             plot([Fcv(t,i);Fcv(t,i)], [objMean(t)-1.96*objStd(t), objMean(t)+1.96*objStd(t)],'-r')
%         end
%     end
%     plot([min(Fcv(:,i)), max(Fcv(:,i))], [min(Fcv(:,i)), max(Fcv(:,i))], '-k')
%     title(sprintf('f_%d', i))
%     xlabel('observations')
%     ylabel('predictions')
%     hold off
% end
% for i = 1:problem.nConst
%     figure, hold on
%     constMean = constPred{i}.mean;
%     constStd = sqrt(constPred{i}.var);
%     for t = 1:size(Ccv,1)
%         if (Ccv(t,i) <= constMean(t)+1.96*constStd(t)) && (Ccv(t,i) >= constMean(t)-1.96*constStd(t))
%             plot(Ccv(t,i),constMean(t),'ob', 'markerfacecolor', 'b', 'markersize', 3)
%             plot([Ccv(t,i);Ccv(t,i)], [constMean(t)-1.96*constStd(t), constMean(t)+1.96*constStd(t)],'-b')
%         else
%             plot(Ccv(t,i),constMean(t),'or', 'markerfacecolor', 'r', 'markersize', 3)
%             plot([Ccv(t,i);Ccv(t,i)], [constMean(t)-1.96*constStd(t), constMean(t)+1.96*constStd(t)],'-r')
%         end
%     end
%     plot([min(Ccv(:,i)), max(Ccv(:,i))], [min(Ccv(:,i)), max(Ccv(:,i))], '-k')
%     title(sprintf('c_%d', i))
%     xlabel('observations')
%     ylabel('predictions')
%     hold off
% end
