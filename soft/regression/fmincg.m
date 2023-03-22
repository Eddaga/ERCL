function [theta, J_history] = fmincg(costFunction, initial_theta, options, gradFunction)
% Minimize a continuous differentialble multivariate function. Starting point
% is given by "initial_theta". Function value at optimal point returned by "theta"
% Cost function has to return a scalar value, and gradient function a vector value

% Initialize some useful values
theta = initial_theta;
old_theta = theta;
gradient = zeros(size(theta));
old_gradient = gradient;
pk = gradient;

if ~isfield(options, 'GradObj')
    options.GradObj = 'off';
end

if strcmpi(options.GradObj,'on')
    disp('Using gradient method.')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             Gradient Method Start (quasi-newton)           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    maxIter = options.MaxIter;
    fval = costFunction(theta);
    J_history = zeros(maxIter, 1);
    J_history(1) = fval;
    for iter = 2:maxIter
        [fval, gradient] = costFunction(theta);
        J_history(iter) = fval;
        if(iter > 2)
            yk = gradient - old_gradient;
            sk = theta - old_theta;
            pk = pk - (yk'*pk)/(yk'*sk)*sk;
            if(pk'*gradient > 0.0)
                pk = -gradient;
            end
        else
            pk = -gradient;
        end
        if(norm(pk) < 1.0e-5)
            break;
        end
        alpha = fminunc(@(alpha)(costFunction(theta+alpha*pk)), 0, optimset('GradObj','off','MaxIter',1));
        old_theta = theta;
        old_gradient = gradient;
        theta = theta + alpha*pk;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               Gradient Method End                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    disp('Using conjugate gradient method.')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          Conjugate Gradient Method Start                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    maxIter = options.MaxIter;
    fval = costFunction(theta);
    J_history = zeros(maxIter, 1);
    J_history(1) = fval;
    for iter = 2:maxIter
        [fval, gradient] = costFunction(theta);
        J_history(iter) = fval;
        if(norm(gradient) < 1.0e-5)
            break;
        end
        if(iter == 2)
            pk = -gradient;
        else
            bk = (gradient'*gradient - gradient'*old_gradient)/(old_gradient'*old_gradient);
            bk = max(0,bk);
            pk = -gradient + bk*pk;
        end
        alpha = fminunc(@(alpha)(costFunction(theta+alpha*pk)), 0, optimset('GradObj','off','MaxIter',1));
        old_theta = theta;
        old_gradient = gradient;
        theta = theta + alpha*pk;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          Conjugate Gradient Method End                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end