classdef model_results
    properties
        %% General data
        RMSE % Vector of Root Mean Square Error of states
        sigma % Noise standard deviation       
        y_hat % Model prediction over testing data       
        y_test % Actual data for testing
        u_test
        t_test      
        y_train % Noisy data for training
        u_train
        t_train   
        y_data % Whole simulation state data
        u_data % Whole simulation input data
        t % Time series
        C % Measurement matrix
        
        %% HAVOK specific data
        A % System matrix
        B % Input matrix
        p
        r
        c
        d
        q
        w
        
    end
    
    methods
        function obj = model_results(RMSE, sigma, y_hat, y_test, u_test, t_test, y_train, u_train, t_train, y_data, u_data, t, C)
            if nargin == 13            
                obj.RMSE = RMSE;
                obj.sigma = sigma;
                obj.y_hat = y_hat;
                obj.y_test = y_test;
                obj.u_test = u_test;
                obj.t_test = t_test;
                obj.y_train = y_train;
                obj.u_train = u_train;
                obj.t_train = t_train;
                obj.y_data = y_data;
                obj.u_data = u_data;
                obj.t = t;
                obj.C = C;
            end
        end
    end
end
