function Output=PredictSurrogate(x,Surrogate,DataInfo)

switch DataInfo.SurrogateType
    case 'RBF'
        Output = rbfinterp(x',Surrogate)';
    case 'SVR'
        Output = predict(Surrogate,x);
    case 'KRIGING'
        [Output, ~] = predictor(x, Surrogate);
end
end