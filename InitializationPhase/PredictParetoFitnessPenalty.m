function Output=PredictParetoFitnessPenalty(x,Struc)
switch Struc.SurrogateType
    case 'RBF'
        Output = rbfinterp(x,Struc.Model);   
    case 'SVR'
        Output = predict(Struc.Model,x');
    case 'KRIGING'
        [Output, ~] = predictor(x', Struc.Model);
end
end