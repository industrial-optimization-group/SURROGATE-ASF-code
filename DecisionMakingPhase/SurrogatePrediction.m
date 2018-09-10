function ASFVal=SurrogatePrediction(x,Struc)
switch Struc.SurrogateType
    case 'RBF'
        ASFVal = rbfinterp([x' Struc.ZRef]',Struc.Model);
    case 'SVR'
        ASFVal = predict(Struc.Model,[x' Struc.ZRef]);
    case 'KRIGING'
        [ASFVal, ~] = predictor([x' Struc.ZRef], Struc.Model);
end
end