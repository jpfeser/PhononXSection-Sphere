function TauInv_NP = Get_TauInv_sphere(k,p,MatParams)

eta_NP = MatParams.eta_NP;
if eta_NP==0
    TauInv_NP = zeros(size(k));
    return
end
vs = MatParams.vs;

sigma_avg = GetSigmaSphere(k,p,MatParams);

TauInv_NP = eta_NP*vs(p).*sigma_avg;

end