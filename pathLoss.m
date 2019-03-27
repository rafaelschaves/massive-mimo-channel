function [L_u,varargout] = pathLoss(cell_parameters)

city   = upper(cell_parameters.city);
H_bs   = cell_parameters.bs_height;
H_user = cell_parameters.user_height;
d      = cell_parameters.user_distance; 
f      = cell_parameters.frequency;

% COST-231 Hata Model

switch city
   
    case 'LARGE'
        if((150 <= f) && (f <= 200))
            C_h = 8.29*(log10(1.54*H_user))^2 - 1.1;
        elseif((200 < f) && (f <= 2000))
            C_h = 3.2*(log10(11.75*H_user))^2 - 4.97;
        else
        end
        C = 3;
    case 'MEDIUM'        
        C_h = 0.8 + (1.1*log10(f) - 0.7)*H_user - 1.56*log10(f);
        C = 0;
    case 'SMALL'
        C_h = 0.8 + (1.1*log10(f) - 0.7)*H_user - 1.56*log10(f);
        C = 0;
    otherwise
    
end

L_u_dB = 46.30 + 33.90*log10(f) - 13.82*log10(H_bs) - C_h + ...
         (44.9 - 6.55*log10(H_bs))*log10(d) + C;
     
varargout{1} = L_u_dB;

L_u = 10^(L_u_dB/10);

end