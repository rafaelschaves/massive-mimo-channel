function [od] = orthogonalityDeficiency(channel_matrix)

od = real(1 - det(channel_matrix'*channel_matrix)/ ...
          prod(vecnorm(channel_matrix).^2));

end

