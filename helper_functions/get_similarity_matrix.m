function [S2] = get_similarity_matrix(S1)

norm_r = sqrt(sum(abs(S1).^2,2)); % same as norm(S1,2,'rows')

idx=find(norm_r==0);
norm_r(idx)=0.01;

S2 = (S1 * S1.') ./ (norm_r * norm_r.');


end

