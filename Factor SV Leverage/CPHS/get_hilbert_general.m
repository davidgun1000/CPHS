function [hil_index]=get_hilbert_general(x,r)

%keyboard

size_var=size(x,2);
for j=1:size_var
    for i=1:r
        store(i,j)=bitget(x(1,j),r+1-i);                       
    end
end
% res_temp=store(:,1);
% for j=2:size_var
%     res_temp=[res_temp,store(:,j)];
% end
res=store(1,:);
for j=2:r
    res=horzcat(res,store(j,:));
end
res=double(res);
[hil_index]=sum(res.*2.^(numel(res)-1:-1:0));