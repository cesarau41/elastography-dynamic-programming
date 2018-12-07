clear all
% close all

%inputs
w=.15; %this is the regularization weight
% g=[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15];
% g_prime=[0;2;3;4;6;7;8;9;10;9;10;11;12;13;14];
load 'rf01.mat'
Im1 = RfDataDouble(1:1700,:);
Im1=Im1/max(Im1(:));
% Im1=[1;2;3;4;5;6;7;8;9;10];
%Im1=imresize(Im1,0.5);

load 'rf03.mat'
Im2 = RfDataDouble(1:1700,:);
Im2=Im2/max(Im2(:));
% Im2=[1;1;2;3;4;5;6;7;8;9];
%Im2=imresize(Im2,0.5);

dmax=0;
dmin=-100;
d=[dmin:1:dmax]; %row vector
max_disp_between_rows=1;

%ensure image is bigger than dmax and dmin.
len_disp=size(d,2);
neg_disp=size(find(d<0),2);
positive_disp=size(find(d>0),2);
last_row=size(Im1,1);
num_cols=size(Im1,2);
num_disp=size(d,2); 
disp_array=[1:1:num_disp];

output_2d=zeros(last_row,num_cols);
disp 'calculating displacement map'
%iterate over all columns on the image
for k=1:num_cols
    k=k
    %assign coulumns to handling variables
    g=Im1(:,k);
    g_prime=Im2(:,k);
    
    %initializations
    C = zeros(last_row,num_disp);
    a = NaN(last_row,num_disp);
    M = NaN(last_row,num_disp);

    %calculate S for efficiency
    %S=zeros(size(d,2),size(d,2));
    % for i=1:size(d,2)
    %     for j=1:size(d,2)
    %         S(i,j)=(d(i)-d(j)).^2;
    %     end
    % end

    %Insert NaN for the rows where the index calculation is invalid for the
    %first rows
    for i=1:neg_disp
        C(i,1:neg_disp+1-i)=NaN(1,neg_disp+1-i);
    end
    %last rows
    for i=1:positive_disp
        C(last_row-i+1,num_disp-(positive_disp-i):num_disp)=NaN(1,positive_disp+1-i);
    end
    
    %Caculate first row. It doesn't use the previous di-1
    tmp=delta(1,d(d>=0), g, g_prime);
    C(1,neg_disp+1:num_disp)=tmp;
    
    %calculate for the second row up to neg_disp. has to handle NaN *****
    for i=2:neg_disp
        %Exlcude NaN from Calculations
        nan_idx=find(isnan(C(i,:)));
        valid_curent_disp=setdiff(disp_array, nan_idx);
        first_allowed_disp=min(valid_curent_disp);
        last_allowed_disp=max(valid_curent_disp);
        
        %evaluate current row
        for j=first_allowed_disp:last_allowed_disp
            tmp=NaN(1,num_disp);
            %get the minimum using values fro previous row
            for p=j-max_disp_between_rows:j+max_disp_between_rows
                if(p>=1 && p<=num_disp) %check borders
%                     tmp2=S(d(j), d(p));
                    tmp(p)=C(i-1,p) + w*S(d(j), d(p));
                end
            end
            [a(i,j),M(i,j)] = min(tmp); %M is the index of the d matrix
%             tmp=delta(i,d(j), g, g_prime)';
            C(i,j)=a(i,j)+delta(i,d(j), g, g_prime)';
        end
        %i=i
    end
    
    %from rows neg_disp + 1 to last_row-positive_disp
    for i=neg_disp+1:last_row-positive_disp
        %row index is previous row and column index is current displacement
        tmp=NaN(num_disp,num_disp); %tmp for current row. 
        %evaluate current row
        for j=1:len_disp
            %get the minimum using values from previous row
            for p=j-max_disp_between_rows:j+max_disp_between_rows
                if(p>=1 && p<=num_disp) %check borders
%                     tmp2=S(d(j), d(p));
                    tmp(p,j)=C(i-1,p) + w*S(d(j), d(p));
                end
            end
        end
        [a(i,:),M(i,:)] = min(tmp); %M is the index of the d matrix
%             tmp=delta(i,d(j), g, g_prime)';
        C(i,:) = a(i,:)+ delta(i,d, g, g_prime)';
        %i=i
    end
    
    %from rows last_row - positive_disp +1 to the end
    for i=last_row-positive_disp+1:last_row
        %evaluate current row
        for j=1:len_disp
            tmp=NaN(1,num_disp);
            %get the minimum using values fro previous row
            for p=j-max_disp_between_rows:j+max_disp_between_rows
                if(p>=1 && p<=num_disp) %check borders
%                     tmp2=S(d(j), d(p));
                    tmp(p)=C(i-1,p) + w*S(d(j), d(p));
                end
            end
            [a(i,j),M(i,j)] = min(tmp); %M is the index of the d matrix
%             tmp=delta(i,d(j), g, g_prime)';
            C(i,j)=a(i,j)+delta(i,d(j), g, g_prime)';
        end
        %i=i
    end

    %Traceback with indexes
    len=size(g,1);
    [val,idx]=min(C(len,:));
    D_idx=zeros(len,1);
    D_idx(len)=idx;
    for i=1:len-1
        j=len-i;
        D_idx(j)=M(j+1,D_idx(j+1));
    end

    D=d(D_idx)';
    output_2d(:,k)=D;

end

figure;imagesc(output_2d); colormap gray; colorbar; title('Axial Displacement');
output_2d_median=medfilt2(output_2d);
figure;imagesc(output_2d_median); colormap gray; colorbar; title('Axial Displacement Median filtered');

disp 'calculating strain. wait a min.'
output_2d_s=strain(output_2d,43);
figure; imagesc(output_2d_s); colormap hot; colorbar;title('Axial Strain');


%*** Anonymous functions ***
function s=S(scalar, d)
%d can be a vector or scalar
s = (scalar-d).^2;
end

function del = delta(i, d, g, g_prime)
%i is a row number and d is displacement column vector
del = abs(g(i)-g_prime(i+d));
end