% *************************************************************************
% v3 some speed issues and strain image calculation integration
% v4 includes lateral displacement instead of only axial.
% v5 resolves speed issues
% v6 speed issues and g and g_prime vriables reassignment
% v7 restriction of search space of i-1 to [-1 0 1] + d -> problems, check
% output
% *************************************************************************

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


%has to include zero!!!!!!!
da_max=0;
da_min=-100;
da=[da_min:1:da_max]; %axial row vector
dl_max=4;
dl_min=-4;
dl=[dl_min:1:dl_max];%lateral row vector

% axial displacement variables
len_disp_a=size(da,2);
neg_disp_a=size(find(da<0),2);
positive_disp_a=size(find(da>0),2);
zero_da_idx=find(da==0);

% lateral displacement variables
len_disp_l=size(dl,2);
neg_disp_l=size(find(dl<0),2);
positive_disp_l=size(find(dl>0),2);
zero_dl_idx=find(dl==0);


% ensure image is bigger than dmax and dmin.
last_row=size(Im1,1);
num_cols=size(Im1,2);
num_disp_a=size(da,2); 
num_disp_l=size(dl,2);
disp_array_a=[1:1:num_disp_a];
disp_array_l=[1:1:num_disp_l];

output_2d_a=zeros(last_row,num_cols);
output_2d_l=zeros(last_row,num_cols);
disp 'calculating displacement map'

%initializations
C = NaN(last_row, 2, num_disp_a, num_disp_l);
C(:,1,:,:)=0;
a = NaN(last_row, num_disp_a, num_disp_l);
M = NaN(last_row, num_disp_a, num_disp_l,2); %1 is for del_a, 2 is for del_l
disp_candidates=[repmat([1:num_disp_a]',num_disp_l,1) reshape(repmat([1:num_disp_l]',1,num_disp_a)',num_disp_l*num_disp_a,1)];


%for all columns:
first_valid_dl_idx=zero_dl_idx+1; %+1 because first column is also calculated.
last_valid_dl_idx=num_disp_l;

tic
%iterate over all kth columns on the image
for k=1:num_cols
    toc
    tic
    k=k
    
    %control lateral displacement
    if(k<=neg_disp_l+1) %begining
        old_first_valid_dl_idx=first_valid_dl_idx;
        first_valid_dl_idx=first_valid_dl_idx-1;
        %assign g and gprime to the valid rows:
        g=[ NaN(last_row,neg_disp_l-k+1) Im1(:,dl(first_valid_dl_idx)+k:dl(last_valid_dl_idx)+k)];
        g_prime=[ NaN(last_row,neg_disp_l-k+1) Im2(:,dl(first_valid_dl_idx)+k:dl(last_valid_dl_idx)+k)]; 
    elseif(k>num_cols-positive_disp_l)
        old_last_valid_dl_idx=last_valid_dl_idx;
        last_valid_dl_idx=last_valid_dl_idx-1;
        %assign g and gprime to the valid rows:
        g=[Im1(:,dl(first_valid_dl_idx)+k:dl(last_valid_dl_idx)+k) NaN(last_row,num_disp_l-last_valid_dl_idx)];
        g_prime=[Im2(:,dl(first_valid_dl_idx)+k:dl(last_valid_dl_idx)+k) NaN(last_row,num_disp_l-last_valid_dl_idx)]; 
    else %no restrictions
        %assign g and gprime to the valid rows:
        g=Im1(:,dl(first_valid_dl_idx)+k:dl(last_valid_dl_idx)+k);
        g_prime=Im2(:,dl(first_valid_dl_idx)+k:dl(last_valid_dl_idx)+k); 
    end
    
    %Caculate first row. It doesn't use the previous di-1
    tmp=delta_2d(1,zero_dl_idx,da(da>=0), dl(first_valid_dl_idx:last_valid_dl_idx), g, g_prime);
    calc_range_a=[neg_disp_a+1:num_disp_a];
    calc_range_l=[first_valid_dl_idx:last_valid_dl_idx];
    C(1,2,calc_range_a, calc_range_l)=tmp;
    
    %for all other rows:
    first_valid_da_idx=zero_da_idx;
    last_valid_da_idx=num_disp_a;
    
    
    %for each row.
    for i=2:last_row
        %lateral displacement restrictions are already accounted for,
        %outside the loop.
        %if the row have restrictions in the axial displacement evaluation
        if(i<=neg_disp_a+1)
            old_first_valid_da_idx=first_valid_da_idx;
            first_valid_da_idx=first_valid_da_idx-1;
        elseif(i>last_row-positive_disp_a)
            old_last_valid_da_idx=last_valid_da_idx;
            last_valid_da_idx=last_valid_da_idx-1;
        end
%         tic
%         i=i

        %calculate cost for each axial and lateral displacement of this row/column
        for r=first_valid_dl_idx:last_valid_dl_idx %lateral displ.
            for j=first_valid_da_idx:last_valid_da_idx %axial displ.
                tmp=NaN(num_disp_a,num_disp_l);
                axial_cand=[-1:1:1]'+j;
                axial_cand(axial_cand>num_disp_a)=[];
                axial_cand(axial_cand<1)=[];
                lateral_cand=[-1:1:1]'+r;
                lateral_cand(lateral_cand>num_disp_l)=[];
                lateral_cand(lateral_cand<1)=[];
                num_axial=size(axial_cand,1);
                num_lateral=size(lateral_cand,1);
                cands=[repmat(axial_cand,num_lateral,1) reshape(repmat(lateral_cand,1,num_axial)',num_axial*num_lateral,1)];
%                 cands=disp_candidates(ismember(disp_candidates(:,1),axial_cand),:);
%                 cands=cands(ismember(cands(:,2),lateral_cand),:);
%                 found_axial=unique(cands(:,1));
%                 found_lateral=unique(cands(:,2));
                tmp2=S_2d(da(j), dl(r), da(cands(:,1)), dl(cands(:,2)))'; %all values of da
%                 tmp2=reshape(S_2d_idx(j, r, disp_array_a, disp_array_l),num_disp_a*num_disp_l,1);
                tmp=reshape((C(i,1,axial_cand,lateral_cand)+C(i-1,2,axial_cand,lateral_cand))/2, size(axial_cand,1)*size(lateral_cand,1), 1) + w*tmp2;

                [a(i,j,r),idx]=min(tmp(:));
                
                del_a=cands(idx,1);
                del_l=cands(idx,2);
                M(i,j,r,:)= [del_a,del_l];
                C(i,2,j,r)=delta_2d(i,zero_dl_idx,da(j),dl(r),g,g_prime)+a(i,j,r);

            end
        end
%         toc
        
    end

    %Traceback with indexes
    disp_m_of_row_col=reshape(C(last_row,2,:,:),num_disp_a,num_disp_l);
    [row_idx,col_dx]=find(disp_m_of_row_col==min(disp_m_of_row_col(:)),1,'first');
    D_idx=zeros(last_row,2); %1st column in axial and second is lateral
    D_idx(last_row,:)=[row_idx,col_dx];
    for i=1:last_row-1
        j=last_row-i;
        D_idx(j,:)=M(j+1,D_idx(j+1,1), D_idx(j+1,2), :);
    end

    Da=da(D_idx(:,1))';
    Dl=dl(D_idx(:,2))';
    output_2d_a(:,k)=Da;
    output_2d_l(:,k)=Dl;
    
    %shift cost function
    C(:,1,:,:)=C(:,2,:,:);
    C(:,2,:,:)=NaN;
end

figure;imagesc(output_2d_a); colormap gray; colorbar; title('Axial Displacement');
output_2d_a_median=medfilt2(output_2d_a);
figure;imagesc(output_2d_a_median); colormap gray; colorbar; title('Axial Displacement Median filtered');

figure;imagesc(output_2d_l); colormap gray; colorbar; title('Lateral Displacement');

disp 'calculating strain. wait a min.'
output_2d_s=strain(output_2d_a,43);
figure; imagesc(output_2d_s); colormap hot; colorbar;title('Axial Strain');



%*** Anonymous functions ***
% function s=S(scalar, d)
% %d can be a vector or scalar
% s = (scalar-d).^2;
% end
% 
% function del = delta(i, d, g, g_prime)
% %i is a row number and d is displacement column vector
% del = abs(g(i)-g_prime(i+d));
% end


%****** 2D Versions **********
function del_2d=delta_2d(i,j,da,dl,g,g_prime)
% i is the row number and j is the A-lines and da and dl are axial and
% lateral displacements
del_2d=abs(g(i,j)-g_prime(i+da,j+dl));
end

function s_2d=S_2d(da_i,dl_i,da_prev,dl_prev)
% da_i is the axial displacement at the ith row-pixel and dl_i is the lateral one
% usually, they're a scalar.
% da_prev and dl_prev are the axial and lateral displacement of the previous row
% let's say we want to calculate 10th row so we need value at the 9th
% row;both axial and lateral.
s_2d=(da_i-da_prev).^2+(dl_i-dl_prev).^2;
end

%calculate S for efficiency
%S=zeros(size(d,2),size(d,2));
% for i=1:size(d,2)
%     for j=1:size(d,2)
%         S(i,j)=(d(i)-d(j)).^2;
%     end
% end
