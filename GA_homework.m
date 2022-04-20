
clear all;
close all;
clc;


kromozoms = kromozoms()
ref_kromozom = kromozoms(1,:);
best_value = 10;

iteration = 1;
while(iteration<=2)


fitness = zeros(6,1);

% % % % % FÝTNESS % % % % %
val=0;
for k=1:6
    val_trans = kromozoms(k,:);
    for i=1:9
        for j=1:10-i            
            if(val_trans(i)>val_trans(i+j))
                val=val+1;
            end
        end
    end
    fitness(k)=val;
    val=0;
end


fprintf('\niterasyon %d:\n',iteration)
for i=1:6
  fprintf('fitness %d: %d\n',i,fitness(i))  
end

% % % % % % RATÝO % % % % % % %
ratio = fitness./sum(fitness);
cumratio = ratio;

fprintf('\niterasyon %d:\n',iteration)
for i=1:6
  fprintf('ratio %d: %d\n',i,ratio(i))  
end

% % % % % % CUMULATÝVE RATÝO
for i=2:6
    cumratio(i)= cumratio(i-1)+ratio(i);
end

fprintf('\niterasyon %d:\n',iteration)
for i=1:6
  fprintf('cumratio %d: %d\n',i,cumratio(i))  
end

% % % % % % % OBTAÝN BEST VALUE AND BEST KROMOZOM % % % % % % % 

if(max(fitness)>best_value)
    fprintf('\niterasyon %d:\n',iteration)
    best_value = max(fitness)    
    fit_idx= find(fitness==best_value);
    fprintf('\niterasyon %d:\n',iteration)
    best_resuld = kromozoms(fit_idx,:)
end



% % % % % % NATUREL SELECTÝON % % % % % % % % % 
tem_rs = unifrnd(0 ,1,[6,1]);
tempopulation=kromozoms;
for i=1:6
    tem_idx = find(tem_rs(i)<cumratio,1);    
    tempopulation(i,:)=kromozoms(tem_idx,:);
end

fprintf('\niterasyon :%d\n',iteration)
fprintf('after naturel selection :\n')
 tempopulation


 


% % % % % % BACKCROSSÝNG % % % % % % % 
cross_idx = randperm(6);

for i=1:3
    parent_idx1=cross_idx(2*i-1);
    parent_idx2=cross_idx(2*i);
    
    parent_1=tempopulation(parent_idx1,:);
    parent_2=tempopulation(parent_idx2,:);
    
    cross_rs = unifrnd(0,1);
    if(cross_rs<0.9)
        cross_point = unidrnd(10);
        if(cross_point>6)
            emty1=zeros(1,10);
            emty2=zeros(1,10);
            emty1(cross_point+1:end)=parent_2(cross_point+1:end);
            emty2(cross_point+1:end)=parent_1(cross_point+1:end);
            
            for j=1:cross_point
                for k=1:10
                    for h=1:10
                        if(parent_1(j)==emty1(h))
                            parent_1(j)=ref_kromozom(k);
                        end
                    end
                end
                emty1(j)=parent_1(j);
            end
            parent_1=emty1;
            
            for j=1:cross_point
                for k=1:10
                    for h=1:10
                        if(parent_2(j)==emty2(h))
                            parent_2(j)=ref_kromozom(k);
                        end
                    end
                end
                emty2(j)=parent_2(j);
            end
            parent_2=emty2;
           
        end
        if(cross_point<6)
            emty1=zeros(1,10);
            emty2=zeros(1,10);
            emty1(cross_point+1:cross_point+4)=parent_2(cross_point+1:cross_point+4);
            emty2(cross_point+1:cross_point+4)=parent_1(cross_point+1:cross_point+4);    
            
            for j=1:cross_point
                for k=1:10
                    for h=1:10
                        if(parent_1(j)==emty1(h))
                            parent_1(j)=ref_kromozom(k);                            
                        end
                    end
                end
                emty1(j)=parent_1(j);
            end
            for j=cross_point+5:10
                for k=1:10
                    for h=1:10
                        if(parent_1(j)==emty1(h))
                            parent_1(j)=ref_kromozom(k);                            
                        end
                    end
                end
                emty1(j)=parent_1(j);
            end
            
            parent_1 = emty1;
            
            for j=1:cross_point
                for k=1:10
                    for h=1:10
                        if(parent_2(j)==emty2(h))
                            parent_2(j)=ref_kromozom(k);                            
                        end
                    end
                end
                emty2(j)=parent_2(j);
            end
            for j=cross_point+5:10
                for k=1:10
                    for h=1:10
                        if(parent_2(j)==emty2(h))
                            parent_2(j)=ref_kromozom(k);                            
                        end
                    end
                end
                emty2(j)=parent_2(j);
            end
            parent_2 = emty2;     
                 
        end
        
                 
        tempopulation(parent_idx1,:)= parent_1;
        tempopulation(parent_idx2,:)= parent_2;
                 
    end
end

fprintf('\niterasyon :%d crosover\n',iteration)
fprintf('after crosover :\n')
tempopulation

% % % % % % % MUTATÝON % % % % % % % % % % 

mut_rs = unifrnd(0,1,[6,10]);
for i=1:6
    for j=1:10
        if(mut_rs<0.05)
            rnd=randperm(10,2); 
            mut_parent = tempopulation(i,:);
            mut_trans= mut_parent(rnd(1));
            mut_parent(rnd(1)) = mut_parent(rnd(2)); 
            mut_parent(rnd(2))=mut_trans;
            tempopulation(i,:)= mut_parent;
        end
    end
end

fprintf('\niterasyon :%d mutation\n',iteration)
fprintf('after mutation :\n')
tempopulation;


kromozoms = tempopulation;

fprintf('\niterasyon :%d last_population\n',iteration)

tempopulation

iteration=iteration+1;
end

