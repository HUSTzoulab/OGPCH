function F=OGPCH
    F.pre_OGPCH_FLAT = @pre_OGPCH_FLAT;
    F.OGPCH_FLAT = @OGPCH_FLAT;
	F.OGPCH_FLAT_sam100 = @OGPCH_FLAT_sam100;
	F.pre_OGPCH_HIE = @pre_OGPCH_HIE;
	F.OGPCH_HIE = @OGPCH_HIE;
	F.OGPCH_HIE_sam100 = @OGPCH_HIE_sam100;
end

function [] = pre_OGPCH_FLAT(dataset,cons)
	load(['./data_files/',dataset,'_data.mat']);
	load(['./data_files/',dataset,'_labels1.mat']);
	load(['./data_files/',dataset,'_names.mat']);
	load(['./data_files/',dataset,'_genes_id.mat']);
	if (strcmp(dataset,'zeisel')) || (strcmp(dataset,'Mouse_E6_75')) || (strcmp(dataset,'Mouse_E9_5'))
		load(['./data_files/',dataset,'_labels2.mat']);
		labels1 = labels2;
	end
	
	if ( strcmp(dataset,'Mouse_E9_5') )
		r1 = 1;
		z2 = histcounts(labels2);
		max_n = 300;
		sam_idx = zeros(1,1);
		for i =1:size(z2,2)
			find_idx = find(labels2==i);
			if z2(i)>300
				rng(r1)
				sam_idx1 = randperm(size(find_idx,1),max_n);
				find_idx = find_idx(sam_idx1);
			end
			sam_idx = [sam_idx ; find_idx];
		end
		sam_idx(1) = [];
		data = data(sam_idx,:);
		labels1 = labels1(sam_idx,:);
		labels2 = labels2(sam_idx,:);
		labels1 = labels2;
		data = full(data);
	end	 	
	labels1 = double(labels1);
	
	tic
	[m,n]=size(data);
	K=max(labels1);
		
	CK=zeros(K,1);
	Y=zeros(K,n);
	for k=1:K
		CK(k)=sum(labels1==k);
		Y(k,:)=mean(data(labels1==k,:));
	end

	Q = zeros(K,n);
	P = ones(1,n);
	for q=1:n
		if sum(data(:,q)>0)<m*0.03
			data(:,q)=0;
			for i=1:K
				Q(i,q)=1;
			end
		else
			for i=1:K
				if sum(data(labels1==i,q)>0)<size(data(labels1==i,q),1)*0.4
					p=1;
				else
					[~,p,~,~] = ttest2(data(labels1==i,q), data(labels1~=i,q));%[h,p,ci,stats] = ttest2(S(C==i,q), S(C==j,q));
					if p==0
						p=1e-323;
					end
				end
				
				Q(i,q)=p;
				if P(q)>p
					P(q) = p;
				end
			end
		end
		P(q)=-log10(P(q));
		P(q)=P(q)/(K*(K-1));
	end

	A3 = zeros(K,n+cons);
	for i =1:K
		for j =1:n
			if Q(i,j)<1e-3
				A3(i,j)=-1;
			end
		end
	end
	
	A4 = zeros(1,n+cons);
	for i =1:n
		A4(i) = sum(A3(:,i));
	end

	lb = [zeros(1,n),-1./zeros(1,cons)];
	ub = [ones(1,n),1./zeros(1,cons)];

	SK=cell(K,1);
	for k = 1:K
		SK{k} = data(labels1==k,:);
	end
	A2=sparse(1,n+cons);
	for h=1:n
		A2(h)=1;
	end
	A11 = zeros(cons,n);
	p = 1;
	for i = 1:cons
		s = randperm(m,1);
		t = randperm(K,1);
		if t == labels1(s)
			t = mod(t+1,K)+1;
		end
		A11(p,:) = abs(data(s,:)-Y(labels1(s),:))-abs(data(s,:)-Y(t,:));
		p = p+1;
	end
    
	A12 = sparse(size(A11,1),size(A11,1));
	for i = 1:size(A11,1)
		A12(i,i) = -1;
	end
	A1 = [A11,A12];
	Aineq = [A1;A2;A3;A4];
	Aineq = full(Aineq);
	bineq = zeros(size(A11,1)+2+K,1);
	for i=size(A11,1)+2:size(A11,1)+1+K
		bineq(i,1)=-1;
	end
	toc

	file=['./pre_matrix/',dataset,'_FLAT.mat'];
	save(file,'-v7.3');
end


function [] = OGPCH_FLAT(dataset,mu)
	file=['./pre_matrix/',dataset,'_FLAT.mat'];
	load(file);
	
	NUM = 9;
	if (strcmp(dataset,'zeisel')) || (strcmp(dataset,'Mouse_E6_75')) || (strcmp(dataset,'Mouse_E9_5'))
		NUM = 20;
	end
	
	for num1=1:NUM
		SGs=num1*5;
		bineq(cons+1)=SGs;
		bineq(cons+2+K)=SGs*10;
		f = [-P(1:n)/n, mu * ones(1, cons)/cons];
		tic
		x=cplexlp(f,Aineq,bineq,[],[],lb,ub);
		toc
		[X,index] = sort(x(1:n),'descend');
		sg(num1,1:SGs) = index(1:SGs);
	end
	sg = sg-1;
	filename=['./markers/',dataset,'/markers_',dataset,'_FLAT_OGPCH.csv'];
	csvwrite(filename ,sg );
end

function [] = OGPCH_FLAT_sam100(dataset,mu,cons)
	load(['./data_files/',dataset,'_data.mat']);
	load(['./data_files/',dataset,'_labels1.mat']);
	load(['./data_files/',dataset,'_names.mat']);
	load(['./data_files/',dataset,'_genes_id.mat']);
	if (strcmp(dataset,'zeisel')) || (strcmp(dataset,'Mouse_E6_75')) || (strcmp(dataset,'Mouse_E9_5')) || (strcmp(dataset,'Mouse_E65'))
		load(['./data_files/',dataset,'_labels2.mat']);
		labels1 = labels2;
	end
	
	if ( strcmp(dataset,'Mouse_E9_5') )
		r1 = 1;
		z2 = histcounts(labels2);
		max_n = 300;
		sam_idx = zeros(1,1);
		for i =1:size(z2,2)
			find_idx = find(labels2==i);
			if z2(i)>300
				rng(r1)
				sam_idx1 = randperm(size(find_idx,1),max_n);
				find_idx = find_idx(sam_idx1);
			end
			sam_idx = [sam_idx ; find_idx];
		end
		sam_idx(1) = [];
		data = data(sam_idx,:);
		labels1 = labels1(sam_idx,:);
		labels2 = labels2(sam_idx,:);
		labels1 = labels2;
		data = full(data);
	end	 	
	labels1 = double(labels1);
	
	tic
	[m,n]=size(data);
	K=max(labels1);
		
	CK=zeros(K,1);
	Y=zeros(K,n);
	for k=1:K
		CK(k)=sum(labels1==k);
		Y(k,:)=mean(data(labels1==k,:));
	end

	Q = zeros(K,n);
	P = ones(1,n);
	for q=1:n
		if sum(data(:,q)>0)<m*0.03
			data(:,q)=0;
			for i=1:K
				Q(i,q)=1;
			end
		else
			for i=1:K
				if sum(data(labels1==i,q)>0)<size(data(labels1==i,q),1)*0.4
					p=1;
				else
					[~,p,~,~] = ttest2(data(labels1==i,q), data(labels1~=i,q));%[h,p,ci,stats] = ttest2(S(C==i,q), S(C==j,q));
					if p==0
						p=1e-323;
					end
				end
				
				Q(i,q)=p;
				if P(q)>p
					P(q) = p;
				end
			end
		end
		P(q)=-log10(P(q));
		P(q)=P(q)/(K*(K-1));
	end

	A3 = zeros(K,n+cons);
	for i =1:K
		for j =1:n
			if Q(i,j)<1e-3
				A3(i,j)=-1;
			end
		end
	end
	
	A4 = zeros(1,n+cons);
	for i =1:n
		A4(i) = sum(A3(:,i));
	end
	
	lb = [zeros(1,n),-1./zeros(1,cons)];
	ub = [ones(1,n),1./zeros(1,cons)];

	SK=cell(K,1);
	for k = 1:K
		SK{k} = data(labels1==k,:);
	end
	A2=sparse(1,n+cons);
	for h=1:n
		A2(h)=1;
	end
	
	for z=1:100
		A11 = zeros(cons,n);
		p = 1;
		for i = 1:cons
			s = randperm(m,1);
			t = randperm(K,1);
			if t == labels1(s)
				t = mod(t+1,K)+1;
			end
			A11(p,:) = abs(data(s,:)-Y(labels1(s),:))-abs(data(s,:)-Y(t,:));
			p = p+1;
		end    
		A12 = sparse(size(A11,1),size(A11,1));
		for i = 1:size(A11,1)
			A12(i,i) = -1;
		end
		A1 = [A11,A12];
		Aineq = [A1;A2;A3;A4];
		Aineq = full(Aineq);
		bineq = zeros(size(A11,1)+1+K,1);
		for i=size(A11,1)+2:size(A11,1)+1+K
			bineq(i,1)=-1;
		end
		
		SGs=25;
		bineq(cons+1)=SGs;
		bineq(cons+2+K)=SGs*10;
		f = [-P(1:n)/n, mu * ones(1, cons)/cons];

		tic
		x=cplexlp(f,Aineq,bineq,[],[],lb,ub);
		toc
		[X,index] = sort(x(1:n),'descend');
		sg(z,1:SGs) = index(1:SGs);		
	end 
		sg = sg-1;
		filename=['./markers/markers_sam/markers_',num2str(cons),'_sam100_',dataset,'_FLAT_OGPCH.csv'];
		csvwrite(filename ,sg );
end

function  [] = pre_OGPCH_HIE(dataset,cons1, cons2)
	load(['./data_files/',dataset,'_data.mat']);
	load(['./data_files/',dataset,'_labels1.mat']);
	load(['./data_files/',dataset,'_names.mat']);
	load(['./data_files/',dataset,'_genes_id.mat']);
	if (strcmp(dataset,'zeisel')) || (strcmp(dataset,'Mouse_E6_75')) || (strcmp(dataset,'Mouse_E9_5'))
		load(['./data_files/',dataset,'_labels2.mat']);
	end
	
	if ( strcmp(dataset,'Mouse_E9_5') )
		r1 = 1;
		z2 = histcounts(labels2);
		max_n = 300;
		sam_idx = zeros(1,1);
		for i =1:size(z2,2)
			find_idx = find(labels2==i);
			if z2(i)>300
				rng(r1)
				sam_idx1 = randperm(size(find_idx,1),max_n);
				find_idx = find_idx(sam_idx1);
			end
			sam_idx = [sam_idx ; find_idx];
		end
		sam_idx(1) = [];
		data = data(sam_idx,:); 
		labels1 = labels1(sam_idx,:);
		labels2 = labels2(sam_idx,:);
		data = full(data);
	end	
	clearvars find_idx max_n sam_idx sam_idx1 z2
	z2 = histcounts(labels2);
	
	tic
	[m,n]=size(data);
	K1=max(labels1);
	K2=max(labels2);
	CK=zeros(K2,1);
	Y=zeros(K2,n);

	for k=1:K2
		Y(k,:)=mean(data(labels2==k,:));
	end

	LABEL=cell(K1,1);
	for i = 1:K1
		LABEL{i} = unique(labels2(labels1==i));
	end

	Y1=zeros(K1,n);
	Z = zeros(1,K1);
	num_class = zeros(1,K1);
	for k=1:K1
		num = 0;
		for j = 1: size(LABEL{k},1)
			Y1(k,:)= Y1(k,:) + sum(data(labels2==LABEL{k}(j),:));
			num = num + z2(LABEL{k}(j));
		end
		Y1(k,:) = Y1(k,:)/num;
		Z(k) = num;
		num_class(k) = size(LABEL{k},1);
	end
	
	k_p = 0;
	for i=1:k
		if num_class(i)>1
			k_p = k_p + num_class(i);
		end
	end

	Q1 = zeros(K1,n);
	Q2 = zeros(k_p,n);
	P = ones(1,n);
	P1 = ones(1,n);
	P2 = ones(1,n);

	for q=1:n
		if sum(data(:,q)>0)<m*0.03
			data(:,q)=0;
			for i=1:K1
				Q1(i,q)=1;
			end
			for i=1:k_p
				Q2(i,q)=1;
			end
		else
			for i=1:K1
				if sum(data(labels1==i,q)>0)<size(data(labels1==i,q),1)*0.4
					p=1;
				else
					[~,p,~,~] = ttest2(data(labels1==i,q),data(labels1~=i,q));
					if p==0
						p=1e-323;
					end
				end
				Q1(i,q)=p;
				if P1(q)>p
					P1(q) = p;
				end
			end
			P1(q)=-log10(P1(q));
			u=1;
			for i=1:K1
				for j=1:size(LABEL{i},1)
					if size(LABEL{i},1)>1
						if sum(data((labels1==i&labels2==LABEL{i}(j)),q)>0)<size(data((labels1==i&labels2==LABEL{i}(j)),q),1)*0.4
							p=1;
						else
							[~,p,~,~] = ttest2(data((labels1==i&labels2==LABEL{i}(j)),q),data((labels1==i&labels2~=LABEL{i}(j)),q));
							if p==0
								p=1e-323;
							end
						end
						Q2(u,q)=p;
						u = u+1;
						if P2(q)>p
							P2(q) = p;
						end
					end
				end
			end
			P2(q)=-log10(P2(q));
		end
		P(q)=(P1(q)+P2(q))/(K2*(K2-1));
	end
	
	Q = [Q1;Q2];
	
	A4 = zeros(K1+k_p,n+cons1+cons2);
	for i =1:K1
		for j =1:n
			if Q1(i,j)<1e-3
				A4(i,j)=-1;
			end
		end
	end

	for i =1:k_p
		for j =1:n
			if Q2(i,j)<1e-3
				A4(i+K1,j)=-1;
			end
		end
	end

	A11  = zeros(cons1,n);
	p = 1;
	for i = 1:cons1
		s = randperm(m,1);
		t = randperm(K1,1);
		if t == labels1(s)
			t = mod(t+1,K1)+1;
		end
		A11(p,:) = abs(data(s,:)-Y1(labels1(s),:))-abs(data(s,:)-Y1(t,:));
		p = p+1;
	end 
	A12 = sparse(cons1,cons1);
	for i = 1:cons1
		A12(i,i) = -1;
	end
	A21 = zeros(cons2,n);
	p = 1;
	for i = 1:cons2
		s = randperm(m,1);
		t = randperm(size(LABEL{labels1(s)},1),1);
		if t == labels2(s)
			t = mod(t+1,size(LABEL{labels1(s)},1))+1;
		end
		A21(p,:) = abs(data(s,:)-Y(labels2(s),:))-abs(data(s,:)-Y(LABEL{labels1(s)}(t),:));
		p = p+1;
	end

	A23 = sparse(cons2,cons2);
	for i = 1:cons2
		A23(i,i) = -1;
	end
	A22 = sparse(cons2,cons1);
	A13 = sparse(cons1,cons2);
	A1 = [A11,A12,A13];
	A2 = [A21,A22,A23];
	A3=[ones(1,n),sparse(1,cons1+cons2)];
	Aineq = [A1;A2;A3;A4];
	Aineq=full(Aineq);
	
	bineq = zeros(cons1+cons2+1+K1+k_p,1);
	for i=cons1+cons2+2:cons1+cons2+1+K1+k_p
		bineq(i,1)=-1;
	end
	
	lb = [zeros(1,n),-1./zeros(1,cons1+cons2)];
	ub = [ones(1,n),1./zeros(1,cons1+cons2)];
	
	clearvars -except dataset Aineq lb ub  m P K1 n A4 k_p genes_id K2 cons1 cons2 bineq Q Q1 Q2 data labels1 labels2 %后面四个之后可以删掉
	
	file=['./pre_matrix/',dataset,'_HIE.mat'];
	save(file,'-v7.3');
	toc
end

function []=OGPCH_HIE(dataset,mu,nu)
	tic
	file=['./pre_matrix/',dataset,'_HIE.mat'];
	load(file);
	toc

	for num1=1:20
		SGs=num1*5;
		bineq(cons1+cons2+1)=SGs;
		f = [(-P(1:n))/n, nu*ones(1,cons1)/cons1,mu * ones(1,cons2)/cons2];
		tic
		x=cplexlp(f,Aineq,bineq,[],[],lb,ub);
		toc
		[X,index] = sort(x(1:n),'descend');
		sg(num1,1:SGs) = index(1:SGs);
	end

	sg = sg-1;
	filename=['./markers/',dataset,'/markers_',dataset,'_HIE_OGPCH.csv'];
	csvwrite(filename ,sg );
end



function  [] = OGPCH_HIE_sam100(dataset,mu,nu,cons1,cons2)
	load(['./data_files/',dataset,'_data.mat']);
	load(['./data_files/',dataset,'_labels1.mat']);
	load(['./data_files/',dataset,'_names.mat']);
	load(['./data_files/',dataset,'_genes_id.mat']);
	if (strcmp(dataset,'zeisel')) || (strcmp(dataset,'Mouse_E6_75')) || (strcmp(dataset,'Mouse_E9_5'))
		load(['./data_files/',dataset,'_labels2.mat']);
	end
	
	if ( strcmp(dataset,'Mouse_E9_5') )
		r1 = 1;
		z2 = histcounts(labels2);
		max_n = 300;
		sam_idx = zeros(1,1);
		for i =1:size(z2,2)
			find_idx = find(labels2==i);
			if z2(i)>300
				rng(r1)
				sam_idx1 = randperm(size(find_idx,1),max_n);
				find_idx = find_idx(sam_idx1);
			end
			sam_idx = [sam_idx ; find_idx];
		end
		sam_idx(1) = [];
		data = data(sam_idx,:);
		labels1 = labels1(sam_idx,:);
		labels2 = labels2(sam_idx,:);
		data = full(data);
	end	
	clearvars find_idx max_n sam_idx sam_idx1 z2
	z2 = histcounts(labels2);
	
	tic
	[m,n]=size(data);
	K1=max(labels1);
	K2=max(labels2);
	CK=zeros(K2,1);
	Y=zeros(K2,n);

	for k=1:K2
			Y(k,:)=mean(data(labels2==k,:));
	end

	LABEL=cell(K1,1);
	for i = 1:K1
		LABEL{i} = unique(labels2(labels1==i));
	end

	Y1=zeros(K1,n);
	Z = zeros(1,K1);
	num_class = zeros(1,K1);
	for k=1:K1
		num = 0;
		for j = 1: size(LABEL{k},1)
			Y1(k,:)= Y1(k,:) + sum(data(labels2==LABEL{k}(j),:));
			num = num + z2(LABEL{k}(j));
		end
		Y1(k,:) = Y1(k,:)/num;
		Z(k) = num;
		num_class(k) = size(LABEL{k},1);
	end
	
	k_p = 0;
	for i=1:k
		if num_class(i)>1
			k_p = k_p + num_class(i);
		end
	end

	Q1 = zeros(K1,n);
	Q2 = zeros(k_p,n);
	P = ones(1,n);
	P1 = ones(1,n);
	P2 = ones(1,n);

	for q=1:n
		if sum(data(:,q)>0)<m*0.03
			data(:,q)=0;
			for i=1:K1
				Q1(i,q)=1;
			end
			for i=1:k_p
				Q2(i,q)=1;
			end
		else
			for i=1:K1
				if sum(data(labels1==i,q)>0)<size(data(labels1==i,q),1)*0.4
					p=1;
				else
					[~,p,~,~] = ttest2(data(labels1==i,q),data(labels1~=i,q));
					if p==0
						p=1e-323;
					end
				end
				Q1(i,q)=p;
				if P1(q)>p
					P1(q) = p;
				end
			end
			P1(q)=-log10(P1(q));
			u=1;
			for i=1:K1
				for j=1:size(LABEL{i},1)
					if size(LABEL{i},1)>1
						if sum(data((labels1==i&labels2==LABEL{i}(j)),q)>0)<size(data((labels1==i&labels2==LABEL{i}(j)),q),1)*0.4
							p=1;
						else
							[~,p,~,~] = ttest2(data((labels1==i&labels2==LABEL{i}(j)),q),data((labels1==i&labels2~=LABEL{i}(j)),q));
							if p==0
								p=1e-323;
							end
						end
						Q2(u,q)=p;
						u = u+1;
						if P2(q)>p
							P2(q) = p;
						end
					end
				end
			end
			P2(q)=-log10(P2(q));
		end
		P(q)=(P1(q)+P2(q))/(K2*(K2-1));
	end
	
	Q = [Q1;Q2];
	
	A4 = zeros(K1+k_p,n+cons1+cons2);
	for i =1:K1
		for j =1:n
			if Q1(i,j)<1e-3
				A4(i,j)=-1;
			end
		end
	end

	for i =1:k_p
		for j =1:n
			if Q2(i,j)<1e-3
				A4(i+K1,j)=-1;
			end
		end
	end
	
	for z = 1:100
		A11  = zeros(cons1,n);
		p = 1;
		for i = 1:cons1
			s = randperm(m,1);
			t = randperm(K1,1);
			if t == labels1(s)
				t = mod(t+1,K1)+1;
			end
			A11(p,:) = abs(data(s,:)-Y1(labels1(s),:))-abs(data(s,:)-Y1(t,:));
			p = p+1;
		end 
		A12 = sparse(cons1,cons1);
		for i = 1:cons1
			A12(i,i) = -1;
		end
		A21 = zeros(cons2,n);
		p = 1;
		for i = 1:cons2
			s = randperm(m,1);
			t = randperm(size(LABEL{labels1(s)},1),1);
			if t == labels2(s)
				t = mod(t+1,size(LABEL{labels1(s)},1))+1;
			end
			A21(p,:) = abs(data(s,:)-Y(labels2(s),:))-abs(data(s,:)-Y(LABEL{labels1(s)}(t),:));
			p = p+1;
		end
		A23 = sparse(cons2,cons2);
		for i = 1:cons2
			A23(i,i) = -1;
		end
		A22 = sparse(cons2,cons1);
		A13 = sparse(cons1,cons2);
		A1 = [A11,A12,A13];
		A2 = [A21,A22,A23];
		A3=[ones(1,n),sparse(1,cons1+cons2)];
		Aineq = [A1;A2;A3;A4];
		Aineq = full(Aineq);
		
		lb = [zeros(1,n),-1./zeros(1,cons1+cons2)];
		ub = [ones(1,n),1./zeros(1,cons1+cons2)];
		SGs=25;
		f = [(-P(1:n))/n, nu*ones(1,cons1)/cons1,mu * ones(1,cons2)/cons2];
		bineq = zeros(cons1+cons2+1+K1+k_p,1);
		for i=cons1+cons2+2:cons1+cons2+1+K1+k_p
			bineq(i,1)=-1;
		end
		bineq(cons1+cons2+1)=SGs;	
		tic
		x=cplexlp(f,Aineq,bineq,[],[],lb,ub);
		toc
			
		[X,index] = sort(x(1:n),'descend');
		sg(z,1:SGs) = index(1:SGs);	
	end

	sg = sg-1;
	filename=['./markers/markers_sam/markers_',num2str(cons1+cons2),'_sam100_',dataset,'_HIE_OGPCH.csv'];
	csvwrite(filename ,sg );
end