%% Compute TLB between pairs of migration networks
% Run in worldbankdata directory

% uniform mass
mA = (1/225)*ones(225,1);

workdir = pwd;
resdir = [pwd,'/res/'];
mkdir(resdir);

addpath('../../tools')

dire=dir('mig_mtl*.mat');
numFiles=length(dire);
p=nchoosek(1:numFiles,2);

for i=1:length(p) %edit this to 1:length(p)

    x=p(i,1);
    y=p(i,2);
    
    ni = sprintf('res_%02d_%02d.mat',x,y);
    fi = [resdir, ni];
    
    if ~exist(fi,'file')
        A = load(dire(x).name);
        B = load(dire(y).name);
        
        A = A.A;
        B = B.A;
        
        [tlb,tlb_out,tlb_in,gamma_out,gamma_in] = emd2RTLB(A,B,mA,mA);
        res = struct('x',x,'y',y,'tlb',tlb,'tlb_out',tlb_out,...
            'tlb_in',tlb_in,'gamma_out',gamma_out,'gamma_in',...
            gamma_in,'name_x',dire(x).name,'name_y',dire(y).name);
 
    writemyfile(fi,res);
    end
end

%% Place the pairwise distances into matrix form

cd(resdir);
n = 10; %//Change this depending on size of database
tlb_out_mat = zeros(n);
tlb_in_mat = zeros(n);
tlb_mat = zeros(n);

res = dir('res_*.mat');
nres = length(res);

for k=1:nres
    rk = load(res(k).name);
rk = rk.var;
    ik = rk.x;
    jk = rk.y;
    tlb_out = rk.tlb_out;
    tlb_in = rk.tlb_in;
    tlb = rk.tlb;
    
    tlb_out_mat(ik,jk)=tlb_out;
    tlb_in_mat(ik,jk)=tlb_in;
    tlb_mat(ik,jk)=tlb;
    %if rem(k,)
end

tlb_out_mat = max(tlb_out_mat,tlb_out_mat');
tlb_in_mat = max(tlb_in_mat,tlb_in_mat');
tlb_mat = max(tlb_mat,tlb_mat');

save('output_tlb.mat','tlb_mat')
save('output_tlb_out.mat','tlb_out_mat')
save('output_tlb_in.mat','tlb_in_mat')

%% create labels
cd(workdir)

L = cell(n,1);
for ii=1:length(dire)
    L{ii} = [dire(ii).name(9),'-',dire(ii).name(end-7:end-4)];
end
save('labels-data','L');

%% Plot dendrogram with labels

load([resdir,'output_tlb.mat'])

figure
M=squareform(tlb_mat);
Y=linkage(M);
dendrogram(Y,0,'Orientation','right','Labels',L);


%Z = cmdscale(b0mat,2);

%% Analysis
% inside res folder
load('res_01_06.mat')

top_in_vals = var.gamma_in;


