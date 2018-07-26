function plotrhod

maxNumCompThreads(1);

Tvec=[40 70 100 130 160 190 220 240 300];

figure;hold on;box on;

alpha=0.568231496731503;
echarge=1.6021766208e-19; 
hbar=1.054571800e-34;   

cons=echarge.^2./hbar;

ng1s=[0:0.5:15 16:1:20 22:2:40 44:4:100]; %in units of 10^10 cm^-2.
ng2s=-ng1s;

for k=1:1
    T=Tvec(k);
load(['monolayersigmas_T' num2str(T) 'K.mat'],'n', 'sigmamono');
nplus=n(2:length(n));
nfull1=[-fliplr(nplus) n];
sigmaplus=sigmamono(2:length(sigmamono));
sigmafull1=[fliplr(sigmaplus) sigmamono];
monocond=@(x) interp1(nfull1,sigmafull1,x,'pchip',NaN);

load(['draggrid-T' num2str(T) '.mat'])
sigmaDinterp=@(n1,n2) interp2(nA,nP,sigmaDgrid,n1,n2,'spline',NaN);%change nimp from 5x15^10 to 10x10^10 

sigmad=zeros(1,length(ng1s));
sigma1=zeros(1,length(ng1s));
sigma2=zeros(1,length(ng2s));

for j=1:length(ng1s)
    sigma1(j)=monocond(ng1s(j));
%     sigma2eff(j)=EMTmono(ng2(j),nrms2,T);
    sigma2(j)=monocond(ng1s(j));
    sigmad(j)=sigmaDinterp(ng1s(j),ng2s(j));
end

rhodfun=@(nq) interp1(ng1s,-sigmad./(sigma1.*sigma2-(sigmad.*4.*alpha.^2.*pi).^2)./cons.*4.*alpha.^2.*pi,nq,'pchip',NaN);
    nplt=(0:0.1:60);
    rhod=rhodfun(nplt);
    plot(nplt,rhod,'LineStyle','-','LineWidth',3)
    save(['rhodhomogdata-T' num2str(T) '.mat'],'nplt','rhod');
end

ng1=(0:1:60);
ng2=-ng1;
    
for k=2:length(Tvec)
    T=Tvec(k);
load(['monolayersigmas_T' num2str(T) 'K.mat'],'n', 'sigmamono');
nplus=n(2:length(n));
nfull1=[-fliplr(nplus) n];
sigmaplus=sigmamono(2:length(sigmamono));
sigmafull1=[fliplr(sigmaplus) sigmamono];
monocond=@(x) interp1(nfull1,sigmafull1,x,'pchip',NaN);

load(['draggrid-T' num2str(T) '.mat'])
sigmaDinterp=@(n1,n2) interp2(nA,nP,sigmaDgrid,n1,n2,'spline',NaN);%change nimp from 5x15^10 to 10x10^10 

sigmad=zeros(1,length(ng1));
sigma1=zeros(1,length(ng1));
sigma2=zeros(1,length(ng2));

for j=1:length(ng1)
    sigma1(j)=monocond(ng1(j));
%     sigma2eff(j)=EMTmono(ng2(j),nrms2,T);
    sigma2(j)=monocond(ng1(j));
    sigmad(j)=sigmaDinterp(ng1(j),ng2(j));
end

rhodfun=@(nq) interp1(ng1,-sigmad./(sigma1.*sigma2-(sigmad.*4.*alpha.^2.*pi).^2)./cons.*4.*alpha.^2.*pi,nq,'pchip',NaN);
    nplt=(0:0.1:60);
    rhod=rhodfun(nplt);
    plot(nplt,rhod,'LineStyle','-','LineWidth',3)
    save(['rhodhomogdata-T' num2str(T) '.mat'],'nplt','rhod');
end

legend({'$T=40K$' '$T=70K$' '$T=100K$' '$T=130K$' '$T=100K$' '$T=190K$' '$T=2200K$' '$T=240K$' '$T=300K$'}, 'Interpreter', 'latex','FontSize',20, 'Location','NorthEast','Orientation','Vertical')

print('-dpdf','rhodvsn.pdf')


end