% function [ Iz ] = quzao( I,miu,bata,esp )
close all
clear all
clc
I = imread('E:\�μ�\С������\С������ͼ��\cameraman.bmp');%��ȡԭͼ�� 
subplot(1,3,1)
imshow(I);
title('ԭʼͼƬ')

I1=imnoise(I,'gaussian',0,0.01);
subplot(1,3,2)
imshow(I1);
title('�����˹����ͼ')
%%
%����Ӧ��ֵС��ȥ��
I1=double(I1);
%%%
miu=0.5;       %��������
bata=0;       %��Ӳ����
esp=0.005;     %������������ 
%%%����Ӧ��ֵ��С��ȥ��
% wp=prep2D_appe(I1,'ghmap');
% w=dec2D_pe(wp,'ghm',1);
% [wh,wz]=size(w);
% for ii=1:4
%     for jj=1:4    
%         Z{ii+(jj-1)*4}=w((ii-1)*wh/4+1:ii*wh/4,(jj-1)*wz/4+1:jj*wz/4);
%     end
% end
[D,Z{1},Z{2},Z{3}]=dwt2(I1,'sym4');
%---------------------------------------------------------------
%---------------------------------------------------------------
namu=[1 10 20];
nu=zeros(1,3);
datanamu=ones(1,3);
% Zb{1}=w(1:wh/4,1:wz/4);
for i=1:3
    while  abs(nu(i)-namu(i))>=esp
    nu(i)=namu(i);
    [m,n]=size(Z{i});
    Zb{i}=zeros(m,n);
    [a{i},b{i}]=find(abs(Z{i})>=namu(i));
    for l=1:length(a{i})
     Zb{i}(a{i}(l),b{i}(l)) =Z{i}(a{i}(l),b{i}(l));
    end

    p1=zeros(m,n);
    p2=zeros(m,n);
    P=zeros(m,n);
    for k=1:m
        for j=1:n
            if abs(Z{i}(k,j))>namu(i)
               p2(k,j)=0;
               if Z{i}(k,j)>namu(i)
                   P(k,j)=(-1+1/(2*bata+1))*(Zb{i}(k,j)-Z{i}(k,j));
               else
                  P(k,j)=(1-1/(2*bata+1))*(Zb{i}(k,j)-Z{i}(k,j));
               end
            else
                P(k,j)=0;
                p2(k,j)=-2*bata*(Zb{i}(k,j)-Z{i}(k,j))^(2*bata)/(namu(i)^(2*bata+1));
            end
           
        end
    end
datanamu(i)=sum(sum(2*P+2*p2));
namu(i)=namu(i)-miu*datanamu(i);
    end
end
%-------------------------------------------------------------
%-------------------------------------------------------------
% for ii=1:4
%     for jj=1:4    
%         w((ii-1)*wh/4+1:ii*wh/4,(jj-1)*wz/4+1:jj*wz/4)=Zb{ii+(jj-1)*4};
%     end
% end
Iz=idwt2(D,Zb{1},Zb{2},Zb{3},'sym4');
% ww=rec2D_pe(w,'ghm',1);
% Iz=postp2D_appe(ww,'ghmap');
subplot(1,3,3)
imshow(Iz,[]);
title('ȥ���ͼƬ')
%%
I=double(I);
[xh, xw]=size(I);
B=8;                %����һ�������ö��ٶ�����λ
MAX=2^B-1;          %ͼ���ж��ٻҶȼ�
MES1=sum(sum((I-I1).^2))/(xh*xw);     %������
PSNR1=20*log10(MAX/sqrt(MES1)) ;       %ԭʼͼ���ֵ�����
MES2=sum(sum((I-Iz).^2))/(xh*xw);     %������
PSNR2=20*log10(MAX/sqrt(MES2)) ;        %ȥ���ķ�ֵ�����