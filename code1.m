close all
clear all
clc
I = imread('E:\�μ�\С������\С������ͼ��\cameraman.bmp');%��ȡԭͼ��
% I=rgb2gray(I);
I=double(I);
% hy = fspecial('sobel');%sobel���� 
% hx = hy'; 
% Iy = imfilter(double(I), hy, 'replicate');%�˲���y�����Ե 
% Ix = imfilter(double(I), hx, 'replicate');%�˲���x�����Ե 
% I= sqrt(Ix.^2 + Iy.^2);%����
% [m,n]=find(I==min(min(I)));
% length(m)
figure
imshow(I,[])
title('ԭʼͼ��')
Iy=I;
%%%%%
%%
% %%%%%%ֱ�ӷ�ˮ��ָ�
% Iw= watershed(I);
% figure
% imshow(Iw,[])
% title('ֱ�ӷ�ˮ��ָ�ͼ��')
% return
%%
% %%%%%%ֱ����ֵ�ָ�
% %��ֵ�ָ�
% If=I;
% [ya1,ya2]=size(If);
% for mm=1:ya1
%     for nn=1:ya2
%         if abs(If(mm,nn))>=100
%             If(mm,nn)=1;
%         else 
%             If(mm,nn)=0;
%         end
%     end
% end
% figure
% imshow(If,[])
% title('ֱ����ֵ�ָ�ͼ��')
% return
%%
%%%%%С���任
[A,B,C,D]=dwt2(I,'sym4');
% [A1,B1,C1,D1]=dwt2(A,'db1');
% [A2,B2,C2,D2]=dwt2(A1,'db1');
%%%%%����
% [C,S]=wavedec2(I,4,'db1');
% for i=1:S(1:2)
% A(:,i)=C((i-1)*S(1:2)+1:i*S(1:2));
% end
% figure
% imshow(A,[])
% Z1=zeros(size(A));
% L = watershed(A);
% [a,b]=find(L==0);
% for j=1:length(a)
%     Z(a(j),b(j))=A(a(j),b(j));
% end
% figure
% imshow(Z,[])
% for i=1:S(1:2)
% C((i-1)*S(1:2)+1:i*S(1:2))=Z1(:,i);
% end
% Iz=waverec2(C,S,'db1');
% figure
% imshow(Iz,[])
% return
%------------------------------
%%%%%%��С���任
% pw=prep2D_appe(I,'ghmap');        %Ԥ�˲�
% w=dec2D_pe(pw,'ghm',1);           %ʹ��ghm��С�����任
% [wh,wz]=size(w);
% figure
% imshow(w,[])
% A=w(1:wh/4,1:wz/4);
%%
%----------------------------------------------------------------
%%%%%%ѡ��
zz=double(A);
%zz=double(Iy);  %%%%%ѡ���ԭͼ��ָ�
I=zz;
%%
%�߽��ں�
[Im,In]=size(I);
L = watershed(I);
[p,q]=find(L==0);
for i=1:length(p)
    if p(i)==1
        a(1)=max(max(I))+1;
        a(2)=abs(I(p(i),q(i))-I(p(i)+1,q(i)));
        a(3)=abs(I(p(i),q(i))-I(p(i),q(i)-1));
        a(4)=abs(I(p(i),q(i))-I(p(i),q(i)+1));
    else if p(i)==Im
            a(2)=max(max(I))+1;
            a(1)=abs(I(p(i),q(i))-I(p(i)-1,q(i)));
            a(3)=abs(I(p(i),q(i))-I(p(i),q(i)-1));
            a(4)=abs(I(p(i),q(i))-I(p(i),q(i)+1));
        else if q(i)==1
                a(3)=max(max(I))+1;
                a(1)=abs(I(p(i),q(i))-I(p(i)-1,q(i)));
                a(2)=abs(I(p(i),q(i))-I(p(i)+1,q(i)));
                a(4)=abs(I(p(i),q(i))-I(p(i),q(i)+1));
                else if q(i)==In
                        a(1)=abs(I(p(i),q(i))-I(p(i)-1,q(i)));
                        a(2)=abs(I(p(i),q(i))-I(p(i)+1,q(i)));
                        a(3)=abs(I(p(i),q(i))-I(p(i),q(i)-1));
                        a(4)=max(max(I))+1;
                    else
                  a(1)=abs(I(p(i),q(i))-I(p(i)-1,q(i)));
                  a(2)=abs(I(p(i),q(i))-I(p(i)+1,q(i)));
                  a(3)=abs(I(p(i),q(i))-I(p(i),q(i)-1));
                  a(4)=abs(I(p(i),q(i))-I(p(i),q(i)+1));
                 
                    end
            end
        end
    end
  while L(p(i),q(i))==0
     ba=min(min(max(a(1),0),max(a(2),0)),min(max(a(3),0),max(a(4),0)));
    k=find(a==ba);
    k=min(k);
   if  k==1;
    L(p(i),q(i))=L(p(i)-1,q(i));
      else if k==2;
           L(p(i),q(i))=L(p(i)+1,q(i));
               else if k==3;
               L(p(i),q(i))=L(p(i),q(i)-1);
                 else if k==4;
               L(p(i),q(i))=L(p(i),q(i)+1);
                     end
                   end  
          end
   end
  if L(p(i),q(i))==0
      a(k)=max(max(I))+1;
  end 
  end

end
%%
%�����ں�             
for i=1:max(max(L))
    [m{i},n{i}]=find(L==i);
    for j=1:length(m{i})
       Ss{i}(j)=I(m{i}(j),n{i}(j));  
    end
     a(i)=sum(Ss{i})/length(m{i});
     b(i)=sum((Ss{i}-a(i)*ones(1,length(m{i}))).^2)/length(m{i});
     c(i)=sum(abs((Ss{i}-a(i)*ones(1,length(m{i}))).^3))/length(m{i});
     d(i)=max(Ss{i});
     e(i)=min(Ss{i});
end
%%
%FCM
data=[a' d' e'];
[center, U, obj_fcm] = FCMClust(data, 2);
maxU=max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2,:) == maxU);
% index3 = find(U(3,:) == maxU);
%%
%��Ƶ��ȡ
for i=1:length(index1)
    for j=1:length(m{index1(i)})
    I(m{index1(i)}(j),n{index1(i)}(j))=0;
    end
end
for i=1:length(index2)
    for j=1:length(m{index2(i)})
    I(m{index2(i)}(j),n{index2(i)}(j))=1;
    end
end
%%
%ƽ��
for i=1:2
    [h{i},g{i}]=find(I==(i-1));
    for j=1:length(h{i})
       Zz{i}(j)=zz(h{i}(j),g{i}(j));  
    end
      Zs(i)=sum(Zz{i})/length(h{i});
    for j=1:length(h{i})
       I(h{i}(j),g{i}(j))=Zs(i);  
    end   
end
%----------------------------------------------------------------
%�۳�3������
% for i=1:length(index3)
%     for j=1:length(m{index3(i)})
%     I(m{index3(i)}(j),n{index3(i)}(j))=2;
%     end
% end
%%
%%%%С����任
% A1=idwt2(A2,B2,C2,D2,'db1');
% A=idwt2(A1,B1,C1,D1,'db1');
If=idwt2(I,B,C,D,'sym4');
%------------
%%%%��С����任
% w(1:wh/4,1:wz/4)=I;
% hw=rec2D_pe(w,'ghm',1);           %���˲�
% If=postp2D_appe(hw,'ghmap');      %��С����任
% figure
% imshow(If,[])
%%
%��ֵ�ָ�
[ya1,ya2]=size(If);
for mm=1:ya1
    for nn=1:ya2
        if abs(If(mm,nn))>=100
            If(mm,nn)=1;
        else 
            If(mm,nn)=0;
        end
    end
end
% If = watershed(If);
% If=waverec2(C,S,'db1');
%%
figure
imshow(I,[])
title('��С���任��Ƶ�ķָ�ͼ')
figure
imshow(If,[])
title('���ڶ�С���任�ķָ�ͼ')

